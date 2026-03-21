from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from pyproj import Transformer

from zstarview.clouddisc import CloudDiscConfig
from zstarview.clouddisc.providers._hima_isatss import (
    DATA_VAR,
    GRID_VAR,
    TemplateMeta,
    generate_sparse_layout,
    stitch_tiles_from_paths,
)
from zstarview.clouddisc.providers.hima import HimaProvider
from zstarview.clouddisc.sampling.bt_sampler import build_bt_sampler

pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:You will likely lose important projection information when converting to a PROJ string:UserWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:Conversion of an array with ndim > 0 to a scalar is deprecated:DeprecationWarning"
    ),
]


def _projection_attrs() -> dict[str, float | str]:
    return {
        "grid_mapping_name": "geostationary",
        "perspective_point_height": 35786023.0,
        "latitude_of_projection_origin": 0.0,
        "longitude_of_projection_origin": 140.7,
        "sweep_angle_axis": "y",
        "semi_major_axis": 6378137.0,
        "semi_minor_axis": 6356752.31414,
    }


def _write_sparse_tile_set(tile_dir: Path) -> list[Path]:
    meta = TemplateMeta(
        bucket="local",
        tile_path=tile_dir / "template.nc",
        tile_height=2,
        tile_width=2,
        product_rows=20,
        product_cols=20,
        tile_count=88,
        x0=-0.001,
        y0=0.001,
        x_step=0.0001,
        y_step=-0.0001,
        geos_scale=35786023.0 * 1e-6,
        crs=None,  # not needed for layout generation
    )
    paths: list[Path] = []
    for record in generate_sparse_layout(meta):
        x = meta.x0 + (record.col_offset + np.arange(meta.tile_width, dtype=np.float64)) * meta.x_step
        y = meta.y0 + (record.row_offset + np.arange(meta.tile_height, dtype=np.float64)) * meta.y_step
        value = np.full((meta.tile_height, meta.tile_width), float(int(record.token)), dtype=np.float32)
        ds = xr.Dataset(
            data_vars={
                DATA_VAR: xr.DataArray(
                    value,
                    dims=("y", "x"),
                    coords={"y": y, "x": x},
                    attrs={"standard_name": "brightness_temperature", "units": "kelvin"},
                ),
                GRID_VAR: xr.DataArray(0, attrs=_projection_attrs()),
            },
            coords={"y": y, "x": x},
            attrs={
                "product_rows": meta.product_rows,
                "product_columns": meta.product_cols,
                "number_product_tiles": meta.tile_count,
                "tile_row_offset": record.row_offset,
                "tile_column_offset": record.col_offset,
                "channel_id": 13,
            },
        )
        path = tile_dir / f"OR_HFD-020-B12-M1C13-T{record.token}_TEST.nc"
        ds.to_netcdf(path)
        ds.close()
        paths.append(path)
    return paths


def test_stitch_isatss_sparse_tiles_into_full_field(tmp_path: Path) -> None:
    paths = _write_sparse_tile_set(tmp_path)
    merged = stitch_tiles_from_paths(paths, source_label=str(tmp_path))
    try:
        assert merged[DATA_VAR].shape == (20, 20)
        assert merged.attrs["stitched_from_tiles"] == 88
        assert float(merged[DATA_VAR].attrs["coverage_fraction"]) == 0.88
        assert np.isnan(merged[DATA_VAR].values[0, 0])
        assert np.isnan(merged[DATA_VAR].values[0, 1])
        assert np.isfinite(merged[DATA_VAR].values[4, 4])
        assert int(merged["tile_coverage_mask"].values.sum()) == 352
    finally:
        merged.close()


def test_hima_provider_local_tiles_attach_area_and_sample(tmp_path: Path) -> None:
    tile_dir = tmp_path / "tiles"
    tile_dir.mkdir()
    paths = _write_sparse_tile_set(tile_dir)
    provider = HimaProvider(CloudDiscConfig(cache_dir=tmp_path / "cache"))
    da, used_time, src_paths = provider.fetch_bt_c13_from_local_dir(
        tile_dir,
        used_time=dt.datetime(2026, 3, 21, 6, 10, tzinfo=dt.timezone.utc),
    )

    assert da.shape == (20, 20)
    assert da.dtype == np.float32
    assert len(src_paths) == len(paths) == 88
    assert used_time == dt.datetime(2026, 3, 21, 6, 10, tzinfo=dt.timezone.utc)
    assert "area" in da.attrs

    area = da.attrs["area"]
    geos_scale = _projection_attrs()["perspective_point_height"] * 1e-6
    assert np.isclose(area.area_extent[0], float(da.x.values[0]) * geos_scale)
    to_lonlat = Transformer.from_crs(area.crs, "EPSG:4326", always_xy=True)
    ix, iy = 10, 10
    lon, lat = to_lonlat.transform(
        float(da.x.values[ix]),
        float(da.y.values[iy]),
    )
    sampled = build_bt_sampler(da)(
        np.array([[lon]], dtype=np.float64),
        np.array([[lat]], dtype=np.float64),
    )[0, 0]
    assert np.isclose(sampled, da.values[iy, ix], equal_nan=False)
