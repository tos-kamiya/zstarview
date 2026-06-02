from __future__ import annotations

import datetime as dt
import threading
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from pyproj import CRS, Transformer

from zstarview.clouddisc import CloudDiscConfig
from zstarview.clouddisc.providers import hima as hima_module
from zstarview.clouddisc.providers._hima_isatss import (
    DATA_VAR,
    GRID_VAR,
    TemplateMeta,
    TileRecord,
    generate_sparse_layout,
    load_template_from_tile,
    select_equator_tiles,
    select_needed_tiles,
    stitch_tiles_from_paths,
)
from zstarview.clouddisc.providers.hima import HimaProvider
from zstarview.clouddisc.sampling.bt_sampler import build_bt_sampler
from zstarview.clouddisc.types import DataNotFoundError, DownloadCancelledError

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
        crs=CRS.from_cf(_projection_attrs()),
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
        ds.to_netcdf(path, engine="h5netcdf")
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


def test_hima_provider_stitch_local_paths_attach_area_and_sample(tmp_path: Path) -> None:
    tile_dir = tmp_path / "tiles"
    tile_dir.mkdir()
    paths = _write_sparse_tile_set(tile_dir)
    provider = HimaProvider(CloudDiscConfig(cache_dir=tmp_path / "cache"))
    da = provider._stitch_local_paths(
        paths,
        source_label=str(tile_dir),
        observer_lat=None,
        observer_lon=None,
    )

    try:
        assert da.shape == (20, 20)
        assert da.dtype == np.float32
        assert len(paths) == 88
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
    finally:
        da.close()


def test_equator_tile_selection_extends_partial_himawari_subset(tmp_path: Path) -> None:
    tile_dir = tmp_path / "tiles"
    tile_dir.mkdir()
    paths = _write_sparse_tile_set(tile_dir)
    provider = HimaProvider(CloudDiscConfig(cache_dir=tmp_path / "cache"))
    meta = load_template_from_tile(sorted(tile_dir.glob("*M1C13*.nc"))[0], bucket="local")

    render_tiles, _poly_x, _poly_y = select_needed_tiles(
        lat_deg=35.483,
        lon_deg=140.7,
        meta=meta,
        cloud_shell_km=6376.0,
        azimuth_samples=1440,
        margin_tiles=1,
    )
    equator_tiles, _eq_x, _eq_y = select_equator_tiles(
        lon_center_deg=140.7,
        meta=meta,
        delta_lon=5.0,
        equator_lat=0.0,
        step_deg=1.0,
        margin_tiles=0,
    )
    render_tokens = {record.token for record in render_tiles}
    equator_tokens = {record.token for record in equator_tiles}
    assert equator_tokens
    assert equator_tokens - render_tokens

    selected_paths = [
        path for path in paths if hima_module.extract_tile_token(path.name) in (render_tokens | equator_tokens)
    ]
    da = provider._stitch_local_paths(
        selected_paths,
        source_label=str(tile_dir),
        observer_lat=35.483,
        observer_lon=140.7,
    )
    try:
        assert len(selected_paths) == len(render_tokens | equator_tokens)
        assert len(selected_paths) > len(render_tokens)
        assert da.attrs["source_key_count"] == len(selected_paths)
    finally:
        da.close()


def test_find_isatss_accepts_incomplete_latest_slot(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    tile_dir = tmp_path / "tiles"
    tile_dir.mkdir()
    template_paths = _write_sparse_tile_set(tile_dir)
    provider = HimaProvider(CloudDiscConfig(cache_dir=tmp_path / "cache"))
    latest_time = dt.datetime(2026, 3, 21, 18, 10, tzinfo=dt.timezone.utc)
    stable_time = dt.datetime(2026, 3, 21, 18, 0, tzinfo=dt.timezone.utc)

    incomplete_keys = [f"AHI-L2-FLDK-ISatSS/2026/03/21/1810/OR_HFD-020-B12-M1C13-T{i:03d}_TEST.nc" for i in range(1, 15)]
    stable_keys = [f"AHI-L2-FLDK-ISatSS/2026/03/21/1800/OR_HFD-020-B12-M1C13-T{i:03d}_TEST.nc" for i in range(1, 89)]

    def fake_find_matching_keys(when_utc: dt.datetime, *, satellite: str, product: str, timeout_s=None) -> tuple[str, list[str]]:
        assert satellite == "HIMAWARI"
        assert product == "ISatSS-B13"
        del timeout_s
        rounded = when_utc.astimezone(dt.timezone.utc).replace(second=0, microsecond=0)
        if rounded == latest_time:
            return "noaa-himawari9", incomplete_keys
        if rounded == stable_time:
            return "noaa-himawari9", stable_keys
        raise FileNotFoundError

    monkeypatch.setattr(hima_module, "find_matching_keys", fake_find_matching_keys)
    monkeypatch.setattr(provider, "_download", lambda bucket, key, abort_event=None: template_paths[0])

    bucket, keys, used_time, expected_tile_count = provider._find_isatss(latest_time)

    assert bucket == "noaa-himawari9"
    assert keys == incomplete_keys
    assert used_time == latest_time
    assert expected_tile_count == 88


def test_fetch_bt_c13_stops_when_cancelled(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    provider = HimaProvider(CloudDiscConfig(cache_dir=tmp_path / "cache"))
    abort_event = threading.Event()
    abort_flag = abort_event
    used_time = dt.datetime(2026, 3, 21, 18, 0, tzinfo=dt.timezone.utc)
    keys = [f"tile-{idx:03d}.nc" for idx in range(1, 4)]
    seen: list[str] = []

    monkeypatch.setattr(provider, "_find_isatss", lambda when_utc, abort_event=None: ("noaa-himawari9", keys, used_time))

    def fake_download(bucket: str, key: str, *, abort_event=None):
        assert bucket == "noaa-himawari9"
        assert abort_event is abort_flag
        seen.append(key)
        if len(seen) == 2:
            raise DownloadCancelledError("Cancelled while downloading")
        return tmp_path / key

    monkeypatch.setattr(provider, "_download", fake_download)

    with pytest.raises(DownloadCancelledError):
        provider.fetch_bt_c13(
            when_utc=used_time,
            abort_event=abort_event,
        )

    assert seen == keys[:2]


def test_select_keys_for_observer_rejects_missing_required_tiles(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tile_dir = tmp_path / "tiles"
    tile_dir.mkdir()
    template_paths = _write_sparse_tile_set(tile_dir)
    provider = HimaProvider(CloudDiscConfig(cache_dir=tmp_path / "cache"))
    near_missing_tile = TileRecord(token="901", row_offset=0, col_offset=0, x_min=0.0, x_max=1.0, y_min=0.0, y_max=1.0)
    present_tile = TileRecord(token="002", row_offset=0, col_offset=2, x_min=2.0, x_max=3.0, y_min=0.0, y_max=1.0)

    keys = [
        f"AHI-L2-FLDK-ISatSS/2026/03/21/1800/{template_paths[1].name}",
    ]

    monkeypatch.setattr(
        hima_module,
        "select_needed_tiles",
        lambda **_kwargs: ([near_missing_tile, present_tile], np.array([], dtype=np.float64), np.array([], dtype=np.float64)),
    )
    monkeypatch.setattr(hima_module, "select_equator_tiles", lambda **_kwargs: ([], np.array([], dtype=np.float64), np.array([], dtype=np.float64)))
    monkeypatch.setattr(hima_module, "tile_distance_km", lambda record, meta, *, observer_lat, observer_lon: 10.0 if record.token == "901" else 60.0)
    monkeypatch.setattr(provider, "_download", lambda bucket, key, abort_event=None: template_paths[0])

    with pytest.raises(DataNotFoundError):
        provider._select_keys_for_observer(
            bucket="noaa-himawari9",
            keys=keys,
            when_utc=dt.datetime(2026, 3, 21, 18, 0, tzinfo=dt.timezone.utc),
            observer_lat=0.0,
            observer_lon=0.0,
            cloud_shell_km=6376.0,
            azimuth_samples=1440,
            margin_tiles=1,
        )


def test_stitch_local_paths_fills_far_missing_render_tiles_as_clear_sky(tmp_path: Path) -> None:
    tile_dir = tmp_path / "tiles"
    tile_dir.mkdir()
    paths = _write_sparse_tile_set(tile_dir)
    provider = HimaProvider(CloudDiscConfig(cache_dir=tmp_path / "cache"))
    meta = load_template_from_tile(sorted(tile_dir.glob("*M1C13*.nc"))[0], bucket="local")
    far_missing_tile = generate_sparse_layout(meta)[0]
    selected_paths = [path for path in paths if hima_module.extract_tile_token(path.name) != "001"]

    da = provider._stitch_local_paths(
        selected_paths,
        source_label=str(tile_dir),
        observer_lat=0.0,
        observer_lon=0.0,
        far_missing_render_tiles=[far_missing_tile],
    )
    try:
        row0 = int(far_missing_tile.row_offset)
        col0 = int(far_missing_tile.col_offset)
        row1 = row0 + 2
        col1 = col0 + 2
        assert np.allclose(da.values[row0:row1, col0:col1], 315.0)
        assert float(da.attrs["coverage_fraction"]) == 1.0
        assert da.attrs["far_clear_fill_tile_count"] == 1
    finally:
        da.close()
