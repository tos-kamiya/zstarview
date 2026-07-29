from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from pyproj import Transformer

from zstarview.clouddisc.providers._goes_abi import (
    DATA_VAR,
    GRID_VAR,
    load_cmi_with_area,
)
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
        "semi_major_axis": 6378137.0,
        "semi_minor_axis": 6356752.31414,
        "inverse_flattening": 298.2572221,
        "latitude_of_projection_origin": 0.0,
        "longitude_of_projection_origin": -75.0,
        "sweep_angle_axis": "x",
    }


def test_load_cmi_with_area_attaches_sampler_compatible_area(tmp_path: Path) -> None:
    x = np.array([-0.001, -0.0005, 0.0, 0.0005], dtype=np.float64)
    y = np.array([0.001, 0.0005, 0.0, -0.0005], dtype=np.float64)
    values = np.arange(16, dtype=np.float32).reshape(4, 4)
    ds = xr.Dataset(
        data_vars={
            DATA_VAR: xr.DataArray(
                values,
                dims=("y", "x"),
                coords={"y": y, "x": x},
                attrs={"units": "K", "grid_mapping": GRID_VAR},
            ),
            GRID_VAR: xr.DataArray(0, attrs=_projection_attrs()),
        },
        coords={"y": y, "x": x},
    )
    path = tmp_path / "goes_cmi.nc"
    ds.to_netcdf(path, engine="h5netcdf")
    ds.close()

    da = load_cmi_with_area(path)
    assert da.dtype == np.float32
    assert "area" in da.attrs
    area = da.attrs["area"]
    scale = _projection_attrs()["perspective_point_height"]
    assert np.isclose(area.area_extent[0], float(x[0]) * scale)

    to_lonlat = Transformer.from_crs(area.crs, "EPSG:4326", always_xy=True)
    ix, iy = 2, 1
    lon, lat = to_lonlat.transform(
        float(da.x.values[ix]) * scale,
        float(da.y.values[iy]) * scale,
    )
    sampled = build_bt_sampler(da)(
        np.array([[float(lon)]], dtype=np.float64),
        np.array([[float(lat)]], dtype=np.float64),
    )[0, 0]
    assert np.isclose(sampled, da.values[iy, ix], equal_nan=False)


def test_load_cmi_with_area_disables_time_decoding(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    path = tmp_path / "goes_cmi.nc"
    opened: dict[str, object] = {}

    class _DummyDataset:
        variables = {DATA_VAR: object(), GRID_VAR: object()}
        x = xr.DataArray(np.array([-0.001, 0.0], dtype=np.float64))
        y = xr.DataArray(np.array([0.001, 0.0], dtype=np.float64))

        def __enter__(self) -> _DummyDataset:
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

        def __getitem__(self, key: str):
            if key == DATA_VAR:
                return xr.DataArray(
                    np.zeros((2, 2), dtype=np.float32),
                    dims=("y", "x"),
                    coords={
                        "y": np.array([0.001, 0.0], dtype=np.float64),
                        "x": np.array([-0.001, 0.0], dtype=np.float64),
                    },
                )
            if key == GRID_VAR:
                return xr.DataArray(0, attrs=_projection_attrs())
            raise KeyError(key)

    def fake_open_dataset(*args, **kwargs):
        opened["args"] = args
        opened["kwargs"] = kwargs
        return _DummyDataset()

    monkeypatch.setattr(xr, "open_dataset", fake_open_dataset)

    da = load_cmi_with_area(path)

    assert opened["args"] == (path,)
    assert opened["kwargs"]["decode_times"] is False
    assert "area" in da.attrs
