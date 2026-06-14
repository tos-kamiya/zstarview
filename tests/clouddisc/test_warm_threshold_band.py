from __future__ import annotations

import numpy as np
import xarray as xr

from zstarview.clouddisc.sampling.estimate_bt_warm_cold import (
    estimate_bt_warm_from_equator_band,
)


class _RecordingArea:
    def __init__(self) -> None:
        self.calls: list[tuple[float, float]] = []

    def get_xy_from_lonlat(self, lon: float, lat: float) -> tuple[int, int]:
        self.calls.append((float(lon), float(lat)))
        return (5, 5)


def test_estimate_bt_warm_samples_plus_minus_five_degree_band() -> None:
    area = _RecordingArea()
    da = xr.DataArray(
        np.full((12, 12), 295.0, dtype=np.float32),
        dims=("y", "x"),
        attrs={"area": area},
    )

    bt_warm, samples = estimate_bt_warm_from_equator_band(
        da,
        lon_center_deg=140.0,
        delta_lon=0.0,
        equator_lat=0.0,
        step_deg=1.0,
        half=1,
        warm_p=97.0,
        equator_lat_half_band_deg=5.0,
    )

    assert float(bt_warm) == 295.0
    assert samples.size == 11
    seen_lats = sorted({round(lat, 3) for _lon, lat in area.calls})
    assert seen_lats == [float(value) for value in range(-5, 6)]
