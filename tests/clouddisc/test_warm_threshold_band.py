from __future__ import annotations

import numpy as np
import xarray as xr

from zstarview.clouddisc.sampling.estimate_bt_warm_cold import (
    estimate_bt_warm_hybrid,
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


def test_estimate_bt_warm_hybrid_prefers_local_clear_scene() -> None:
    bt_view = np.linspace(296.0, 297.0, 64, dtype=np.float32).reshape(8, 8)
    mask_inside = np.ones_like(bt_view, dtype=bool)
    eq_samples = np.full(12, 310.0, dtype=np.float32)

    bt_warm = estimate_bt_warm_hybrid(bt_view, mask_inside, eq_samples)

    assert 305.0 <= float(bt_warm) <= 308.0
    assert float(bt_warm) < 309.0


def test_estimate_bt_warm_hybrid_falls_back_when_no_samples() -> None:
    bt_view = np.full((4, 4), np.nan, dtype=np.float32)
    mask_inside = np.zeros((4, 4), dtype=bool)
    eq_samples = np.array([], dtype=np.float32)

    bt_warm = estimate_bt_warm_hybrid(bt_view, mask_inside, eq_samples, fallback_bt_warm=287.5)

    assert float(bt_warm) == 287.5
