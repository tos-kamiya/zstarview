"""Tests for alt/az grid debug visualisation."""

from __future__ import annotations

import datetime as dt

import numpy as np
import pytest

from zstarview.clouddisc.altaz_debug import save_altaz_grid_debug_image
from zstarview.clouddisc.altaz_grid import CloudAltAzGrid
from zstarview.clouddisc.types import SourceKey


def _make_grid():
    amount = np.zeros((90, 720), dtype=np.float32)
    amount[45, 180] = 1.0
    missing = np.zeros((90, 720), dtype=np.uint8)
    source_key = SourceKey(
        satellite="G19",
        provider="GOES",
        timeslot_utc=dt.datetime.now(dt.timezone.utc),
        sat_priority="AUTO",
    )
    return CloudAltAzGrid(
        amount=amount,
        missing_mask=missing,
        alt_min_deg=0.0,
        alt_max_deg=90.0,
        az_min_deg=0.0,
        az_max_deg=360.0,
        observer_lat=35.0,
        observer_lon=135.0,
        satellite="G19",
        product="CMIPF-C13",
        time_utc=dt.datetime.now(dt.timezone.utc),
        shells_km=(6374.0, 6376.0, 6378.0),
        source_key=source_key,
        coverage_ratio=1.0,
        source_completeness_ratio=0.95,
    )


def test_save_altaz_grid_debug_image_writes_png(tmp_path):
    grid = _make_grid()
    out_path = save_altaz_grid_debug_image(grid, tmp_path, view_center=(45.0, 180.0))
    if out_path is None:
        pytest.skip("matplotlib is not available")
    assert out_path.exists()
    assert out_path.stat().st_size > 0
    assert out_path.suffix == ".png"


def test_save_altaz_grid_debug_image_path_includes_satellite(tmp_path):
    grid = _make_grid()
    out_path = save_altaz_grid_debug_image(grid, tmp_path)
    if out_path is None:
        pytest.skip("matplotlib is not available")
    assert grid.satellite in out_path.name
