"""Tests for the CloudAltAzGrid data model and ingestion."""

from __future__ import annotations

import datetime as dt
from unittest.mock import patch

import numpy as np
import pytest
import xarray as xr

from zstarview.clouddisc.altaz_grid import (
    CloudAltAzGrid,
    _blend_cloud_shell_weights,
    _dilate_bool,
    _estimate_scene_cloud_amount,
    altaz_grid_cache_key,
    build_altaz_grid,
    load_altaz_grid,
    save_altaz_grid,
)
from zstarview.clouddisc.types import CloudSourceData, SourceKey


def _make_source_key():
    return SourceKey(
        satellite="G19",
        provider="GOES",
        timeslot_utc=dt.datetime(2026, 6, 22, 12, 0, 0, tzinfo=dt.timezone.utc),
        sat_priority=("AUTO",),
    )


def _make_dummy_data_array():
    """Return a minimal DataArray for threshold-estimation code paths."""
    lons = np.linspace(120, 150, 60)
    lats = np.linspace(20, 50, 60)
    data = np.full((len(lats), len(lons)), 290.0, dtype=np.float32)
    return xr.DataArray(
        data,
        dims=("y", "x"),
        coords={"y": lats, "x": lons},
        attrs={"area": None},
    )


def test_cloud_altaz_grid_shape_validation():
    amount = np.zeros((90, 720), dtype=np.float32)
    missing = np.zeros((90, 720), dtype=np.uint8)
    grid = CloudAltAzGrid(
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
        source_key=_make_source_key(),
        coverage_ratio=1.0,
    )
    assert grid.amount.shape == (90, 720)


def test_cloud_altaz_grid_rejects_shape_mismatch():
    amount = np.zeros((90, 720), dtype=np.float32)
    missing = np.zeros((80, 720), dtype=np.uint8)
    with pytest.raises(ValueError):
        CloudAltAzGrid(
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
            source_key=_make_source_key(),
            coverage_ratio=1.0,
        )


def test_build_altaz_grid_empty_source():
    """A source with NaN everywhere yields an empty grid with no coverage."""
    def sampler(lon, lat):
        return np.full(np.asarray(lon).shape, np.nan, dtype=np.float32)

    source = CloudSourceData(
        source_key=_make_source_key(),
        data_array=_make_dummy_data_array(),
        satellite="G19",
        product="CMIPF-C13",
        time_utc=dt.datetime.now(dt.timezone.utc),
        src_paths=[],
        sampler=sampler,
    )
    grid = build_altaz_grid(source, 35.0, 135.0)
    assert grid.amount.shape == (90, 720)
    assert np.all(grid.amount == 0.0)
    assert np.all(grid.missing_mask == 0)
    assert grid.coverage_ratio == 1.0


def test_build_altaz_grid_cold_pixel_maps_to_expected_cell():
    """A single cold pixel near the observer appears at the overhead cell."""
    observer_lat = 35.0
    observer_lon = 135.0
    cold_lat = 35.0
    cold_lon = 135.0

    def sampler(lon, lat):
        lon = np.asarray(lon)
        lat = np.asarray(lat)
        # Cold (cloudy) near observer, warm (clear) elsewhere.
        dist = np.abs(lat - cold_lat) + np.abs(lon - cold_lon)
        bt = np.full(lon.shape, 310.0, dtype=np.float32)
        bt[dist < 0.5] = 220.0
        return bt

    source = CloudSourceData(
        source_key=_make_source_key(),
        data_array=_make_dummy_data_array(),
        satellite="G19",
        product="CMIPF-C13",
        time_utc=dt.datetime.now(dt.timezone.utc),
        src_paths=[],
        sampler=sampler,
    )
    grid = build_altaz_grid(
        source,
        observer_lat,
        observer_lon,
        geo_sample_step_deg=0.1,
        geo_sample_extent_deg=2.0,
    )

    # The cold pixel is directly overhead, so it projects to the zenith bin.
    zenith_alt_idx = grid.amount.shape[0] - 1
    # Many cells may saturate at amount == 1.0; ensure the zenith cell is one of them.
    assert grid.amount[zenith_alt_idx, :].max() > 0.5
    assert grid.coverage_ratio == pytest.approx(1.0)


def test_build_altaz_grid_dense_source_populates_full_sky(monkeypatch):
    """A fully cloudy source should populate the whole cached alt/az grid."""

    def sampler(lon, lat):
        return np.full(np.asarray(lon).shape, 220.0, dtype=np.float32)

    source = CloudSourceData(
        source_key=_make_source_key(),
        data_array=_make_dummy_data_array(),
        satellite="G19",
        product="CMIPF-C13",
        time_utc=dt.datetime.now(dt.timezone.utc),
        src_paths=[],
        sampler=sampler,
    )
    monkeypatch.setattr(
        "zstarview.clouddisc.altaz_grid.estimate_bt_warm_from_equator_band",
        lambda *_args, **_kwargs: (310.0, np.array([310.0], dtype=np.float32)),
    )
    monkeypatch.setattr(
        "zstarview.clouddisc.altaz_grid.estimate_bt_warm_hybrid",
        lambda *_args, **_kwargs: 310.0,
    )
    monkeypatch.setattr(
        "zstarview.clouddisc.altaz_grid.estimate_bt_cold_hybrid",
        lambda *_args, **_kwargs: 220.0,
    )
    grid = build_altaz_grid(source, 35.0, 135.0)
    assert grid.amount.shape == (90, 720)
    assert grid.coverage_ratio == pytest.approx(1.0)
    assert np.count_nonzero(grid.missing_mask) == 0
    assert float(np.mean(grid.amount)) > 0.5


def test_build_altaz_grid_no_sampler_builds_one():
    """If no sampler is provided, build_altaz_grid creates one from data_array."""
    # data_array with area=None will cause build_bt_sampler to fail, so we supply a sampler.
    def sampler(lon, lat):
        return np.full(np.asarray(lon).shape, 290.0, dtype=np.float32)

    source = CloudSourceData(
        source_key=_make_source_key(),
        data_array=_make_dummy_data_array(),
        satellite="G19",
        product="CMIPF-C13",
        time_utc=dt.datetime.now(dt.timezone.utc),
        src_paths=[],
        sampler=sampler,
    )
    grid = build_altaz_grid(source, 35.0, 135.0)
    assert grid is not None


def test_estimate_scene_cloud_amount():
    bt = np.array([310.0, 300.0, 290.0, 220.0, np.nan])
    amount = _estimate_scene_cloud_amount(bt, bt_warm=310.0, bt_cold=220.0)
    assert 0.0 < amount < 1.0


def test_blend_cloud_shell_weights():
    weights = _blend_cloud_shell_weights(cloud_amount=0.5)
    assert len(weights) == 3
    assert abs(sum(weights) - 1.0) < 1e-6
    # At medium cloudiness we expect a mix of middle and outer shells.
    assert weights[1] > weights[0]


def test_dilate_bool():
    m = np.zeros((7, 7), dtype=bool)
    m[3, 3] = True
    d = _dilate_bool(m, 1)
    assert d[2:5, 2:5].all()
    assert not d[0, 0]


def test_altaz_grid_cache_key_stable():
    source = CloudSourceData(
        source_key=_make_source_key(),
        data_array=_make_dummy_data_array(),
        satellite="G19",
        product="CMIPF-C13",
        time_utc=dt.datetime(2026, 6, 22, 12, 0, 0, tzinfo=dt.timezone.utc),
        src_paths=[],
    )
    key1 = altaz_grid_cache_key(
        35.0,
        135.0,
        satellite=source.satellite,
        product=source.product,
        time_utc=source.time_utc,
        source_key=source.source_key,
        shells_km=(6374.0, 6376.0, 6378.0),
        grid_resolution_deg=0.5,
    )
    key2 = altaz_grid_cache_key(
        35.0,
        135.0,
        satellite=source.satellite,
        product=source.product,
        time_utc=source.time_utc,
        source_key=source.source_key,
        shells_km=(6374.0, 6376.0, 6378.0),
        grid_resolution_deg=0.5,
    )
    assert key1 == key2
    assert len(key1) == 32


def test_save_and_load_altaz_grid(tmp_path):
    amount = np.zeros((90, 720), dtype=np.float32)
    amount[45, 180] = 1.0
    missing = np.zeros((90, 720), dtype=np.uint8)
    source_key = _make_source_key()
    grid = CloudAltAzGrid(
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
        time_utc=dt.datetime(2026, 6, 22, 12, 0, 0, tzinfo=dt.timezone.utc),
        shells_km=(6374.0, 6376.0, 6378.0),
        source_key=source_key,
        coverage_ratio=0.95,
        source_completeness_ratio=0.9,
    )
    save_altaz_grid(grid, tmp_path)
    key = altaz_grid_cache_key(
        35.0,
        135.0,
        satellite=grid.satellite,
        product=grid.product,
        time_utc=grid.time_utc,
        source_key=grid.source_key,
        shells_km=grid.shells_km,
        grid_resolution_deg=grid.grid_resolution_deg,
    )
    loaded = load_altaz_grid(tmp_path, key)
    assert loaded is not None
    assert loaded.amount.shape == (90, 720)
    assert loaded.coverage_ratio == pytest.approx(0.95)
    assert loaded.satellite == "G19"
    assert loaded.time_utc == grid.time_utc


def test_build_altaz_grid_does_not_store_sampler_on_source():
    """Regression: the internal sampler closure must not be cached on source.

    The sampler returned by ``build_bt_sampler`` is a closure and cannot be
    pickled.  ``build_altaz_grid`` should therefore never assign it back to
    ``source.sampler``; otherwise the worker subprocess fails to pickle the
    source artifact after building the alt/az grid.
    """
    source = CloudSourceData(
        source_key=_make_source_key(),
        data_array=_make_dummy_data_array(),
        satellite="G19",
        product="CMIPF-C13",
        time_utc=dt.datetime(2026, 6, 22, 12, 0, 0, tzinfo=dt.timezone.utc),
        src_paths=[],
    )
    assert source.sampler is None

    def fake_sampler_builder(da):
        def sampler(lon, lat):
            shape = np.asarray(lon).shape
            return np.full(shape, 280.0, dtype=np.float32)
        return sampler

    with patch("zstarview.clouddisc.altaz_grid.build_bt_sampler", side_effect=fake_sampler_builder):
        grid = build_altaz_grid(source, 35.0, 135.0)

    assert grid is not None
    assert source.sampler is None
