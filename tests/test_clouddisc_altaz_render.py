# -*- coding: utf-8 -*-
"""Tests for the alt/az grid circle renderer."""

from __future__ import annotations

import datetime as dt

import numpy as np

from zstarview.clouddisc.altaz_grid import CloudAltAzGrid
from zstarview.clouddisc.altaz_render import (
    render_altaz_grid_circles,
    render_altaz_missing_mask,
)
from zstarview.clouddisc.types import SourceKey


def _make_grid(amount_pattern=None):
    amount = np.zeros((90, 720), dtype=np.float32)
    if amount_pattern is not None:
        amount_pattern(amount)
    missing = np.zeros((90, 720), dtype=np.uint8)
    source_key = SourceKey(
        satellite="G19",
        provider="GOES",
        timeslot_utc=dt.datetime(2026, 6, 22, 12, 0, 0, tzinfo=dt.timezone.utc),
        sat_priority=("AUTO",),
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
        time_utc=dt.datetime(2026, 6, 22, 12, 0, 0, tzinfo=dt.timezone.utc),
        shells_km=(6374.0, 6376.0, 6378.0),
        source_key=source_key,
        coverage_ratio=1.0,
    )


def test_render_empty_grid_is_transparent():
    grid = _make_grid()
    img = render_altaz_grid_circles(
        grid, 600, 600, center_alt_deg=45.0, center_az_deg=180.0, edge_fov_deg=90.0
    )
    assert img.shape == (600, 600, 4)
    assert np.all(img[..., 3] == 0)


def test_render_single_cell_produces_white_circle():
    def pattern(amount):
        amount[45, 360] = 1.0  # mid-altitude, south-ish direction

    grid = _make_grid(pattern)
    img = render_altaz_grid_circles(
        grid, 600, 600, center_alt_deg=45.0, center_az_deg=180.0, edge_fov_deg=90.0
    )
    assert img.shape == (600, 600, 4)
    assert np.any(img[..., 3] > 0)
    # RGB should be white wherever alpha is positive.
    positive = img[..., 3] > 0
    assert np.all(img[..., :3][positive] == 255)


def test_render_amount_scales_alpha():
    def pattern(amount):
        amount[30, 360] = 0.25
        amount[60, 360] = 1.0

    grid = _make_grid(pattern)
    img = render_altaz_grid_circles(
        grid, 600, 600, center_alt_deg=45.0, center_az_deg=180.0, edge_fov_deg=90.0
    )
    alpha = img[..., 3]
    nonzero = alpha > 0
    assert nonzero.any()
    assert alpha.max() > 100  # strong cell should be visible


def test_render_missing_mask_shape_and_values():
    def pattern(amount):
        amount[45, 360] = 1.0

    def missing_pattern(missing):
        missing[45, 360] = 255
        missing[45, 361] = 255

    grid = _make_grid(pattern)
    grid.missing_mask[:] = 0
    missing_pattern(grid.missing_mask)

    mask = render_altaz_missing_mask(
        grid, 600, 600, center_alt_deg=45.0, center_az_deg=180.0, edge_fov_deg=90.0
    )
    assert mask.shape == (600, 600)
    assert mask.dtype == np.uint8
    assert set(np.unique(mask)).issubset({0, 255})


def test_render_respects_mask_fov():
    def pattern(amount):
        amount[45, 360] = 1.0

    grid = _make_grid(pattern)
    # View centre is far from the cell direction; with a narrow FOV the cell
    # should not be drawn.
    img = render_altaz_grid_circles(
        grid,
        600,
        600,
        center_alt_deg=80.0,
        center_az_deg=0.0,
        edge_fov_deg=10.0,
        mask_fov_deg=10.0,
    )
    assert np.all(img[..., 3] == 0)


def test_missing_mask_handles_cells_projected_off_screen():
    """Regression: missing cells near the image edge must not raise broadcast errors."""
    amount = np.zeros((90, 720), dtype=np.float32)
    missing = np.zeros((90, 720), dtype=np.uint8)
    # Mark a cell that projects beyond the right image edge with the chosen FOV.
    missing[45, 0] = 255

    source_key = SourceKey(
        satellite="G19",
        provider="GOES",
        timeslot_utc=dt.datetime(2026, 6, 22, 12, 0, 0, tzinfo=dt.timezone.utc),
        sat_priority=("AUTO",),
    )
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
        coverage_ratio=1.0,
    )
    # Small image + wide mask FOV lets the missing cell project outside the canvas.
    mask = render_altaz_missing_mask(
        grid,
        20,
        20,
        center_alt_deg=45.0,
        center_az_deg=0.0,
        edge_fov_deg=10.0,
        mask_fov_deg=90.0,
    )
    assert mask.shape == (20, 20)
    assert mask.dtype == np.uint8
