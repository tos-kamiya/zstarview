# -*- coding: utf-8 -*-
"""Tests for halftone cloud sizing."""

from __future__ import annotations

from zstarview.gui.composite import _halftone_grid_delta, _halftone_level_diameters


def test_halftone_grid_delta_has_minimum_spacing() -> None:
    assert _halftone_grid_delta(600.0, 24) == 25.0
    assert _halftone_grid_delta(100.0, 24) == 20.0


def test_halftone_level_diameters_scale_with_grid_spacing() -> None:
    diameters = _halftone_level_diameters(8.0, 1.0)
    assert diameters[0] == 0.0
    assert diameters[-1] > diameters[1]
    small = _halftone_level_diameters(8.0, 1.0)
    large = _halftone_level_diameters(120.0, 1.0)
    assert large[-1] > small[-1]
