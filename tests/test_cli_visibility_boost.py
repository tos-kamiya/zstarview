from __future__ import annotations

import math
import sys
from datetime import timedelta

from zstarview.cli.args import parse_args
from zstarview.cloud_voxel_style import CloudVoxelMarkStyle
from zstarview.gui.window_inputs import (
    SkyWindowRuntimeOptions,
    prepare_window_runtime_options,
)


def test_parse_args_visibility_boost_default(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    args = parse_args()

    assert math.isclose(float(args.visibility_boost), 1.0, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_visibility_boost_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--visibility-boost", "1.5"])
    args = parse_args()

    assert math.isclose(float(args.visibility_boost), 1.5, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_asterism_opacity(monkeypatch) -> None:
    monkeypatch.setattr(sys, "argv", ["zstarview", "--asterism-opacity", "0.25"])
    args = parse_args()

    assert math.isclose(float(args.asterism_opacity), 0.25, rel_tol=0.0, abs_tol=1e-9)


def test_prepare_window_runtime_options_keeps_cloud_mark_style() -> None:
    style = CloudVoxelMarkStyle("cutout", 8)
    options = prepare_window_runtime_options(
        delta_t=timedelta(0),
        sky_update_interval=60,
        urban_outline_radius_km=2.5,
        urban_outline_skyscraper_radius_km=60.0,
        urban_outline_min_height_m=0.0,
        urban_outline_max_candidates=5000,
        urban_outline_feature_type="both",
        urban_outline_skyscraper_only=False,
        cloud_voxel_mark_style=style,
        visibility_boost=2.0,
        star_render_expected_width=600,
        content_fov_deg=110.0,
        window_geometry_arg=None,
        window_frame_mode="frameless",
    )

    assert options.cloud_voxel_mark_style == style


def test_sky_window_runtime_options_cloud_mark_default_matches_cli(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    assert SkyWindowRuntimeOptions().cloud_voxel_mark_style == CloudVoxelMarkStyle.from_cli(
        parse_args().cloud_stripe
    )
