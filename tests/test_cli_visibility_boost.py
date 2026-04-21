from __future__ import annotations

import math
from datetime import timedelta

from zstarview.cli.args import parse_args
from zstarview.gui.window_inputs import prepare_window_runtime_options


def test_parse_args_visibility_boost_default(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    args = parse_args()

    assert math.isclose(float(args.visibility_boost), 1.0, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_visibility_boost_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--visibility-boost", "1.5"])
    args = parse_args()

    assert math.isclose(float(args.visibility_boost), 1.5, rel_tol=0.0, abs_tol=1e-9)


def test_prepare_window_runtime_options_boosts_cloud_missing_tint() -> None:
    options = prepare_window_runtime_options(
        delta_t=timedelta(0),
        sky_update_interval=60,
        urban_outline_radius_km=2.5,
        urban_outline_skyscraper_radius_km=60.0,
        urban_outline_min_height_m=0.0,
        urban_outline_feature_type="both",
        urban_outline_skyscraper_only=False,
        cloud_stripe_style=(50, 0.85),
        cloud_stripe_mode="width",
        cloud_missing_tint_opacity=0.2,
        visibility_boost=2.0,
        star_render_expected_width=600,
        content_fov_deg=115.0,
        window_geometry_arg=None,
        window_frame_mode="frameless",
    )

    assert math.isclose(float(options.cloud_missing_tint_opacity), 0.4, rel_tol=0.0, abs_tol=1e-9)
