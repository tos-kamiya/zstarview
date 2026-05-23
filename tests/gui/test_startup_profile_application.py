from __future__ import annotations

from types import SimpleNamespace

from zstarview.gui.viewer import _apply_gui_profile_to_args


def test_apply_gui_profile_to_args_parses_profile_values() -> None:
    args = SimpleNamespace()
    profile = {
        "city": "Tokyo",
        "cloud_stripe": "alpha,12,0.5",
        "window_geometry": "restore",
        "view_center_alt": 45.0,
        "view_center_az": 200.0,
    }

    _apply_gui_profile_to_args(args, profile)

    assert args.city == "Tokyo"
    assert args.cloud_stripe == ("alpha", 12, 0.5)
    assert args.window_geometry == "restore"
    assert args.view_center_alt_specified is True
    assert args.view_center_az_specified is True
