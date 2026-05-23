from __future__ import annotations

from types import SimpleNamespace

from zstarview.gui.viewer import _apply_gui_profile_to_args


def test_apply_gui_profile_to_args_parses_profile_values() -> None:
    args = SimpleNamespace()
    profile = {
        "city": "Tokyo",
        "place": "",
        "cloud_stripe": "alpha,12,0.5",
        "window_geometry": "restore",
        "view_center_alt": 45.0,
        "view_center_az": 200.0,
    }

    _apply_gui_profile_to_args(args, profile)

    assert args.city == "Tokyo"
    assert args.place is None
    assert args.cloud_stripe == ("alpha", 12, 0.5)
    assert args.window_geometry == "restore"
    assert args.view_center_alt_specified is True
    assert args.view_center_az_specified is True


def test_apply_gui_profile_to_args_ignores_structured_city_payload_for_cli_args() -> None:
    args = SimpleNamespace()
    profile = {
        "city": {
            "resolver": "nominatim",
            "query": "Matsue Station",
            "result": {
                "name": "Matsue Station, Matsue, Shimane, Japan",
                "lat": 35.4641778,
                "lon": 133.0628539,
            },
        },
        "place": "Matsue Station",
        "place_countrycode": "jp",
        "place_lang": "",
        "window_geometry": "restore",
    }

    _apply_gui_profile_to_args(args, profile)

    assert args.city == ""
    assert args.place is None
    assert args.place_countrycode is None
    assert args.place_lang == "en"
