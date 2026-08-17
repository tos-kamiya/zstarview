from __future__ import annotations

from types import SimpleNamespace

from zstarview.gui import viewer


def test_apply_gui_profile_to_args_parses_profile_values() -> None:
    args = SimpleNamespace()
    profile = {
        "city": "Tokyo",
        "place": "",
        "cloud_stripe": "alpha,12,0.5",
        "precipitation_opacity": 0.65,
        "window_geometry": "restore",
        "view_center_alt": 45.0,
        "view_center_az": 200.0,
    }

    viewer._apply_gui_profile_to_args(args, profile)

    assert args.city == "Tokyo"
    assert args.place is None
    assert args.cloud_stripe == ("alpha", 12, 0.5)
    assert args.precipitation_opacity == 0.65
    assert args.window_geometry == "restore"
    assert args.view_center_alt_specified is True
    assert args.view_center_az_specified is True


def test_apply_gui_profile_to_args_ignores_legacy_enlarge_moon() -> None:
    args = SimpleNamespace(enlarge_moon=False, moon_style="marker", moon_scale=1)

    viewer._apply_gui_profile_to_args(
        args,
        {"enlarge_moon": True, "moon_style": "image", "moon_scale": 4},
    )

    assert args.enlarge_moon is False
    assert args.moon_style == "image"
    assert args.moon_scale == 4


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

    viewer._apply_gui_profile_to_args(args, profile)

    assert args.city == ""
    assert args.place is None
    assert args.place_countrycode is None
    assert args.place_lang == "en"


def test_apply_gui_profile_to_args_preserves_explicit_cli_values(monkeypatch) -> None:
    args = SimpleNamespace(
        aircraft_opacity=0.0,
        city="CLI city",
        place="CLI place",
        place_countrycode="us",
        place_lang="ja",
    )
    profile = {
        "aircraft_opacity": 0.4,
        "city": {
            "resolver": "nominatim",
            "query": "Matsue Station",
            "result": {
                "name": "Matsue Station, Matsue, Shimane, Japan",
                "lat": 35.4641778,
                "lon": 133.0628539,
            },
        },
        "place": "Profile place",
        "place_countrycode": "jp",
        "place_lang": "en",
    }

    monkeypatch.setattr(
        viewer,
        "default_gui_launch_profile",
        lambda: {
            "aircraft_opacity": 0.4,
            "city": "",
            "place": None,
            "place_countrycode": None,
            "place_lang": "en",
        },
    )

    viewer._apply_gui_profile_to_args(args, profile)

    assert args.aircraft_opacity == 0.0
    assert args.city == "CLI city"
    assert args.place == "CLI place"
    assert args.place_countrycode == "us"
    assert args.place_lang == "ja"
