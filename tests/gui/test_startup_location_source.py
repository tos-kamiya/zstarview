from __future__ import annotations

import json
from types import SimpleNamespace

import pytest
from PySide6.QtWidgets import QApplication

from zstarview.gui import launch_profile
from zstarview.gui import viewer
from zstarview.gui.startup_dialog import StartupDialog

_app = QApplication.instance() or QApplication([])

_SEARCH_RESULT = {
    "resolver": "nominatim",
    "query": "Matsue Station",
    "display_name": "Matsue Station",
    "result": {"name": "Matsue Station", "lat": 35.4641778, "lon": 133.0628539},
}


@pytest.mark.parametrize(
    "profile,expected",
    [
        ({"location_source": "last", "city": "Tokyo"}, "last"),
        ({"city": _SEARCH_RESULT}, "place"),
        ({"city": {"resolver": "auto", "display_name": "JP/Matsue"}}, "city"),
        ({"city": "Matsue"}, "city"),
        ({"city": ""}, "auto"),
        ({}, "auto"),
        ({"location_source": "bogus", "city": "Matsue"}, "city"),
    ],
)
def test_location_source_from_profile(profile: dict, expected: str) -> None:
    assert launch_profile.location_source_from_profile(profile) == expected


def _write_profile(tmp_path, monkeypatch, payload: dict) -> None:
    monkeypatch.setattr(launch_profile, "user_config_dir", lambda *_args: str(tmp_path))
    (tmp_path / launch_profile.GUI_LAUNCH_PROFILE_FILENAME).write_text(
        json.dumps(payload), encoding="utf-8"
    )


def test_pre_v2_profile_drops_old_view_center_alt_default(tmp_path, monkeypatch) -> None:
    _write_profile(tmp_path, monkeypatch, {"view_center_alt": 90.0, "theme": "day"})

    loaded = launch_profile.load_gui_launch_profile()

    assert "view_center_alt" not in loaded
    assert "theme" not in loaded
    assert "profile_version" not in loaded


def test_pre_v2_profile_keeps_custom_view_center_alt(tmp_path, monkeypatch) -> None:
    _write_profile(tmp_path, monkeypatch, {"view_center_alt": 30.0, "theme": "night"})

    loaded = launch_profile.load_gui_launch_profile()

    assert loaded["view_center_alt"] == 30.0
    assert loaded["theme"] == "night"


def test_v2_profile_keeps_explicit_zenith_view(tmp_path, monkeypatch) -> None:
    _write_profile(
        tmp_path,
        monkeypatch,
        {"profile_version": 2, "view_center_alt": 90.0},
    )

    assert launch_profile.load_gui_launch_profile()["view_center_alt"] == 90.0


def test_saved_profile_records_version(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(launch_profile, "user_config_dir", lambda *_args: str(tmp_path))

    launch_profile.save_gui_launch_profile({"view_center_alt": 90.0})

    raw = json.loads(
        (tmp_path / launch_profile.GUI_LAUNCH_PROFILE_FILENAME).read_text(encoding="utf-8")
    )
    assert raw["profile_version"] == launch_profile.GUI_LAUNCH_PROFILE_VERSION
    assert launch_profile.load_gui_launch_profile()["view_center_alt"] == 90.0


def test_default_profile_uses_auto_location_and_low_view() -> None:
    defaults = launch_profile.default_gui_launch_profile()

    assert defaults["location_source"] == "auto"
    assert defaults["view_center_alt"] == 5.0


@pytest.mark.parametrize(
    "profile,expected_city",
    [
        ({"location_source": "auto", "city": "Matsue"}, ""),
        ({"location_source": "auto", "city": _SEARCH_RESULT}, ""),
        ({"location_source": "last", "city": "Matsue"}, "last"),
        ({"location_source": "city", "city": "Matsue"}, "Matsue"),
        ({"location_source": "place", "city": _SEARCH_RESULT}, "last"),
        ({"city": "Matsue"}, "Matsue"),
    ],
)
def test_apply_gui_profile_maps_location_source(profile: dict, expected_city: str) -> None:
    args = SimpleNamespace()

    viewer._apply_gui_profile_to_args(args, dict(profile))

    assert args.city == expected_city
    assert not hasattr(args, "location_source")


def test_apply_gui_profile_keeps_explicit_cli_location(monkeypatch) -> None:
    monkeypatch.setattr(
        viewer, "default_gui_launch_profile", lambda: {"city": "", "place": None}
    )
    args = SimpleNamespace(city="Osaka", place=None)

    viewer._apply_gui_profile_to_args(args, {"location_source": "auto", "city": "Matsue"})

    assert args.city == "Osaka"


def _dialog(profile: dict) -> StartupDialog:
    services = SimpleNamespace(submit=lambda target, *a, **k: target(*a, **k))
    return StartupDialog(profile=profile, services=services)  # type: ignore[arg-type]


@pytest.mark.parametrize("mode", ["auto", "last"])
def test_dialog_auto_and_last_keep_stored_location(mode: str) -> None:
    dialog = _dialog({"location_source": mode, "city": _SEARCH_RESULT})

    selected = dialog.selected_profile()

    assert selected["location_source"] == mode
    assert selected["city"] == _SEARCH_RESULT
    assert dialog._widgets["city"].isEnabled() is False


def test_dialog_switching_to_city_uses_typed_text() -> None:
    dialog = _dialog({"location_source": "auto", "city": ""})
    assert dialog._location_city_radio is not None
    dialog._location_city_radio.setChecked(True)
    dialog._widgets["city"].setPlainText("Matsue")  # type: ignore[attr-defined]

    selected = dialog.selected_profile()

    assert selected["location_source"] == "city"
    assert selected["city"] == "Matsue"


def test_old_profile_keeps_its_location_source(tmp_path, monkeypatch) -> None:
    _write_profile(tmp_path, monkeypatch, {"city": _SEARCH_RESULT})
    profile = dict(launch_profile.default_gui_launch_profile())
    profile.update(launch_profile.load_gui_launch_profile())

    assert profile["location_source"] == "place"
    dialog = _dialog(profile)
    assert dialog._location_place_radio is not None
    assert dialog._location_place_radio.isChecked() is True


def test_old_profile_with_city_text_stays_on_city(tmp_path, monkeypatch) -> None:
    _write_profile(tmp_path, monkeypatch, {"city": "Matsue"})

    assert launch_profile.load_gui_launch_profile()["location_source"] == "city"


@pytest.mark.parametrize("point_size", [10, 16])
def test_dialog_width_fits_tab_contents(point_size: int) -> None:
    from PySide6.QtWidgets import QScrollArea

    class _LargeScreenDialog(StartupDialog):
        def screen(self):  # type: ignore[override]
            return None  # no screen-width cap, like a large desktop

    app = QApplication.instance()
    assert isinstance(app, QApplication)
    original_font = app.font()
    font = app.font()
    font.setPointSize(point_size)
    app.setFont(font)
    try:
        services = SimpleNamespace(submit=lambda target, *a, **k: target(*a, **k))
        dialog = _LargeScreenDialog(profile={}, services=services)  # type: ignore[arg-type]
        dialog.show()
        app.processEvents()
        for index in range(dialog._tabs.count()):
            dialog._tabs.setCurrentIndex(index)
            app.processEvents()
            scroll_area = dialog._tabs.widget(index)
            assert isinstance(scroll_area, QScrollArea)
            assert scroll_area.horizontalScrollBar().isVisible() is False
        dialog.close()
    finally:
        app.setFont(original_font)
