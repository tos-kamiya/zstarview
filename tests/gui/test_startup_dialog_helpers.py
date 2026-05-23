from __future__ import annotations

import pytest

from PySide6.QtWidgets import QDoubleSpinBox
from PySide6.QtWidgets import QApplication

import zstarview.gui.startup_dialog as startup_dialog_module
from zstarview.gui.startup_dialog import _coerce_float_value, _coerce_int_value
from zstarview.gui.startup_dialog import StartupDialog
from zstarview.location_resolver import ResolvedLocation

_app = QApplication.instance() or QApplication([])


def test_startup_dialog_numeric_coercion_uses_fallbacks() -> None:
    assert _coerce_float_value("observer_height_m", None, 0.0) == 1.7
    assert _coerce_float_value("other", None, 2.5) == 2.5
    assert _coerce_int_value("other", None, 7) == 7


def test_startup_dialog_tabs_follow_requested_order() -> None:
    dialog = StartupDialog()

    assert dialog.width() == 640
    assert dialog.height() == 380
    assert dialog._tabs.count() == 6
    assert dialog._tabs.tabText(0) == "Location"
    assert dialog._tabs.tabText(1) == "Time"
    assert dialog._tabs.tabText(2) == "Stars"
    assert dialog._tabs.tabText(3) == "Overlays"
    assert dialog._tabs.tabText(4) == "General"
    assert dialog._tabs.tabText(5) == "Search Objects at Startup"
    assert dialog._widgets["time_shift_heading"].text() == "<b>Time shift</b>"
    assert dialog._widgets["time_absolute_heading"].text() == "<b>Absolute time</b>"
    assert set(dialog._overlay_sections) == {
        "Sky",
        "Clouds",
        "Aircraft and Satellites",
        "Ground and Guides",
        "Urban Outline",
    }
    assert dialog._overlay_sections["Sky"].is_expanded() is True
    dialog._overlay_sections["Sky"]._button.setChecked(False)
    assert dialog._overlay_sections["Sky"].is_expanded() is False
    assert "overlay_font_size" in dialog._widgets
    assert "overlay_font_size" not in dialog._overlay_section_by_key
    assert dialog._reset_button.text() == "Reset to Default Values"
    terrain_widget = dialog._widgets["terrain_horizon_opacity"]
    assert isinstance(terrain_widget, QDoubleSpinBox)
    assert terrain_widget.decimals() == 3
    assert terrain_widget.value() == 0.003


def test_startup_dialog_time_modes_are_mutually_exclusive() -> None:
    dialog = StartupDialog()
    base_profile = {
        "window_geometry": "restore",
        "cloud_stripe": "width,50,0.85",
    }

    with pytest.raises(ValueError, match="Relative shift and absolute time"):
        dialog._validate_profile(
            {
                **base_profile,
                "hours": 1.0,
                "datetime": "2026-05-23 12:00",
            }
        )

    with pytest.raises(ValueError, match="Timezone requires Date/time"):
        dialog._validate_profile(
            {
                **base_profile,
                "timezone": "Asia/Tokyo",
            }
        )


def test_startup_dialog_city_auto_button_fills_city(monkeypatch) -> None:
    monkeypatch.setattr(
        startup_dialog_module,
        "submit_gui_work",
        lambda target, *args, **kwargs: target(*args, **kwargs),
    )

    calls: list[str] = []

    def fake_auto_resolver() -> ResolvedLocation:
        calls.append("called")
        return ResolvedLocation(
            display_name="JP/Matsue",
            lat=35.4723,
            lon=133.0505,
            tz="Asia/Tokyo",
            persistence_key="auto:35.472300,133.050500",
            observer_height_m=1.7,
            kind="auto",
            cc="JP",
        )

    dialog = StartupDialog(auto_location_resolver=fake_auto_resolver)
    city_widget = dialog._widgets["city"]
    city_widget.setText("Tokyo")

    dialog._city_auto_button.click()

    assert calls == ["called"]
    assert city_widget.text() == "JP/Matsue"
    assert dialog._city_auto_button.isEnabled() is True
