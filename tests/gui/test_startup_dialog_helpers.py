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

    assert dialog.width() == 480
    assert dialog.height() == 456
    assert dialog._tabs.count() == 6
    assert dialog._tabs.tabText(0) == "Location"
    assert dialog._tabs.tabText(1) == "View"
    assert dialog._tabs.tabText(2) == "Time"
    assert dialog._tabs.tabText(3) == "Stars"
    assert dialog._tabs.tabText(4) == "Overlays"
    assert dialog._tabs.tabText(5) == "General"
    assert dialog._location_city_radio.text() == "City"
    assert dialog._location_place_radio.text() == "Place query"
    assert dialog._location_city_radio.isChecked() is True
    assert dialog._location_place_radio.isChecked() is False
    assert dialog._widgets["city"].isEnabled() is True
    assert dialog._widgets["place"].isEnabled() is False
    assert dialog._widgets["view_center_alt"].isEnabled() is True
    assert dialog._view_center_alt_hint_label.text() == "Alt value: 0 is horizontal, 90 is zenith."
    assert dialog._widgets["edge_fov_deg"].isEnabled() is True
    assert set(dialog._view_center_az_buttons) == {"N", "E", "S", "W"}
    assert all(button.width() == 28 for button in dialog._view_center_az_buttons.values())
    assert dialog._time_current_radio.text() == "Current time"
    assert dialog._time_relative_radio.text() == "Relative time"
    assert dialog._time_absolute_radio.text() == "Absolute time"
    assert dialog._time_current_radio.isChecked() is True
    assert dialog._time_relative_radio.isChecked() is False
    assert dialog._time_absolute_radio.isChecked() is False
    assert dialog._widgets["hours"].isEnabled() is False
    assert dialog._widgets["datetime"].isEnabled() is False
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


def test_startup_dialog_location_mode_checkboxes_are_exclusive() -> None:
    dialog = StartupDialog()

    dialog._location_place_radio.setChecked(True)

    assert dialog._location_place_radio.isChecked() is True
    assert dialog._location_city_radio.isChecked() is False
    assert dialog._widgets["city"].isEnabled() is False
    assert dialog._widgets["place"].isEnabled() is True

    dialog._location_city_radio.setChecked(True)

    assert dialog._location_city_radio.isChecked() is True
    assert dialog._location_place_radio.isChecked() is False
    assert dialog._widgets["city"].isEnabled() is True
    assert dialog._widgets["place"].isEnabled() is False


def test_startup_dialog_view_center_az_buttons_set_cardinal_angles() -> None:
    dialog = StartupDialog()
    az_widget = dialog._widgets["view_center_az"]

    dialog._view_center_az_buttons["N"].click()
    assert az_widget.value() == 0.0

    dialog._view_center_az_buttons["E"].click()
    assert az_widget.value() == 90.0

    dialog._view_center_az_buttons["S"].click()
    assert az_widget.value() == 180.0

    dialog._view_center_az_buttons["W"].click()
    assert az_widget.value() == 270.0


def test_startup_dialog_time_source_radios_are_exclusive() -> None:
    dialog = StartupDialog()

    dialog._time_relative_radio.setChecked(True)

    assert dialog._time_current_radio.isChecked() is False
    assert dialog._time_relative_radio.isChecked() is True
    assert dialog._time_absolute_radio.isChecked() is False
    assert dialog._widgets["hours"].isEnabled() is True
    assert dialog._widgets["days"].isEnabled() is True
    assert dialog._widgets["datetime"].isEnabled() is False
    assert dialog._widgets["timezone"].isEnabled() is False

    dialog._time_absolute_radio.setChecked(True)

    assert dialog._time_current_radio.isChecked() is False
    assert dialog._time_relative_radio.isChecked() is False
    assert dialog._time_absolute_radio.isChecked() is True
    assert dialog._widgets["hours"].isEnabled() is False
    assert dialog._widgets["days"].isEnabled() is False
    assert dialog._widgets["datetime"].isEnabled() is True
    assert dialog._widgets["timezone"].isEnabled() is True

    dialog._time_current_radio.setChecked(True)

    assert dialog._time_current_radio.isChecked() is True
    assert dialog._time_relative_radio.isChecked() is False
    assert dialog._time_absolute_radio.isChecked() is False
    assert dialog._widgets["hours"].isEnabled() is False
    assert dialog._widgets["days"].isEnabled() is False
    assert dialog._widgets["datetime"].isEnabled() is False
    assert dialog._widgets["timezone"].isEnabled() is False


def test_startup_dialog_absolute_time_requires_datetime() -> None:
    dialog = StartupDialog()
    base_profile = {
        "window_geometry": "restore",
        "cloud_stripe": "width,50,0.85",
    }

    dialog._time_absolute_radio.setChecked(True)
    with pytest.raises(ValueError, match="Date/time is required when Absolute time is selected"):
        dialog._validate_profile(
            {
                **base_profile,
                "timezone": "Asia/Tokyo",
            }
        )


def test_startup_dialog_relative_time_requires_shift() -> None:
    dialog = StartupDialog()
    base_profile = {
        "window_geometry": "restore",
        "cloud_stripe": "width,50,0.85",
    }

    dialog._time_relative_radio.setChecked(True)
    with pytest.raises(ValueError, match="Relative time requires Hours or Days"):
        dialog._validate_profile(base_profile)


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
