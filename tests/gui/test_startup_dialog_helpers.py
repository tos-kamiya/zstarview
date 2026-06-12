from __future__ import annotations

import pytest

from PySide6.QtWidgets import QFormLayout
from PySide6.QtWidgets import QDoubleSpinBox
from PySide6.QtWidgets import QCheckBox
from PySide6.QtWidgets import QApplication

import zstarview.gui.startup_dialog as startup_dialog_module
from zstarview.gui.startup_dialog import _coerce_float_value, _coerce_int_value
from zstarview.gui.startup_dialog import StartupDialog
from zstarview.location_resolver import ResolvedLocation
from zstarview.location_resolver import PlaceSearchCandidate
from zstarview.search.models import SearchJumpTarget

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
    assert dialog._location_place_radio.text() == "Search results"
    assert dialog._location_city_radio.isChecked() is True
    assert dialog._location_place_radio.isChecked() is False
    assert dialog._widgets["city"].isEnabled() is True
    assert dialog._widgets["city"].minimumHeight() == 40
    assert "place" not in dialog._widgets
    assert "place_countrycode" not in dialog._widgets
    assert "place_lang" not in dialog._widgets
    assert dialog._city_auto_button.text() == "Auto Search"
    assert dialog._city_search_button.text() == "Search ..."
    assert dialog._city_search_button.isEnabled() is False
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
    assert list(dialog._overlay_sections) == [
        "Sky",
        "Clouds",
        "Tropical Cyclone",
        "Aircraft and Satellites",
        "Ground and Guides",
        "Urban Outline",
    ]
    assert dialog._overlay_sections["Sky"].is_expanded() is True
    assert "tropical_cyclone_opacity" in dialog._widgets
    assert dialog._overlay_section_by_key["tropical_cyclone_opacity"] == "Tropical Cyclone"
    geo_widget = dialog._widgets["geo_satellite"]
    assert isinstance(geo_widget, QCheckBox)
    assert geo_widget.isChecked() is False
    general_layout = dialog._tab_layouts["General"]
    assert [
        general_layout.itemAt(index, QFormLayout.ItemRole.LabelRole).widget().text()
        for index in range(general_layout.rowCount())
    ] == [
        "Theme",
        "Window geometry",
        "Window frame",
        "Observation info",
        "Visibility boost",
        "Overlay font size",
        "Geo-satellite",
    ]
    dialog._overlay_sections["Sky"]._button.setChecked(False)
    assert dialog._overlay_sections["Sky"].is_expanded() is False
    assert "overlay_font_size" in dialog._widgets
    assert "overlay_font_size" not in dialog._overlay_section_by_key
    assert dialog._reset_button.text() == "Reset to Default Values"
    terrain_widget = dialog._widgets["terrain_horizon_opacity"]
    assert isinstance(terrain_widget, QDoubleSpinBox)
    assert terrain_widget.decimals() == 3
    assert terrain_widget.value() == 0.003
    ground_layout = dialog._overlay_layouts["Ground and Guides"]
    assert isinstance(ground_layout, QFormLayout)
    assert [
        ground_layout.itemAt(index, QFormLayout.ItemRole.LabelRole).widget().text()
        for index in range(ground_layout.rowCount())
        ] == [
            "Terrain horizon opacity",
            "Ridge glow opacity",
            "Earth guide opacity",
            "Ground tint opacity",
            "Water surface opacity",
        "Night light opacity",
    ]
    cyclone_widget = dialog._widgets["tropical_cyclone_opacity"]
    assert isinstance(cyclone_widget, QDoubleSpinBox)
    assert 0.0 <= cyclone_widget.value() <= 1.0


def test_startup_dialog_location_mode_checkboxes_are_exclusive() -> None:
    dialog = StartupDialog()

    dialog._location_place_radio.setChecked(True)

    assert dialog._location_place_radio.isChecked() is True
    assert dialog._location_city_radio.isChecked() is False
    assert dialog._widgets["city"].isEnabled() is False
    assert dialog._city_auto_button.isEnabled() is False
    assert dialog._city_search_button.isEnabled() is True

    dialog._location_city_radio.setChecked(True)

    assert dialog._location_city_radio.isChecked() is True
    assert dialog._location_place_radio.isChecked() is False
    assert dialog._widgets["city"].isEnabled() is True
    assert dialog._city_auto_button.isEnabled() is True
    assert dialog._city_search_button.isEnabled() is False


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
    city_widget.setPlainText("Tokyo")

    dialog._city_auto_button.click()

    assert calls == ["called"]
    assert city_widget.toPlainText() == "JP/Matsue"
    assert dialog._city_auto_button.isEnabled() is True


def test_startup_dialog_place_search_uses_search_dialog_and_saves_selection(monkeypatch) -> None:
    recorded_calls: list[tuple[str, str | None, str]] = []

    def fake_search_place_candidates(
        query: str,
        *,
        countrycode: str | None = None,
        language: str = "en",
        limit: int = 5,
        user_agent: str = "",
    ) -> tuple[PlaceSearchCandidate, ...]:
        recorded_calls.append((query, countrycode, language))
        del limit, user_agent
        return (
            PlaceSearchCandidate(
                name="Matsue Station",
                display_name="Matsue Station, Matsue, Shimane, Japan",
                latitude_deg=35.4641778,
                longitude_deg=133.0628539,
                category="railway",
                type_name="station",
                importance=0.9,
            ),
        )

    class FakePlaceSearchDialog:
        next_query = "Matsue Station"
        next_countrycode = None
        next_language = "en"

        def __init__(
            self,
            place_search_callback,
            parent=None,
            *,
            initial_query: str = "",
            initial_countrycode: str = "",
            initial_language: str = "en",
            show_cli_view_center_checkbox: bool = True,
            **_kwargs,
        ) -> None:
            del parent, initial_query, initial_countrycode, initial_language, show_cli_view_center_checkbox
            self._place_search_callback = place_search_callback
            self._targets: list[SearchJumpTarget] = []

        def setWindowTitle(self, _title: str) -> None:
            pass

        def exec(self) -> int:
            self._targets = list(
                self._place_search_callback(
                    self.next_query,
                    self.next_countrycode,
                    self.next_language,
                )
            )
            return 1

        def selected_target(self) -> SearchJumpTarget | None:
            return self._targets[0] if self._targets else None

        def search_query(self) -> str:
            return self.next_query

        def search_countrycode(self) -> str | None:
            return self.next_countrycode

        def search_language(self) -> str:
            return self.next_language

    monkeypatch.setattr(startup_dialog_module, "PlaceSearchDialog", FakePlaceSearchDialog)
    monkeypatch.setattr(startup_dialog_module, "search_place_candidates", fake_search_place_candidates)

    dialog = StartupDialog(
        profile={
            "window_geometry": "restore",
            "cloud_stripe": "width,50,0.85",
        }
    )
    dialog._location_place_radio.setChecked(True)

    dialog._city_search_button.click()

    assert recorded_calls == [("Matsue Station", None, "en")]
    assert dialog._widgets["city"].toPlainText() == "Matsue Station, Matsue, Shimane, Japan"
    profile = dialog.selected_profile()
    assert isinstance(profile["city"], dict)
    assert profile["city"]["resolver"] == "nominatim"
    assert profile["city"]["query"] == "Matsue Station"
    assert dialog._place_selected_payload is not None


def test_startup_dialog_preserves_tropical_cyclone_opacity_in_profile() -> None:
    dialog = StartupDialog(
        profile={
            "window_geometry": "restore",
            "cloud_stripe": "width,50,0.85",
            "tropical_cyclone_opacity": 0.27,
        }
    )

    cyclone_widget = dialog._widgets["tropical_cyclone_opacity"]
    assert isinstance(cyclone_widget, QDoubleSpinBox)
    assert cyclone_widget.value() == pytest.approx(0.27)

    profile = dialog.selected_profile()
    assert profile["tropical_cyclone_opacity"] == pytest.approx(0.27)


def test_startup_dialog_restores_geo_satellite_checkbox_from_profile() -> None:
    dialog = StartupDialog(
        profile={
            "window_geometry": "restore",
            "cloud_stripe": "width,50,0.85",
            "geo_satellite": True,
        }
    )

    geo_widget = dialog._widgets["geo_satellite"]
    assert isinstance(geo_widget, QCheckBox)
    assert geo_widget.isChecked() is True

    profile = dialog.selected_profile()
    assert profile["geo_satellite"] is True
