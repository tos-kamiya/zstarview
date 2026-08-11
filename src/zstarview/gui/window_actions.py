from __future__ import annotations

from collections.abc import Callable
from datetime import datetime, timedelta, timezone

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QAction, QDesktopServices, QKeySequence
from PySide6.QtWidgets import QMenu

from ..__about__ import __version__
from ..overlay_time import overlay_availability_for_delta
from ..paths import (
    CLOUD_UPDATE_INTERVAL,
    OBSERVER_MAX_ALT_DEG,
    OBSERVER_MIN_ALT_DEG,
)
from .license_dialog import LicenseDialog

OPEN_METEO_TERMS_URL = "https://open-meteo.com/en/terms"


def open_open_meteo_terms() -> None:
    QDesktopServices.openUrl(QUrl(OPEN_METEO_TERMS_URL))


class SkyWindowActionsMixin:
    def _build_window_menu(self) -> None:
        """Build window actions and the popup menu shared by chrome implementations."""
        from . import window as window_module

        menu_class = window_module.QMenu
        self.menu = menu_class("Menu", self)
        self.file_menu = menu_class("File", self)
        self.search_menu = menu_class("Search", self)
        self.observer_view_menu = menu_class("View Direction", self)
        self.display_menu = menu_class("Layers", self)
        self.help_menu = menu_class("Help", self)
        if not self._frameless_window:
            self.menu.addMenu(self.file_menu)
        self.menu.addMenu(self.search_menu)
        self.menu.addMenu(self.display_menu)
        self.menu.addMenu(self.observer_view_menu)
        self.menu.addMenu(self.help_menu)

        self._add_menu_action(
            self.observer_view_menu,
            f"Rotate Left (-{self.state.rotation_step:.0f}°)",
            triggered=lambda: self._rotate_view(
                d_az=-self.state.rotation_step, interactive_viewport=True
            ),
        )
        self._add_menu_action(
            self.observer_view_menu,
            f"Rotate Right (+{self.state.rotation_step:.0f}°)",
            triggered=lambda: self._rotate_view(
                d_az=+self.state.rotation_step, interactive_viewport=True
            ),
        )
        self._action_raise_view = self._add_menu_action(
            self.observer_view_menu,
            f"Raise View (+{self.state.rotation_step:.0f}° alt)",
            triggered=lambda: self._rotate_view(
                d_alt=+self.state.rotation_step, interactive_viewport=True
            ),
        )
        self._action_lower_view = self._add_menu_action(
            self.observer_view_menu,
            f"Lower View (-{self.state.rotation_step:.0f}° alt)",
            triggered=lambda: self._rotate_view(
                d_alt=-self.state.rotation_step, interactive_viewport=True
            ),
        )

        self.observer_view_menu.addSeparator()
        self._add_menu_action(
            self.observer_view_menu,
            "Set View Center...",
            triggered=self._open_view_direction_dialog,
        )

        self.observer_view_menu.addSeparator()
        self._add_menu_action(
            self.search_menu,
            "Jump to Named Star...",
            shortcut=QKeySequence("Ctrl+J"),
            triggered=self._open_named_star_jump_dialog,
        )
        self._add_menu_action(
            self.search_menu,
            "Search Objects...",
            shortcut=QKeySequence("Ctrl+F"),
            triggered=self._open_named_star_search_dialog,
        )
        self._add_menu_action(
            self.search_menu,
            "Search Places...",
            triggered=self._open_place_search_dialog,
        )

        self.display_menu.addSeparator()
        self._action_enlarge_moon = self._add_checkable_menu_action(
            self.display_menu,
            "Enlarge Moon",
            checked=self.enlarge_moon,
            shortcut=QKeySequence(Qt.Key.Key_M),
            triggered=self.toggle_enlarge_moon,
        )
        self._action_toggle_dso = self._add_checkable_menu_action(
            self.display_menu,
            "DSO",
            checked=self.show_dso,
            enabled=self.dso_catalog_np is not None,
            shortcut=QKeySequence(Qt.Key.Key_D),
            triggered=self.toggle_dso,
        )
        if hasattr(self, "_akari_ir_bands_toggle_supported"):
            self._action_toggle_akari_ir_bands = self._add_checkable_menu_action(
                self.display_menu,
                "AKARI IR bands",
                checked=float(getattr(self, "akari_ir_bands_opacity", 0.10)) > 0.0,
                enabled=self._akari_ir_bands_toggle_supported,
                shortcut=QKeySequence(Qt.Key.Key_K),
                triggered=getattr(self, "toggle_akari_ir_bands", lambda: None),
            )
        else:
            self._action_toggle_akari_ir_bands = None
        self._action_toggle_asterisms = self._add_checkable_menu_action(
            self.display_menu,
            "Asterisms",
            checked=self.show_asterisms,
            shortcut=QKeySequence(Qt.Key.Key_A),
            triggered=self.toggle_asterisms,
        )
        self._action_toggle_guidelines = self._add_checkable_menu_action(
            self.display_menu,
            "Sky Guides",
            checked=self.show_guidelines,
            shortcut=QKeySequence(Qt.Key.Key_G),
            triggered=self.toggle_guidelines,
        )
        self._action_toggle_observation_info = self._add_checkable_menu_action(
            self.display_menu,
            "Observation Info",
            checked=self.show_observation_info,
            enabled=self.observation_info_mode != "off",
            triggered=self.toggle_observation_info,
        )

        self.display_menu.addSeparator()
        self._action_toggle_sky_disc = self._add_checkable_menu_action(
            self.display_menu,
            "Sky Color",
            checked=self.sky_disc_alpha > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_S),
            triggered=self.toggle_sky_disc,
        )
        self._action_toggle_clouds = self._add_checkable_menu_action(
            self.display_menu,
            "Clouds",
            checked=self.cloud_disc_alpha > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_C),
            triggered=self.toggle_clouds,
        )
        self._action_toggle_geo_satellite = self._add_checkable_menu_action(
            self.display_menu,
            "Geo-satellite",
            checked=self._geo_satellite_enabled
            and self._geo_satellite_toggle_supported(),
            enabled=self._geo_satellite_toggle_supported(),
            triggered=self.toggle_geo_satellite,
        )
        self._action_toggle_satellites = self._add_checkable_menu_action(
            self.display_menu,
            "Satellites",
            checked=self.satellite_opacity > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_I),
            triggered=self.toggle_satellites,
        )
        self._action_toggle_aircraft = self._add_checkable_menu_action(
            self.display_menu,
            "Aircraft",
            checked=self.aircraft_opacity > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_P),
            triggered=self.toggle_aircraft,
        )
        self._action_toggle_tropical_cyclone = self._add_checkable_menu_action(
            self.display_menu,
            "Typhoon / Cyclone",
            checked=self.show_tropical_cyclone_overlay,
            enabled=self._tropical_cyclone_toggle_supported,
            triggered=self.toggle_tropical_cyclone_overlay,
        )
        self._action_toggle_precipitation = self._add_checkable_menu_action(
            self.display_menu,
            "Forecast Precipitation",
            checked=float(getattr(self, "precipitation_opacity", 0.0)) > 0.0,
            enabled=bool(getattr(self, "_precipitation_toggle_supported", False)),
            triggered=getattr(self, "toggle_precipitation", lambda: None),
        )
        self.display_menu.addSeparator()
        self._action_toggle_night_lights = self._add_checkable_menu_action(
            self.display_menu,
            "Night Lights",
            checked=self.night_light_opacity > 0.0,
            enabled=bool(self._night_light_toggle_supported),
            shortcut=QKeySequence(Qt.Key.Key_L),
            triggered=self.toggle_night_lights,
        )
        self._action_toggle_road_lights = self._add_checkable_menu_action(
            self.display_menu,
            "Road Lights",
            checked=self.road_night_lights_opacity > 0.0,
            enabled=self._road_light_toggle_supported,
            shortcut=QKeySequence(Qt.Key.Key_R),
            triggered=self.toggle_road_lights,
        )
        self._action_toggle_urban_outline = self._add_checkable_menu_action(
            self.display_menu,
            "Urban Outline",
            checked=self.urban_outline_opacity > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_U),
            triggered=self.toggle_urban_outline,
        )
        self._action_toggle_terrain_horizon = self._add_checkable_menu_action(
            self.display_menu,
            "Terrain Horizon",
            checked=self.terrain_horizon_opacity > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_T),
            triggered=self.toggle_terrain_horizon,
        )
        self._action_toggle_water_overlay = self._add_checkable_menu_action(
            self.display_menu,
            "Water Surface",
            checked=self.water_overlay_opacity > 0.0,
            enabled=self._water_overlay_action_enabled(),
            shortcut=QKeySequence(Qt.Key.Key_W),
            triggered=self.toggle_water_overlay,
        )
        self._action_toggle_earth_guide = self._add_checkable_menu_action(
            self.display_menu,
            "Earth Guide",
            checked=self.earth_guide_opacity > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_E),
            triggered=self.toggle_earth_guide,
        )

        square_window_action = self._add_menu_action(
            self.file_menu,
            "Square Window",
            triggered=self.toggle_square_window,
        )
        self._action_square_window = square_window_action
        restore_default_size_action = self._add_menu_action(
            self.file_menu,
            "Default Window Size",
            triggered=self._restore_default_window_size,
        )
        fit_to_screen_action = self._add_menu_action(
            self.file_menu,
            "Fit to Screen",
            triggered=self._fit_client_area_to_screen,
        )
        fullscreen_action = self._add_menu_action(
            self.file_menu,
            "Fullscreen",
            shortcut=QKeySequence(Qt.Key.Key_F11),
            triggered=self.toggle_fullscreen,
        )

        exit_action = self._add_menu_action(
            self.file_menu,
            "Exit",
            shortcut=QKeySequence(Qt.Key.Key_Q),
            triggered=self._request_application_quit,
        )

        self.display_menu.addSeparator()
        vmag_limit_action = self._add_menu_action(
            self.display_menu,
            self._vmag_limit_menu_text(),
        )
        vmag_limit_action.setEnabled(False)

        self._add_menu_action(
            self.help_menu,
            "Open-Meteo Terms...",
            triggered=open_open_meteo_terms,
        )
        self._add_menu_action(
            self.help_menu,
            "Licenses and Data Sources...",
            triggered=lambda: LicenseDialog(self).exec(),
        )
        version_action = self._add_menu_action(
            self.help_menu,
            f"Version {__version__}",
        )
        version_action.setEnabled(False)

        if self._frameless_window:
            self.menu.addSeparator()
            self.menu.addAction(square_window_action)
            self.menu.addAction(restore_default_size_action)
            self.menu.addAction(fit_to_screen_action)
            self.menu.addAction(fullscreen_action)
            self.menu.addSeparator()
            self.menu.addAction(exit_action)

    def _create_menu_action(self, menu: QMenu, text: str) -> QAction:
        """Create a menu action and add it to the target menu."""
        action = QAction(text, self)
        menu.addAction(action)
        return action

    def _add_menu_action(
        self,
        menu: QMenu,
        text: str,
        *,
        shortcut: QKeySequence | str | int | None = None,
        enabled: bool = True,
        triggered: Callable[..., object] | None = None,
    ) -> QAction:
        """Create a non-checkable menu action."""
        action = self._create_menu_action(menu, text)
        if shortcut is not None:
            action.setShortcut(QKeySequence(shortcut))
            action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
            self.addAction(action)
        action.setEnabled(enabled)
        if triggered is not None:
            action.triggered.connect(triggered)
        return action

    def _add_checkable_menu_action(
        self,
        menu: QMenu,
        text: str,
        *,
        checked: bool = False,
        shortcut: QKeySequence | str | int | None = None,
        enabled: bool = True,
        triggered: Callable[..., object] | None = None,
    ) -> QAction:
        """Create a checkable menu action."""
        action = self._create_menu_action(menu, text)
        action.setCheckable(True)
        action.setChecked(checked)
        if shortcut is not None:
            action.setShortcut(QKeySequence(shortcut))
            action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
            self.addAction(action)
        action.setEnabled(enabled)
        if triggered is not None:
            action.triggered.connect(triggered)
        return action
    def _vmag_limit_menu_text(self) -> str:
        return f"Vmag limit {self.vmag_limit:.1f}"
    def _sync_view_altitude_actions(self) -> None:
        alt, _ = self.viewer_data.view_center
        if self._action_raise_view is not None:
            self._action_raise_view.setEnabled(float(alt) < OBSERVER_MAX_ALT_DEG)
        if self._action_lower_view is not None:
            self._action_lower_view.setEnabled(float(alt) > OBSERVER_MIN_ALT_DEG)
    def _sync_geo_satellite_action_state(self) -> None:
        if self._action_toggle_geo_satellite is None:
            return
        supported = self._geo_satellite_toggle_supported()
        self._action_toggle_geo_satellite.setEnabled(supported)
        self._action_toggle_geo_satellite.setChecked(
            bool(self._geo_satellite_enabled) and supported
        )

    def _sync_cloud_action_state(self) -> None:
        if self._action_toggle_clouds is None:
            return
        supported = bool(
            self._cloud_toggle_supported
            and self._cloud_gui_allowed
            and (self._clouddisc is not None or self._geo_satellite_enabled)
        )
        self._action_toggle_clouds.setEnabled(supported)
        self._action_toggle_clouds.setChecked(
            supported and float(self.cloud_disc_alpha) > 0.0
        )
    def toggle_enlarge_moon(self) -> None:
        self.enlarge_moon = not self.enlarge_moon
        if (
            self._action_enlarge_moon is not None
            and self._action_enlarge_moon.isChecked() != self.enlarge_moon
        ):
            self._action_enlarge_moon.setChecked(self.enlarge_moon)
        self.request_client_update()

    def toggle_clouds(self) -> None:
        if not self._cloud_toggle_supported or not self._cloud_gui_allowed:
            self._sync_cloud_action_state()
            return

        enable_clouds = self.cloud_disc_alpha <= 0.0
        self.cloud_disc_alpha = self._cloud_alpha_when_enabled if enable_clouds else 0.0
        self._sync_cloud_action_state()

        if enable_clouds:
            self.request_cloud_projection_update(reason="toggle-on")
            self.state.cloud_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
                seconds=CLOUD_UPDATE_INTERVAL
            )
        else:
            self.state.cloud_next_refresh_utc = None

        self.request_client_update()

    def toggle_geo_satellite(self) -> None:
        if not self._geo_satellite_toggle_supported():
            self._sync_geo_satellite_action_state()
            return
        self._geo_satellite_enabled = not bool(self._geo_satellite_enabled)
        self._sync_geo_satellite_action_state()
        overlay_availability = overlay_availability_for_delta(self.delta_t)
        self._cloud_toggle_supported = overlay_availability.cloud and (
            self._clouddisc is not None or self._geo_satellite_enabled
        )
        if not self._cloud_toggle_supported:
            self.cloud_disc_alpha = 0.0
        self._sync_cloud_action_state()

        self.request_cloud_projection_update(reason="toggle-geo-satellite")

        self.state.cloud_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=CLOUD_UPDATE_INTERVAL
        )
        self.request_client_update()

    def toggle_satellites(self) -> None:
        if not self._satellite_toggle_supported or not self._satellite_gui_allowed:
            if self._action_toggle_satellites is not None:
                self._action_toggle_satellites.setChecked(self.satellite_opacity > 0.0)
            return

        enable_satellites = self.satellite_opacity <= 0.0
        self.satellite_opacity = (
            self._satellite_opacity_when_enabled if enable_satellites else 0.0
        )
        if (
            self._action_toggle_satellites is not None
            and self._action_toggle_satellites.isChecked() != enable_satellites
        ):
            self._action_toggle_satellites.setChecked(enable_satellites)

        if enable_satellites:
            self._enable_satellite_layer(reason="toggle-on")
        else:
            self._stop_satellite_timers()

        self.request_client_update()

    def toggle_aircraft(self) -> None:
        if not self._aircraft_toggle_supported or not self._aircraft_gui_allowed:
            if self._action_toggle_aircraft is not None:
                self._action_toggle_aircraft.setChecked(self.aircraft_opacity > 0.0)
            return

        enable_aircraft = self.aircraft_opacity <= 0.0
        self.aircraft_opacity = (
            self._aircraft_opacity_when_enabled if enable_aircraft else 0.0
        )
        if (
            self._action_toggle_aircraft is not None
            and self._action_toggle_aircraft.isChecked() != enable_aircraft
        ):
            self._action_toggle_aircraft.setChecked(enable_aircraft)

        if enable_aircraft:
            self._enable_aircraft_layer(reason="toggle-on")
        else:
            self._stop_aircraft_timers()

        self.request_client_update()

    def toggle_tropical_cyclone_overlay(self) -> None:
        if not self._tropical_cyclone_toggle_supported:
            if self._action_toggle_tropical_cyclone is not None:
                self._action_toggle_tropical_cyclone.setChecked(False)
            return

        if self.tropical_cyclone_opacity > 0.0:
            self._tropical_cyclone_opacity_when_enabled = self.tropical_cyclone_opacity
        enable_tropical_cyclone = self.tropical_cyclone_opacity <= 0.0
        self.tropical_cyclone_opacity = (
            self._tropical_cyclone_opacity_when_enabled
            if enable_tropical_cyclone
            else 0.0
        )
        self.show_tropical_cyclone_overlay = bool(self.tropical_cyclone_opacity > 0.0)
        if (
            self._action_toggle_tropical_cyclone is not None
            and self._action_toggle_tropical_cyclone.isChecked()
            != self.show_tropical_cyclone_overlay
        ):
            self._action_toggle_tropical_cyclone.setChecked(
                self.show_tropical_cyclone_overlay
            )
        cyclone_snapshots = self.tropical_cyclone_state.snapshots
        if self.show_tropical_cyclone_overlay and not cyclone_snapshots:
            self.start_background_tropical_cyclone_update(reason="toggle-on")
        self.request_client_update()

    def toggle_precipitation(self) -> None:
        if not self._precipitation_toggle_supported:
            if self._action_toggle_precipitation is not None:
                self._action_toggle_precipitation.setChecked(False)
            return

        enable_precipitation = self.precipitation_opacity <= 0.0
        self.precipitation_opacity = (
            self._precipitation_opacity_when_enabled
            if enable_precipitation
            else 0.0
        )
        if (
            self._action_toggle_precipitation is not None
            and self._action_toggle_precipitation.isChecked() != enable_precipitation
        ):
            self._action_toggle_precipitation.setChecked(enable_precipitation)
        if enable_precipitation:
            self.start_background_precipitation_update(reason="toggle-on")
        self._compositor.invalidate()
        self.request_client_update()

    def toggle_dso(self) -> None:
        if self.dso_catalog_np is None:
            self.show_dso = False
            if self._action_toggle_dso is not None:
                self._action_toggle_dso.setChecked(False)
            return
        self.show_dso = not self.show_dso
        if (
            self._action_toggle_dso is not None
            and self._action_toggle_dso.isChecked() != self.show_dso
        ):
            self._action_toggle_dso.setChecked(self.show_dso)
        self.request_client_update()

    def toggle_asterisms(self) -> None:
        self.show_asterisms = not self.show_asterisms
        if (
            self._action_toggle_asterisms is not None
            and self._action_toggle_asterisms.isChecked() != self.show_asterisms
        ):
            self._action_toggle_asterisms.setChecked(self.show_asterisms)
        self.request_client_update()

    def toggle_guidelines(self) -> None:
        self.show_guidelines = not self.show_guidelines
        if (
            self._action_toggle_guidelines is not None
            and self._action_toggle_guidelines.isChecked() != self.show_guidelines
        ):
            self._action_toggle_guidelines.setChecked(self.show_guidelines)
        self.request_client_update()

    def toggle_observation_info(self) -> None:
        # If CLI explicitly disabled the overlay, the menu toggle should be inert.
        if self.observation_info_mode == "off":
            self.show_observation_info = False
            if (
                self._action_toggle_observation_info is not None
                and self._action_toggle_observation_info.isChecked()
                != self.show_observation_info
            ):
                self._action_toggle_observation_info.setChecked(
                    self.show_observation_info
                )
            return

        # Toggle visibility. When re-enabling and the overlay is pinned, ensure
        # the pinned position is enforced.
        self.show_observation_info = not self.show_observation_info

        if self.show_observation_info and self.observation_info_pinned:
            mode = self.observation_info_mode
            try:
                if mode == "bottom":
                    self.state.overlay_info_bottom_left = True
                elif mode == "top":
                    self.state.overlay_info_bottom_left = False
            except Exception:
                pass

        if (
            self._action_toggle_observation_info is not None
            and self._action_toggle_observation_info.isChecked()
            != self.show_observation_info
        ):
            self._action_toggle_observation_info.setChecked(self.show_observation_info)
        self.request_client_update()

    def toggle_sky_disc(self) -> None:
        if not self._sky_disc_gui_allowed:
            if self._action_toggle_sky_disc is not None:
                self._action_toggle_sky_disc.setChecked(self.sky_disc_alpha > 0.0)
            return
        enable_sky_disc_gradient = self.sky_disc_alpha <= 0.0
        self.sky_disc_alpha = (
            self._sky_disc_alpha_when_enabled if enable_sky_disc_gradient else 0.0
        )
        if (
            self._action_toggle_sky_disc is not None
            and self._action_toggle_sky_disc.isChecked() != enable_sky_disc_gradient
        ):
            self._action_toggle_sky_disc.setChecked(enable_sky_disc_gradient)
        self._compositor.invalidate()
        self.request_sky_data_update()
        self.request_client_update()

    def toggle_terrain_horizon(self) -> None:
        if not self._terrain_horizon_gui_allowed:
            if self._action_toggle_terrain_horizon is not None:
                self._action_toggle_terrain_horizon.setChecked(
                    self.terrain_horizon_opacity > 0.0
                )
            self._sync_water_overlay_action_enabled()
            return

        enable_terrain = self.terrain_horizon_opacity <= 0.0
        self.terrain_horizon_opacity = (
            self._terrain_horizon_opacity_when_enabled if enable_terrain else 0.0
        )
        if (
            self._action_toggle_terrain_horizon is not None
            and self._action_toggle_terrain_horizon.isChecked() != enable_terrain
        ):
            self._action_toggle_terrain_horizon.setChecked(enable_terrain)
        self._refresh_water_overlay_active_dots()
        self._sync_water_overlay_action_enabled()
        self._compositor.invalidate()
        if enable_terrain:
            self.start_background_terrain_horizon_update(reason="toggle-on")
        self.request_client_update()

    def toggle_water_overlay(self) -> None:
        if not self._water_overlay_action_enabled():
            if self._action_toggle_water_overlay is not None:
                self._action_toggle_water_overlay.setChecked(
                    self.water_overlay_opacity > 0.0
                )
            return

        enable_water_overlay = self.water_overlay_opacity <= 0.0
        self.water_overlay_opacity = (
            self._water_overlay_opacity_when_enabled if enable_water_overlay else 0.0
        )
        self.show_water_overlay_layer = self.water_overlay_opacity > 0.0
        if (
            self._action_toggle_water_overlay is not None
            and self._action_toggle_water_overlay.isChecked() != enable_water_overlay
        ):
            self._action_toggle_water_overlay.setChecked(enable_water_overlay)
        self.request_client_update()

    def toggle_earth_guide(self) -> None:
        if not self._earth_guide_gui_allowed:
            if self._action_toggle_earth_guide is not None:
                self._action_toggle_earth_guide.setChecked(
                    self.earth_guide_opacity > 0.0
                )
            return

        enable_earth_guide = self.earth_guide_opacity <= 0.0
        self.earth_guide_opacity = (
            self._earth_guide_opacity_when_enabled if enable_earth_guide else 0.0
        )
        if (
            self._action_toggle_earth_guide is not None
            and self._action_toggle_earth_guide.isChecked() != enable_earth_guide
        ):
            self._action_toggle_earth_guide.setChecked(enable_earth_guide)
        self._compositor.invalidate()
        self.request_client_update()

    def toggle_night_lights(self) -> None:
        if not bool(self._night_light_toggle_supported):
            if self._action_toggle_night_lights is not None:
                self._action_toggle_night_lights.setChecked(
                    self.night_light_opacity > 0.0
                )
            return

        enable_night_lights = self.night_light_opacity <= 0.0
        self.night_light_opacity = (
            self._night_light_opacity_when_enabled if enable_night_lights else 0.0
        )
        if (
            self._action_toggle_night_lights is not None
            and self._action_toggle_night_lights.isChecked() != enable_night_lights
        ):
            self._action_toggle_night_lights.setChecked(enable_night_lights)
        self.request_client_update()

    def toggle_road_lights(self) -> None:
        if not bool(self._road_light_toggle_supported):
            if self._action_toggle_road_lights is not None:
                self._action_toggle_road_lights.setChecked(
                    self.road_night_lights_opacity > 0.0
                )
            return

        enable_road_lights = self.road_night_lights_opacity <= 0.0
        self.road_night_lights_opacity = (
            self._road_light_opacity_when_enabled if enable_road_lights else 0.0
        )
        if (
            self._action_toggle_road_lights is not None
            and self._action_toggle_road_lights.isChecked() != enable_road_lights
        ):
            self._action_toggle_road_lights.setChecked(enable_road_lights)
        if enable_road_lights:
            self.start_background_road_night_lights_update(reason="toggle-on")
        self.request_client_update()

    def toggle_akari_ir_bands(self) -> None:
        if not self._akari_ir_bands_toggle_supported:
            if self._action_toggle_akari_ir_bands is not None:
                self._action_toggle_akari_ir_bands.setChecked(
                    self.akari_ir_bands_opacity > 0.0
                )
            return

        enable_akari_ir_bands = self.akari_ir_bands_opacity <= 0.0
        self.akari_ir_bands_opacity = (
            self._akari_ir_bands_opacity_when_enabled
            if enable_akari_ir_bands
            else 0.0
        )
        if (
            self._action_toggle_akari_ir_bands is not None
            and self._action_toggle_akari_ir_bands.isChecked() != enable_akari_ir_bands
        ):
            self._action_toggle_akari_ir_bands.setChecked(enable_akari_ir_bands)
        self.request_client_update()

    def toggle_urban_outline(self) -> None:
        if not self._urban_outline_gui_allowed:
            if self._action_toggle_urban_outline is not None:
                self._action_toggle_urban_outline.setChecked(
                    self.urban_outline_opacity > 0.0
                )
            return
        enable_urban_outline = self.urban_outline_opacity <= 0.0
        self.urban_outline_opacity = (
            self._urban_outline_opacity_when_enabled if enable_urban_outline else 0.0
        )
        self.show_urban_outline_layer = self.urban_outline_opacity > 0.0
        if (
            self._action_toggle_urban_outline is not None
            and self._action_toggle_urban_outline.isChecked() != enable_urban_outline
        ):
            self._action_toggle_urban_outline.setChecked(enable_urban_outline)
        if enable_urban_outline:
            self.start_background_urban_outline_update(reason="toggle-on")
        self.request_client_update()

    def toggle_square_window(self) -> None:
        self.square_client_area()
        self.request_client_update()

    def toggle_fullscreen(self) -> None:
        if self.isFullScreen():
            self.showNormal()
        else:
            self.showFullScreen()
