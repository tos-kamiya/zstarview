from __future__ import annotations

from dataclasses import replace

from PySide6.QtCore import QEvent, Qt
from PySide6.QtGui import QKeyEvent, QMouseEvent
from PySide6.QtWidgets import QApplication

from ..types import ViewerData
from .view_direction import clamp_view_center_alt_az, resolve_view_direction_step

HOVER_REPAINT_INTERVAL_MS = 33


def _replace_viewer_data(viewer_data: ViewerData, /, **changes: object) -> ViewerData:
    return replace(viewer_data, **changes)

class SkyWindowInputMixin:
    def eventFilter(self, watched: object, event: QEvent) -> bool:
        if isinstance(watched, QApplication):
            if event.type() in (QEvent.Type.KeyPress, QEvent.Type.KeyRelease):
                if self._startup_input_blocked():
                    return True
                key_event = event if isinstance(event, QKeyEvent) else None
                if key_event is not None and key_event.key() in {
                    Qt.Key.Key_Left,
                    Qt.Key.Key_Right,
                    Qt.Key.Key_Up,
                    Qt.Key.Key_Down,
                }:
                    if (
                        QApplication.activePopupWidget() is None
                        and self.isActiveWindow()
                    ):
                        if event.type() == QEvent.Type.KeyPress:
                            self._handle_client_key_press(key_event)
                        else:
                            self._handle_client_key_release(key_event)
                        return True
        return super().eventFilter(watched, event)
    def _begin_interaction_mode(self) -> None:
        self.state.interaction_mode = True
        self._interaction_idle_timer.start()

    def _end_interaction_mode(self) -> None:
        if not self.state.interaction_mode:
            return
        self.state.interaction_mode = False
        self.request_sky_data_update(reason="interaction-idle")
        self.request_cloud_projection_update(reason="view-change-idle")
        self.start_background_terrain_horizon_update(reason="view-change-idle")

    def _begin_viewport_interaction_mode(
        self,
        preserve_cloud_buffers: bool = False,
        *,
        start_idle_timer: bool = True,
    ) -> None:
        if not self._startup_initial_load_started:
            return
        self.state.viewport_interaction_mode = True
        self.state.viewport_interaction_release_pending = False
        cloud_controller = self._cloud_controller
        from . import window as window_module

        # Geo-satellite projection is camera-independent after the source image
        # has been converted to an alt/az grid.  Keep that grid available while
        # the user is rotating the view; otherwise a slow or failed follow-up
        # request leaves the viewport with an empty cloud layer.
        preserve_geo_satellite_buffers = False
        geo_mode_active = getattr(self, "_geo_satellite_mode_active", None)
        if callable(geo_mode_active):
            preserve_geo_satellite_buffers = bool(geo_mode_active())
        cleared_cloud = window_module.SkyWindow._clear_cloud_render_buffers(
            self,
            preserve_cloud_buffers=preserve_cloud_buffers,
            preserve_geo_satellite_buffers=preserve_geo_satellite_buffers,
        )
        if cleared_cloud:
            self._compositor.invalidate()
        if cloud_controller is not None:
            invalidate = cloud_controller.invalidate_pending_render_results
            if callable(invalidate):
                invalidate()
        window_module.SkyWindow._sync_viewport_interaction_chrome_visibility(self)
        if start_idle_timer:
            timer = self._viewport_interaction_idle_timer
            if timer.isActive():
                timer.stop()
            timer.start()

    def _update_viewport_interaction_stars(self) -> None:
        if self.state.celestial_data is None:
            self.state.viewport_interaction_stars = None
            return
        from . import window as window_module

        stars, _location = window_module.calculate_visible_stars(
            self.star_catalog_np,
            self.viewer_data.location[0],
            self.viewer_data.location[1],
            self.viewer_data.observer_height_m,
            self._current_time_obj(),
            self.state.render_view_center,
            max_vmag=4.0,
            subset_indices=self.star_catalog_lod6_indices,
            star_data_policy=getattr(
                self,
                "star_data_policy",
                "scenic_view_scoped",
            ),
        )
        self.state.viewport_interaction_stars = stars

    def _end_viewport_interaction_mode(
        self,
        reason: str = "viewport-interaction-idle",
    ) -> None:
        from . import window as window_module

        if not self._startup_initial_load_started:
            return
        if (
            reason == "viewport-interaction-idle"
            and not getattr(self, "_startup_initial_data_loaded", True)
            and bool(getattr(self, "_startup_input_blocked_state", False))
        ):
            return
        if not self.state.viewport_interaction_mode:
            return
        sky_update_started = self.request_sky_data_update(
            reason=reason,
            allow_during_viewport_interaction=True,
        )
        if reason.endswith("release"):
            if not sky_update_started:
                self.state.viewport_interaction_mode = False
                self.state.viewport_interaction_stars = None
                self.state.viewport_interaction_completion_reason = None
                window_module.SkyWindow._sync_viewport_interaction_chrome_visibility(
                    self
                )
                self.reproject_tropical_cyclone_overlay()
                self.request_client_update()
                return
            self.state.viewport_interaction_release_pending = True
            self.state.viewport_interaction_completion_reason = "view-change-release"
            self.reproject_tropical_cyclone_overlay(
                allow_during_viewport_interaction=True
            )
            return
        self.state.viewport_interaction_mode = False
        self.state.viewport_interaction_stars = None
        self.state.viewport_interaction_completion_reason = None
        window_module.SkyWindow._sync_viewport_interaction_chrome_visibility(self)
        refresh_reason = (
            "view-change-release" if reason.endswith("release") else "view-change-idle"
        )
        self.request_cloud_projection_update(reason=refresh_reason)
        self.start_background_terrain_horizon_update(reason=refresh_reason)
        self.reproject_tropical_cyclone_overlay()
        self.request_client_update()
    def _handle_client_leave(self, event: QEvent) -> None:
        self.state.mouse_pos = None
        if self._hover_repaint_timer.isActive():
            self._hover_repaint_timer.stop()
        self.request_client_update()
        event.accept()

    def _handle_client_mouse_move(self, event: QMouseEvent) -> None:
        if self._startup_input_blocked():
            event.accept()
            return
        self.state.mouse_pos = event.pos()
        if not self._hover_repaint_timer.isActive():
            self._hover_repaint_timer.start()
        event.accept()

    def _rotate_view(
        self,
        d_alt: float = 0.0,
        d_az: float = 0.0,
        *,
        interactive_viewport: bool = False,
        start_viewport_idle_timer: bool = True,
    ) -> None:
        alt, az = self.viewer_data.view_center
        from . import window as window_module

        window_module.SkyWindow._set_view_center(
            self,
            alt + d_alt,
            az + d_az,
            interactive_viewport=interactive_viewport,
            start_viewport_idle_timer=start_viewport_idle_timer,
        )

    def _set_view_center(
        self,
        alt_deg: float,
        az_deg: float,
        *,
        interactive_viewport: bool = False,
        start_viewport_idle_timer: bool = True,
    ) -> None:
        if interactive_viewport:
            if not bool(self.state.viewport_interaction_mode):
                self._begin_viewport_interaction_mode(
                    start_idle_timer=start_viewport_idle_timer
                )
        new_alt, new_az = clamp_view_center_alt_az(alt_deg, az_deg)
        self.viewer_data = _replace_viewer_data(
            self.viewer_data, view_center=(new_alt, new_az)
        )
        self.state.render_view_center = (new_alt, new_az)
        self._sync_view_altitude_actions()
        if interactive_viewport:
            self._update_viewport_interaction_stars()
            self.request_client_update()
            return
        if self.state.viewport_interaction_mode:
            self._end_viewport_interaction_mode(reason="view-change-idle")
            return
        self._begin_interaction_mode()
        self.request_sky_data_update()
        self.request_client_update()

    def _viewport_rotation_keys(self) -> set[int]:
        return self._viewport_rotation_keys_down

    def _handle_client_key_press(self, event: QKeyEvent) -> None:
        if self._startup_input_blocked():
            event.accept()
            return
        key = event.key()
        modifiers = event.modifiers()

        # --- View Control ---
        if key == Qt.Key.Key_Left:
            step_deg = resolve_view_direction_step(modifiers, self.state.rotation_step)
            if not event.isAutoRepeat():
                self._viewport_rotation_keys().add(key)
            self._rotate_view(
                d_az=-step_deg,
                interactive_viewport=True,
                start_viewport_idle_timer=False,
            )
            event.accept()
        elif key == Qt.Key.Key_Right:
            step_deg = resolve_view_direction_step(modifiers, self.state.rotation_step)
            if not event.isAutoRepeat():
                self._viewport_rotation_keys().add(key)
            self._rotate_view(
                d_az=step_deg,
                interactive_viewport=True,
                start_viewport_idle_timer=False,
            )
            event.accept()
        elif key == Qt.Key.Key_Up:
            step_deg = resolve_view_direction_step(modifiers, self.state.rotation_step)
            if not event.isAutoRepeat():
                self._viewport_rotation_keys().add(key)
            self._rotate_view(
                d_alt=step_deg,
                interactive_viewport=True,
                start_viewport_idle_timer=False,
            )
            event.accept()
        elif key == Qt.Key.Key_Down:
            step_deg = resolve_view_direction_step(modifiers, self.state.rotation_step)
            if not event.isAutoRepeat():
                self._viewport_rotation_keys().add(key)
            self._rotate_view(
                d_alt=-step_deg,
                interactive_viewport=True,
                start_viewport_idle_timer=False,
            )
            event.accept()

        # --- Toggles ---
        elif key == Qt.Key.Key_Space:
            if not event.isAutoRepeat():
                self.toggle_simplified_view()
            event.accept()
        elif key == Qt.Key.Key_M:
            self.toggle_moon_option()
            event.accept()
        elif key == Qt.Key.Key_C:
            self.toggle_clouds()
            event.accept()
        elif key == Qt.Key.Key_A:
            self.toggle_aircraft()
            event.accept()
        elif key == Qt.Key.Key_T:
            self.toggle_terrain_horizon()
            event.accept()
        elif key == Qt.Key.Key_W:
            self.toggle_water_overlay()
            event.accept()
        elif key == Qt.Key.Key_E:
            self.toggle_earth_guide()
            event.accept()
        elif key == Qt.Key.Key_L:
            self.toggle_night_lights()
            event.accept()
        elif key == Qt.Key.Key_U:
            self.toggle_urban_outline()
            event.accept()
        elif key == Qt.Key.Key_D:
            self.toggle_dso()
            event.accept()
        elif key == Qt.Key.Key_P:
            self.toggle_asterisms()
            event.accept()
        elif key == Qt.Key.Key_G:
            self.toggle_guidelines()
            event.accept()
        elif key == Qt.Key.Key_S:
            self.toggle_sky_disc()
            event.accept()

        # --- Window Control ---
        elif key == Qt.Key.Key_F11:
            self.toggle_fullscreen()
            event.accept()
        elif key == Qt.Key.Key_Escape:
            if self.isFullScreen():
                self.showNormal()
            event.accept()
        elif key == Qt.Key.Key_Q:
            self._request_application_quit()
            event.accept()
        else:
            super().keyPressEvent(event)

    def _handle_client_key_release(self, event: QKeyEvent) -> None:
        if self._startup_input_blocked():
            event.accept()
            return
        key = event.key()
        if key not in {
            Qt.Key.Key_Left,
            Qt.Key.Key_Right,
            Qt.Key.Key_Up,
            Qt.Key.Key_Down,
        }:
            if key == Qt.Key.Key_Space:
                event.accept()
            return
        if event.isAutoRepeat():
            event.accept()
            return
        keys_down = self._viewport_rotation_keys()
        keys_down.discard(key)
        if not keys_down:
            self._end_viewport_interaction_mode(reason="viewport-interaction-release")
        event.accept()
