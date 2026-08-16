"""GUI update coordinator.

``SkyWindowUpdatesMixin`` remains the public entry used by ``window.py`` and
tests. Status text, sky-data loading, cloud/geo-satellite updates, and other
overlay handlers live in sibling ``window_update_*.py`` modules.
"""

from __future__ import annotations

import random
import time
from datetime import datetime, timedelta, timezone

import numpy as np
from astropy.time import Time

from ..aircraft_constants import AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS
from ..paths import CLOUD_UPDATE_INTERVAL
from ..render.twinkle import (
    TWINKLE_TARGET_COUNT,
    nearest_twinkle_star_rows,
    sample_twinkle_direction,
    twinkle_alpha,
)
from ..satellite_constants import SATELLITE_POSITION_REFRESH_INTERVAL_SECONDS
from .window_update_cloud import SkyWindowCloudUpdatesMixin
from .window_update_common import (
    _clear_viewport_interaction_wait,
    _CloudProjectionUpdateOwner,
    _extract_sun_altitude_deg,
    _initial_data_load_active,
    _request_cloud_projection_update,
    _startup_night_light_requires_warmup,
    _ViewportInteractionState,
    _ViewportInteractionWaitOwner,
)
from .window_update_overlays import SkyWindowOverlayUpdatesMixin
from .window_update_sky import SkyWindowSkyUpdatesMixin
from .window_update_status import (
    SkyWindowStatusUpdatesMixin,
    _cloud_satellite_group,
    _cloud_source_label,
    _status_segment,
    _strip_status_prefix,
    _urban_outline_source_name,
)
_SCINTILLATION_RNG = random.Random()

# Re-export helpers that tests and callers historically imported from this module.
__all__ = [
    "SkyWindowUpdatesMixin",
    "_CloudProjectionUpdateOwner",
    "_ViewportInteractionState",
    "_ViewportInteractionWaitOwner",
    "_clear_viewport_interaction_wait",
    "_cloud_satellite_group",
    "_cloud_source_label",
    "_extract_sun_altitude_deg",
    "_initial_data_load_active",
    "_request_cloud_projection_update",
    "_startup_night_light_requires_warmup",
    "_status_segment",
    "_strip_status_prefix",
    "_urban_outline_source_name",
]


class SkyWindowUpdatesMixin(
    SkyWindowStatusUpdatesMixin,
    SkyWindowSkyUpdatesMixin,
    SkyWindowCloudUpdatesMixin,
    SkyWindowOverlayUpdatesMixin,
):
    def _background_updates_busy(self) -> bool:
        controllers = (
            self._sky_worker,
            self._cloud_controller,
            self._geosatellite_controller,
            self._satellite_controller,
            self._aircraft_controller,
            getattr(self, "_meteor_controller", None),
            self._tropical_cyclone_controller,
            self._jpl_small_body_controller,
            self._terrain_horizon_controller,
            self._water_overlay_controller,
            getattr(self, "_precipitation_controller", None),
            getattr(self, "_road_night_lights_controller", None),
            self._urban_outline_controller,
        )
        for controller in controllers:
            if controller is not None and controller.has_in_flight_update():
                return True
        return False

    def _on_scheduler_tick(self) -> None:
        if self._is_shutting_down:
            return
        now_ms = time.monotonic() * 1000.0
        now_utc = datetime.now(timezone.utc)
        if (
            self.state.jump_highlight_name is not None
            and self.state.jump_highlight_until_ms > 0.0
            and now_ms >= self.state.jump_highlight_until_ms
        ):
            self.state.jump_highlight_name = None
            self.state.jump_highlight_altaz = None
            self.state.jump_highlight_until_ms = 0.0
            self.request_client_update()
        if self._viewport_interaction_active():
            return
        current_time = self._current_time_obj()
        current_second = int(float(current_time.unix))
        if current_second % 2 == 0:
            self._publish_dynamic_display_tick(current_second)
        background_updates_busy = self._background_updates_busy()
        next_even_second = (current_second // 2 + 1) * 2
        self._request_dynamic_planet_update(
            Time(next_even_second, format="unix")
        )

        sky_next_refresh = self.state.sky_next_refresh_utc
        if not background_updates_busy and self.state.sky_update_pending:
            started = self.start_background_sky_data_update(
                star_vmag_limit=self.state.pending_star_vmag_limit,
                reason="scheduler",
                allow_during_viewport_interaction=False,
            )
            if started:
                self.state.sky_update_pending = False
                self.state.pending_star_vmag_limit = None
                self.state.sky_next_refresh_utc = now_utc + timedelta(
                    seconds=self.sky_update_interval
                )
                return

        if (
            not background_updates_busy
            and isinstance(sky_next_refresh, datetime)
            and now_utc >= sky_next_refresh
        ):
            started = self.start_background_sky_data_update(
                reason="scheduler",
                allow_during_viewport_interaction=False,
            )
            if started:
                self.state.sky_next_refresh_utc = now_utc + timedelta(
                    seconds=self.sky_update_interval
                )
                return

        sky_disc_next_refresh = self.state.sky_disc_next_refresh_utc
        if (
            not background_updates_busy
            and isinstance(sky_disc_next_refresh, datetime)
            and now_utc >= sky_disc_next_refresh
        ):
            started = self.start_background_sky_disc_update(reason="scheduler")
            if started:
                self.state.sky_disc_next_refresh_utc = now_utc + timedelta(
                    seconds=self._sky_disc_update_interval()
                )
                return

        persistent_next_refresh = self.state.persistent_search_next_refresh_utc
        if (
            not background_updates_busy
            and isinstance(persistent_next_refresh, datetime)
            and now_utc >= persistent_next_refresh
        ):
            started = self._start_persistent_search_refresh(reason="scheduler")
            if started:
                self.state.persistent_search_next_refresh_utc = (
                    persistent_next_refresh + timedelta(hours=1)
                )
                return

        cloud_next_refresh = self.state.cloud_next_refresh_utc
        if (
            not background_updates_busy
            and isinstance(cloud_next_refresh, datetime)
            and now_utc >= cloud_next_refresh
            and self._cloud_layer_enabled()
        ):
            controller = (
                self._geosatellite_controller
                if self._geo_satellite_mode_active()
                else self._cloud_controller
            )
            if controller is None or controller.has_in_flight_update():
                return
            self.start_background_cloud_update(reason="scheduler")
            if controller.has_in_flight_update():
                self.state.cloud_next_refresh_utc = now_utc + timedelta(
                    seconds=CLOUD_UPDATE_INTERVAL
                )
                return

        precipitation_next_refresh = getattr(
            self.state, "precipitation_next_refresh_utc", None
        )
        if (
            not background_updates_busy
            and self._precipitation_layer_enabled()
            and isinstance(precipitation_next_refresh, datetime)
            and now_utc >= precipitation_next_refresh
        ):
            if self.start_background_precipitation_update(reason="scheduler"):
                return

        satellite_next_refresh = self.state.satellite_next_refresh_utc
        if (
            not background_updates_busy
            and self._satellite_layer_enabled()
            and isinstance(satellite_next_refresh, datetime)
            and now_utc >= satellite_next_refresh
        ):
            started = self.start_background_satellite_update(reason="scheduler")
            if started:
                self.state.satellite_next_refresh_utc = now_utc + timedelta(
                    seconds=SATELLITE_POSITION_REFRESH_INTERVAL_SECONDS
                )
                return

        aircraft_next_refresh = self.state.aircraft_next_refresh_utc
        if (
            not background_updates_busy
            and self._aircraft_layer_enabled()
            and isinstance(aircraft_next_refresh, datetime)
            and now_utc >= aircraft_next_refresh
        ):
            started = self.start_background_aircraft_update(reason="scheduler")
            if started:
                self.state.aircraft_next_refresh_utc = now_utc + timedelta(
                    seconds=AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS
                )
                return

        meteor_next_refresh = getattr(self.state, "meteor_next_refresh_utc", None)
        if (
            not background_updates_busy
            and float(getattr(self, "meteor_opacity", 0.0)) > 0.0
            and isinstance(meteor_next_refresh, datetime)
            and now_utc >= meteor_next_refresh
        ):
            if self.start_background_meteor_update(reason="scheduler"):
                self.state.meteor_next_refresh_utc = now_utc + timedelta(hours=1)
                return

        cyclone_state = self.tropical_cyclone_state
        cyclone_next_check = cyclone_state.next_check_utc
        cyclone_next_refresh = cyclone_state.next_refresh_utc
        if (
            not background_updates_busy
            and self._tropical_cyclone_layer_enabled()
            and (
                (
                    isinstance(cyclone_next_refresh, datetime)
                    and now_utc >= cyclone_next_refresh
                )
                or (
                    isinstance(cyclone_next_check, datetime)
                    and now_utc >= cyclone_next_check
                )
            )
        ):
            controller = self._tropical_cyclone_controller
            if controller is not None and not controller.has_in_flight_update():
                started = self.start_background_tropical_cyclone_update(
                    reason="scheduler"
                )
                if started:
                    return

        # Lowest-priority idle work for overlays that are not part of the display tick.
        if (
            self._tropical_cyclone_layer_enabled()
            and self._tropical_cyclone_projection_next_refresh_delay_ms() == 0
        ):
            self.reproject_tropical_cyclone_overlay()
            return

        if (
            self._cloud_layer_enabled()
            and self._cloud_projection_next_refresh_delay_ms() == 0
        ):
            self._start_cloud_projection_update(reason="scheduler")
            return

    def _publish_dynamic_display_tick(self, current_second: int) -> None:
        state = self.state
        if state.dynamic_display_second == current_second:
            return
        state.dynamic_display_second = current_second
        state.dynamic_display_time = Time(current_second, format="unix")
        self._prepare_star_interpolation_mesh_cache()
        planet_bucket = current_second // 2
        if state.prepared_dynamic_planet_bucket == planet_bucket:
            state.dynamic_planets = state.prepared_dynamic_planets
            state.dynamic_planet_bucket = planet_bucket
            state.prepared_dynamic_planets = None
            state.prepared_dynamic_planet_bucket = None
        self._update_twinkle(time_obj=state.dynamic_display_time, request_repaint=False)
        now_utc = datetime.now(timezone.utc)
        if self._aircraft_layer_enabled():
            state.aircraft_projection_next_refresh_utc = now_utc + timedelta(
                seconds=AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS
            )
        if self._satellite_layer_enabled():
            state.satellite_projection_next_refresh_utc = now_utc + timedelta(
                seconds=SATELLITE_POSITION_REFRESH_INTERVAL_SECONDS
            )
        self.request_client_update()

    def _update_twinkle(
        self, *, time_obj: Time | None = None, request_repaint: bool = True
    ) -> None:
        """Choose the transient faint-star dimming target for this 2-second bucket."""
        state = self.state
        if (
            str(getattr(self, "presentation_id", "scenic")).strip().lower() != "scenic"
            or not bool(getattr(self, "twinkle_enabled", True))
            or bool(getattr(state, "simplified_view_enabled", False))
            or bool(getattr(state, "viewport_interaction_mode", False))
            or state.celestial_data is None
        ):
            state.twinkle_targets = ()
            return
        try:
            effective_time = time_obj or self._current_time_obj()
            time_bucket = int(float(effective_time.unix) // 2.0)
        except Exception:
            state.twinkle_targets = ()
            return
        if state.twinkle_bucket == time_bucket:
            return
        state.twinkle_bucket = time_bucket
        state.twinkle_targets = ()
        viewer = self._viewer_data_for_render()
        stars = state.celestial_data.stars
        selected_rows: set[int] = set()
        selected_targets: list[tuple[int, float]] = []
        twinkle_count = max(
            0,
            int(getattr(self, "twinkle_count", TWINKLE_TARGET_COUNT)),
        )
        directions = [
            sample_twinkle_direction(
                rng=_SCINTILLATION_RNG,
            )
            for _ in range(twinkle_count)
        ]
        if directions:
            target_alt, target_az = np.asarray(directions, dtype=float).T
            nearest_rows = nearest_twinkle_star_rows(
                stars,
                target_alt_deg=target_alt,
                target_az_deg=target_az,
                view_center=viewer.view_center,
                content_fov_deg=viewer.edge_fov_deg,
                vmag_limit=float(self.vmag_limit),
            )
        else:
            nearest_rows = np.array([], dtype=np.int32)
        for row in nearest_rows:
            star_row = int(row)
            if star_row < 0 or star_row in selected_rows:
                continue
            selected_rows.add(star_row)
            selected_targets.append(
                (star_row, twinkle_alpha(float(stars["alt"][star_row])))
            )
        state.twinkle_targets = tuple(selected_targets)
        if request_repaint:
            self.request_client_update()

    def _on_periodic_debug_snapshot_timer(self) -> None:
        """Queue the current frame at the periodic debug interval."""
        if self._resolve_periodic_debug_snapshot_dir() is None:
            return
        self._queue_periodic_debug_snapshot(
            {
                "refreshed_at_utc": datetime.now(timezone.utc),
                "source": "periodic",
            }
        )
        self.request_client_update()
