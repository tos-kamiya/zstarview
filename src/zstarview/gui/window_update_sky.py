from __future__ import annotations

import logging
from datetime import datetime, timedelta, timezone
from typing import cast

import astropy.time

from ..render import geometry as render_geometry
from ..render import sky_disc as render_sky_disc
from .window_update_common import (
    _clear_viewport_interaction_wait,
    _CloudProjectionUpdateOwner,
    _initial_data_load_active,
    _request_cloud_projection_update,
    _startup_night_light_requires_warmup,
    _ViewportInteractionWaitOwner,
)
from .window_update_overlays import precipitation_layer_enabled

logger = logging.getLogger(__name__)


class SkyWindowSkyUpdatesMixin:
    def _viewport_interaction_active(self) -> bool:
        return bool(self.state.viewport_interaction_mode)

    def _on_sky_disc_calculated(self, payload: dict) -> None:
        current_generation = int(getattr(self, "_disc_generation", 0))
        current_geometry = render_geometry.get_screen_geometry(
            max(2, int(self.client_width())),
            max(2, int(self.client_height())),
            self.viewer_data.view_alt_deg,
            edge_fov_deg=self.viewer_data.edge_fov_deg,
            content_fov_deg=self.viewer_data.content_fov_deg,
        )
        if (
            int(payload.get("render_generation", current_generation))
            != current_generation
            or payload.get("geometry") != current_geometry
            or tuple(payload.get("view_center", ()))
            != tuple(self.viewer_data.view_center)
        ):
            return
        sky_disc_image = payload.get("sky_disc")
        if sky_disc_image is None:
            return
        self.state.sky_disc_image = sky_disc_image
        self.state.sky_disc_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=render_sky_disc.sky_disc_update_interval(payload.get("sun_alt_deg"))
        )
        self._compositor.invalidate()
        self.request_client_update()

    def _on_sky_data_calculated(self, payload: dict) -> None:
        current_generation = int(getattr(self, "_disc_generation", 0))
        payload_generation = int(payload.get("render_generation", current_generation))
        payload_geometry = payload.get("geometry")
        current_geometry = render_geometry.get_screen_geometry(
            max(2, int(self.client_width())),
            max(2, int(self.client_height())),
            self.viewer_data.view_alt_deg,
            edge_fov_deg=self.viewer_data.edge_fov_deg,
            content_fov_deg=self.viewer_data.content_fov_deg,
        )
        if (
            payload_generation != current_generation
            or payload_geometry != current_geometry
        ):
            logger.debug(
                "Discard stale sky payload generation=%s current=%s geometry=%s current_geometry=%s",
                payload_generation,
                current_generation,
                payload_geometry,
                current_geometry,
            )
            if not self._is_shutting_down:
                if self.state.viewport_interaction_release_pending:
                    _clear_viewport_interaction_wait(
                        cast(_ViewportInteractionWaitOwner, self)
                    )
                self.request_sky_data_update(
                    reason="stale-render",
                )
            return
        view_center = payload.get("view_center", self.viewer_data.view_center)
        current_view_center = tuple(self.viewer_data.view_center)
        payload_matches_current_view = (
            isinstance(view_center, (tuple, list))
            and len(view_center) >= 2
            and abs(float(view_center[0]) - float(current_view_center[0])) < 1e-9
            and abs(float(view_center[1]) - float(current_view_center[1])) < 1e-9
        )
        if not payload_matches_current_view and (
            not self.state.viewport_interaction_mode
            or self.state.viewport_interaction_release_pending
        ):
            logger.debug(
                "Discard stale sky payload view_center=%s current=%s generation=%s",
                view_center,
                current_view_center,
                payload_generation,
            )
            if not self._is_shutting_down:
                if self.state.viewport_interaction_release_pending:
                    _clear_viewport_interaction_wait(
                        cast(_ViewportInteractionWaitOwner, self)
                    )
                self.request_sky_data_update(
                    reason="stale-view-center",
                )
            return
        if not self.state.viewport_interaction_mode:
            if isinstance(view_center, (tuple, list)) and len(view_center) >= 2:
                self.state.render_view_center = (
                    float(view_center[0]),
                    float(view_center[1]),
                )
            else:
                self.state.render_view_center = tuple(self.viewer_data.view_center)
        self.state.celestial_data = payload["celestial"]
        self.state.twinkle_bucket = None
        self.state.twinkle_targets = ()
        self.state.dynamic_planets = None
        self.state.dynamic_planet_bucket = None
        self.state.prepared_dynamic_planets = None
        self.state.prepared_dynamic_planet_bucket = None
        self.state.dynamic_planet_requested_bucket = None
        self.state.sky_disc_image = payload["sky_disc"]
        self.state.night_light_glow_profile = payload.get("night_light_glow_profile")
        self.state.sky_disc_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=self._sky_disc_update_interval()
        )

        if self.state.viewport_interaction_release_pending:
            refresh_reason = (
                self.state.viewport_interaction_completion_reason
                or "view-change-release"
            )
            self.state.viewport_interaction_release_pending = False
            self.state.viewport_interaction_completion_reason = None
            self.state.viewport_interaction_mode = False
            self.state.viewport_interaction_stars = None
            self._sync_viewport_interaction_chrome_visibility()
            if not self._is_shutting_down:
                _request_cloud_projection_update(
                    cast(_CloudProjectionUpdateOwner, self), reason=refresh_reason
                )
                self.start_background_terrain_horizon_update(reason=refresh_reason)
                self.reproject_tropical_cyclone_overlay()

        self._compositor.invalidate()
        self.request_client_update()

        if _initial_data_load_active(self):
            self._startup_initial_sky_loaded = True
            if (
                not _startup_night_light_requires_warmup(self, payload)
                or self.state.night_light_glow_profile is not None
            ):
                self._startup_initial_night_light_loaded = True
            self._continue_initial_data_load()
            return

        self.state.sky_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=self.sky_update_interval
        )

        if self.state.sky_update_pending and not self._is_shutting_down:
            self.request_sky_data_update(
                self.state.pending_star_vmag_limit,
                reason="pending-follow-up",
            )

        if self.state.cloud_repaint_deferred and not self.state.interaction_mode:
            self.state.cloud_repaint_deferred = False
            self._safe_request_cloud_repaint()

    def _continue_initial_data_load(self) -> None:
        if self._is_shutting_down:
            return
        if not _initial_data_load_active(self):
            return
        if not self._startup_initial_sky_loaded:
            return
        if hasattr(self, "start_background_road_night_lights_update"):
            self.start_background_road_night_lights_update(reason="initial")
        if precipitation_layer_enabled(self):
            self.start_background_precipitation_update(reason="initial")
        if (
            self.terrain_horizon_opacity > 0.0
            and not self._startup_initial_terrain_loaded
        ):
            if self.start_background_terrain_horizon_update(reason="initial"):
                return
            return
        if self.urban_outline_opacity > 0.0 and not self._startup_initial_urban_loaded:
            if self.start_background_urban_outline_update(reason="initial"):
                return
            return
        self._finish_initial_data_load()

    def _finish_initial_data_load(self) -> None:
        if self._startup_initial_data_loaded:
            return
        self._startup_initial_data_loaded = True
        self.initial_data_loaded.emit()
        if getattr(self, "_startup_resize_pending", False):
            self._startup_resize_pending = False
            self.request_sky_data_update(reason="startup-resize")

    def request_sky_data_update(
        self,
        star_vmag_limit: float | None = None,
        *,
        reason: str = "manual",
        allow_during_viewport_interaction: bool = False,
    ) -> bool:
        if (
            self._viewport_interaction_active()
            and not allow_during_viewport_interaction
        ):
            self.state.sky_update_pending = True
            self.state.pending_star_vmag_limit = star_vmag_limit
            logger.debug(
                "Sky data update deferred during viewport interaction (reason=%s).",
                reason,
            )
            return False
        if self.start_background_sky_data_update(
            star_vmag_limit=star_vmag_limit,
            reason=reason,
            allow_during_viewport_interaction=allow_during_viewport_interaction,
        ):
            self.state.sky_update_pending = False
            self.state.pending_star_vmag_limit = None
            return True
        self.state.sky_update_pending = True
        self.state.pending_star_vmag_limit = star_vmag_limit
        logger.debug(
            "Sky data update deferred; worker is busy (reason=%s).",
            reason,
        )
        return False

    def _request_dynamic_planet_update(
        self, target_time: astropy.time.Time | None = None
    ) -> None:
        if str(getattr(self, "presentation_id", "scenic")).strip().lower() != "scenic":
            return
        if self.state.celestial_data is None or self._viewport_interaction_active():
            return
        if target_time is None:
            current_unix = float(self._current_time_obj().unix)
            target_unix = (int(current_unix) // 2 + 1) * 2
            target_time = astropy.time.Time(target_unix, format="unix")
        bucket = int(float(target_time.unix) // 2.0)
        if bucket in {
            self.state.dynamic_planet_bucket,
            self.state.prepared_dynamic_planet_bucket,
            self.state.dynamic_planet_requested_bucket,
        }:
            return
        ephemeris = self._ephemeris
        if ephemeris is None:
            ephemeris = self._services.ephemeris.load()
        if self._sky_worker.update_planets(
            ephemeris=ephemeris,
            viewer_data=self._viewer_data_for_render(),
            time_obj=target_time,
            render_generation=int(getattr(self, "_disc_generation", 0)),
        ):
            self.state.dynamic_planet_requested_bucket = bucket

    def _on_planet_data_calculated(self, payload: dict) -> None:
        if self._is_shutting_down:
            return
        planets = payload.get("planets")
        if not isinstance(planets, list):
            return
        bucket = int(float(payload.get("time_unix", 0.0)) // 2.0)
        if int(payload.get("render_generation", -1)) != int(
            getattr(self, "_disc_generation", 0)
        ):
            if self.state.dynamic_planet_requested_bucket == bucket:
                self.state.dynamic_planet_requested_bucket = None
            return
        self.state.prepared_dynamic_planets = planets
        self.state.prepared_dynamic_planet_bucket = bucket
        if self.state.dynamic_planet_requested_bucket == bucket:
            self.state.dynamic_planet_requested_bucket = None

    def _current_sun_alt_deg(self) -> float | None:
        planets = self.state.dynamic_planets
        if planets is None and self.state.celestial_data is not None:
            planets = self.state.celestial_data.planets
        if planets is None:
            return None
        for body in planets:
            if body.name == "sun":
                try:
                    return float(body.alt)
                except (TypeError, ValueError):
                    return None
        return None

    def _sky_disc_update_interval(self) -> int:
        return render_sky_disc.sky_disc_update_interval(self._current_sun_alt_deg())

    def start_background_sky_disc_update(self, reason: str = "manual") -> bool:
        del reason
        if self._viewport_interaction_active():
            return False
        ephemeris = self._ephemeris
        if ephemeris is None:
            ephemeris = self._services.ephemeris.load()
        return self._sky_worker.update_sky_disc(
            ephemeris=ephemeris,
            viewer_data=self.viewer_data,
            geometry=render_geometry.get_screen_geometry(
                max(2, int(self.client_width())),
                max(2, int(self.client_height())),
                self.viewer_data.view_alt_deg,
                edge_fov_deg=self.viewer_data.edge_fov_deg,
                content_fov_deg=self.viewer_data.content_fov_deg,
            ),
            delta_t=self.delta_t,
            sky_disc_alpha=self.sky_disc_alpha,
            theme=self.theme,
            image_size=(
                max(2, int(self.client_width())),
                max(2, int(self.client_height())),
            ),
            render_generation=int(self._disc_generation),
        )

    def start_background_sky_data_update(
        self,
        is_initial_load: bool = False,
        star_vmag_limit: float | None = None,
        reason: str = "manual",
        allow_during_viewport_interaction: bool = False,
    ) -> bool:
        if (
            self._viewport_interaction_active()
            and not allow_during_viewport_interaction
        ):
            return False
        _, _ = self.viewer_data.location
        use_lod6_catalog = star_vmag_limit is not None and float(star_vmag_limit) <= 6.0
        star_catalog = self.star_catalog_np
        star_subset_indices = (
            self.star_catalog_lod6_indices if use_lod6_catalog else None
        )
        worker_star_vmag_limit = None if use_lod6_catalog else star_vmag_limit
        ephemeris = self._ephemeris
        if ephemeris is None:
            ephemeris = self._services.ephemeris.load()
        started = self._sky_worker.update(
            ephemeris=ephemeris,
            viewer_data=self.viewer_data,
            geometry=render_geometry.get_screen_geometry(
                max(2, int(self.client_width())),
                max(2, int(self.client_height())),
                self.viewer_data.view_alt_deg,
                edge_fov_deg=self.viewer_data.edge_fov_deg,
                content_fov_deg=self.viewer_data.content_fov_deg,
            ),
            star_catalog=star_catalog,
            dso_catalog=self.dso_catalog_np,
            star_vmag_limit=worker_star_vmag_limit,
            star_subset_indices=star_subset_indices,
            delta_t=self.delta_t,
            sky_update_interval=self.sky_update_interval,
            sky_disc_alpha=self.sky_disc_alpha,
            theme=self.theme,
            star_catalog_meta=self.star_catalog_meta,
            image_size=(
                max(2, int(self.client_width())),
                max(2, int(self.client_height())),
            ),
            terrain_horizon_profile_altaz=self.terrain_horizon_state.profile_altaz,
            terrain_horizon_profile_distances_m=self.terrain_horizon_state.profile_distances_m,
            terrain_secondary_ridges_altaz_layers=self.terrain_horizon_state.secondary_ridges_altaz_layers,
            terrain_secondary_ridges_distances_m_layers=self.terrain_horizon_state.secondary_ridges_distances_m_layers,
            terrain_sample_distances_m=self.terrain_horizon_state.sample_distances_m,
            terrain_sample_terrain_elevation_m=self.terrain_horizon_state.sample_terrain_elevation_m,
            night_light_glow_profile=self.state.night_light_glow_profile,
            night_light_opacity=float(getattr(self, "night_light_opacity", 0.0)),
            render_generation=int(self._disc_generation),
        )
        if started:
            if is_initial_load:
                logger.info("Calculating initial sky data (reason=%s)...", reason)
            else:
                logger.info("Updating sky data (reason=%s)...", reason)
        return started
