from __future__ import annotations

import logging
from dataclasses import replace
from datetime import datetime, timedelta, timezone

from ..aircraft_constants import AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS
from ..overlay_time import overlay_availability_for_delta
from ..precipitation import PRECIPITATION_REFRESH_SECONDS
from ..satellite_constants import SATELLITE_POSITION_REFRESH_INTERVAL_SECONDS
from ..search.jpl import project_jpl_target_altaz_from_state_vector
from ..tropical_cyclones.cache import (
    TROPICAL_CYCLONE_CACHE_TTL_SECONDS,
    TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS,
)
from .window_update_common import _initial_data_load_active

logger = logging.getLogger(__name__)


def precipitation_layer_enabled(obj: object) -> bool:
    return bool(
        overlay_availability_for_delta(
            getattr(obj, "delta_t", timedelta(0))
        ).precipitation
        and float(getattr(obj, "precipitation_opacity", 0.0)) > 0.0
        and getattr(obj, "_precipitation_controller", None) is not None
    )


class SkyWindowOverlayUpdatesMixin:
    def _water_overlay_action_enabled(self) -> bool:
        if not bool(self._water_overlay_gui_allowed):
            return False
        return float(self.terrain_horizon_opacity) > 0.0

    def _sync_water_overlay_action_enabled(self) -> None:
        action = self._action_toggle_water_overlay
        if action is not None:
            action.setEnabled(self._water_overlay_action_enabled())

    def _tropical_cyclone_layer_enabled(self) -> bool:
        return bool(
            self.show_tropical_cyclone_overlay
            and float(self.tropical_cyclone_opacity) > 0.0
            and self._tropical_cyclone_controller is not None
        )

    def _precipitation_layer_enabled(self) -> bool:
        return precipitation_layer_enabled(self)

    def _satellite_projection_next_refresh_delay_ms(self) -> int | None:
        next_refresh_utc = self.state.satellite_projection_next_refresh_utc
        if next_refresh_utc is None:
            return None
        return max(
            0,
            int(
                round(
                    (next_refresh_utc - datetime.now(timezone.utc)).total_seconds()
                    * 1000.0
                )
            ),
        )

    def _aircraft_projection_next_refresh_delay_ms(self) -> int | None:
        next_refresh_utc = self.state.aircraft_projection_next_refresh_utc
        if next_refresh_utc is None:
            return None
        return max(
            0,
            int(
                round(
                    (next_refresh_utc - datetime.now(timezone.utc)).total_seconds()
                    * 1000.0
                )
            ),
        )

    def _tropical_cyclone_projection_next_refresh_delay_ms(self) -> int | None:
        next_refresh_utc = self.state.tropical_cyclone_projection_next_refresh_utc
        if next_refresh_utc is None:
            return None
        return max(
            0,
            int(
                round(
                    (next_refresh_utc - datetime.now(timezone.utc)).total_seconds()
                    * 1000.0
                )
            ),
        )

    def start_background_meteor_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down or float(self.meteor_opacity) <= 0.0:
            return False
        controller = self._meteor_controller
        if controller is None:
            return False
        lat, lon = self.viewer_data.location
        return controller.update(
            display_time_utc=self._current_time_obj().to_datetime(timezone=timezone.utc),
            observer_lat=lat,
            observer_lon=lon,
            observer_height_m=self.viewer_data.observer_height_m,
            max_display_trails=self.meteor_trails_max_candidates,
            reason=reason,
        )

    def _on_meteor_started(self, payload: dict) -> None:
        self.meteor_state.set_banner(str(payload.get("banner", "GMN meteors: loading...")))
        self.request_client_update()

    def _on_meteor_ready(self, payload: dict) -> None:
        result = payload.get("result")
        if result is not None:
            self.meteor_state.set_result(result)
            self.state.meteor_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
                hours=1
            )
        self.request_client_update()

    def _on_meteor_failed(self, payload: dict) -> None:
        self.meteor_state.set_banner(str(payload.get("banner", "GMN meteors: unavailable")))
        self.state.meteor_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            minutes=10
        )
        self.request_client_update()

    def start_background_satellite_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down:
            return False
        if self._viewport_interaction_active():
            return False
        if float(self.satellite_opacity) <= 0.0:
            return False
        if self._satellite_controller is None:
            return False
        lat, lon = self.viewer_data.location
        return self._satellite_controller.update(
            observer_lat=lat,
            observer_lon=lon,
            observer_height_m=self.viewer_data.observer_height_m,
            time_obj=self._current_time_obj(),
            enabled_groups=tuple(self._enabled_satellite_groups),
            reason=reason,
        )

    def reproject_satellite_overlay(self) -> None:
        if float(self.satellite_opacity) <= 0.0:
            return
        if self._viewport_interaction_active():
            return
        validity_remaining_ms = self._satellite_validity_remaining_ms()
        if validity_remaining_ms is not None and validity_remaining_ms <= 0:
            self.state.satellite_projection_next_refresh_utc = None
            self.request_client_update()
            self.start_background_satellite_update(reason="time-window-shift")
            return
        records_by_group = self.satellite_state.records_by_group or {}
        if not records_by_group:
            self.state.satellite_projection_next_refresh_utc = None
            self.request_client_update()
            return
        self.state.satellite_projection_next_refresh_utc = datetime.now(
            timezone.utc
        ) + timedelta(seconds=SATELLITE_POSITION_REFRESH_INTERVAL_SECONDS)
        self.request_client_update()

    def refresh_projected_persistent_search_target(self) -> None:
        target = self.state.persistent_search_target
        if target is None:
            return
        if self._viewport_interaction_active():
            return
        if not bool(target.persistent_keep_marker):
            return
        if target.kind not in {"jpl_small_body", "jpl_body"}:
            return
        projected_altaz = project_jpl_target_altaz_from_state_vector(
            target,
            observer_lat=float(self.viewer_data.location[0]),
            observer_lon=float(self.viewer_data.location[1]),
            observer_height_m=float(self.viewer_data.observer_height_m),
            time_obj=self._current_time_obj(),
        )
        if projected_altaz is None:
            return
        alt_deg, az_deg = projected_altaz
        if (
            target.alt_deg is not None
            and target.az_deg is not None
            and abs(float(target.alt_deg) - float(alt_deg)) < 1e-9
            and abs((float(target.az_deg) % 360.0) - float(az_deg)) < 1e-9
        ):
            return
        updated_target = replace(
            target,
            alt_deg=float(alt_deg),
            az_deg=float(az_deg) % 360.0,
        )
        self.state.persistent_search_target = updated_target
        self.request_client_update()

    def _on_satellite_started(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.satellite_state.set_banner(banner)
            self.request_client_update()

    def _on_satellite_ready(self, payload: dict) -> None:
        element_epoch = payload.get("element_epoch_utc")
        if not isinstance(element_epoch, datetime):
            element_epoch = datetime.now(timezone.utc)
        banner = str(payload.get("banner", "")).strip()
        self.satellite_state.set_result(
            payload.get("records_by_group", {}),
            element_epoch_utc=element_epoch,
            refreshed_at_utc=payload.get("refreshed_at_utc"),
        )
        if banner:
            self.satellite_state.set_banner(banner)
        requested_update = False
        if float(self.satellite_opacity) > 0.0:
            self._schedule_next_satellite_refresh()
            self.reproject_satellite_overlay()
            requested_update = True
        if not requested_update:
            self.request_client_update()

    def _on_satellite_failed(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.satellite_state.set_error_banner(banner)
        if float(self.satellite_opacity) > 0.0:
            self._schedule_satellite_retry_after_failure()
        self.request_client_update()

    def start_background_terrain_horizon_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down:
            return False
        if self._viewport_interaction_active():
            return False
        if (
            self.terrain_horizon_opacity <= 0.0
            or self._terrain_horizon_controller is None
        ):
            return False
        lat, lon = self.viewer_data.location
        return self._terrain_horizon_controller.update(
            lat=lat,
            lon=lon,
            observer_height_m=self.viewer_data.observer_height_m,
            reason=reason,
        )

    def start_background_water_overlay_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down:
            return False
        if self.water_overlay_opacity <= 0.0:
            return False
        if self.terrain_horizon_opacity <= 0.0:
            return False
        if self.terrain_horizon_state.profile_altaz is None:
            return False
        if self._water_overlay_controller is None:
            return False
        surface_size_px = (
            max(1, int(self.client_width())),
            max(1, int(self.client_height())),
        )
        return self._water_overlay_controller.update(
            viewer_data=self.viewer_data,
            observer_ground_m=self._water_overlay_ground_elevation_m(),
            use_dem_ground=self._water_overlay_use_dem_ground(),
            reason=reason,
            fast_mode=self._viewport_interaction_active(),
            surface_size_px=surface_size_px,
            terrain_horizon_profile_altaz=self.state.terrain_horizon_profile,
            terrain_horizon_profile_distances_m=self.state.terrain_horizon_profile_distances_m,
        )

    def start_background_road_night_lights_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down or self.road_night_lights_opacity <= 0.0:
            return False
        if self._road_night_lights_controller is None:
            return False
        return self._road_night_lights_controller.update(
            viewer_data=self.viewer_data,
            reason=reason,
        )

    def start_background_precipitation_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down or not self._precipitation_layer_enabled():
            return False
        controller = self._precipitation_controller
        if controller is None:
            return False
        return controller.update(viewer_data=self.viewer_data, reason=reason)

    def _water_overlay_ground_elevation_m(self) -> float:
        ground_m = self.terrain_horizon_state.ground_elevation_m
        if ground_m is not None:
            return max(0.0, float(ground_m))
        viewer_ground_m = float(self.viewer_data.ground_elevation_m or 0.0)
        return max(0.0, viewer_ground_m)

    def _water_overlay_use_dem_ground(self) -> bool:
        return self.terrain_horizon_state.ground_elevation_m is not None

    def _refresh_water_overlay_active_dots(self) -> None:
        use_dem = self._water_overlay_use_dem_ground()
        self.water_overlay_state.select_active_dots(use_dem=use_dem)
        self.state.water_overlay_dots = self.water_overlay_state.dots

    def start_background_urban_outline_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down:
            return False
        if self._viewport_interaction_active():
            return False
        if self.urban_outline_opacity <= 0.0 or self._urban_outline_controller is None:
            return False
        return self._urban_outline_controller.update(
            viewer_data=self.viewer_data,
            reason=reason,
        )

    def start_background_aircraft_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down:
            return False
        if self._viewport_interaction_active():
            return False
        if float(self.aircraft_opacity) <= 0.0:
            return False
        if self._aircraft_controller is None:
            return False
        lat, lon = self.viewer_data.location
        return self._aircraft_controller.update(
            observer_lat=lat,
            observer_lon=lon,
            reason=reason,
        )

    def reproject_aircraft_overlay(self) -> None:
        if float(self.aircraft_opacity) <= 0.0:
            return
        if self._viewport_interaction_active():
            return
        snapshots = self.aircraft_state.snapshots
        if not snapshots:
            self.state.aircraft_projection_next_refresh_utc = None
            self.request_client_update()
            return
        self.state.aircraft_projection_next_refresh_utc = datetime.now(
            timezone.utc
        ) + timedelta(seconds=AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS)
        self.request_client_update()

    def _on_aircraft_started(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.aircraft_state.set_banner(banner)
            self.request_client_update()

    def _on_aircraft_ready(self, payload: dict) -> None:
        source = str(payload.get("source", "")).strip()
        refreshed_at = payload.get("refreshed_at_utc")
        if source == "rate-limited-skip":
            refreshed_at = None
        elif not isinstance(refreshed_at, datetime):
            refreshed_at = datetime.now(timezone.utc)
        banner = str(payload.get("banner", "")).strip()
        self.aircraft_state.set_result(
            payload.get("snapshots", []),
            bbox=payload.get("bbox"),
            refreshed_at_utc=refreshed_at,
        )
        if banner:
            self.aircraft_state.set_banner(banner)
        requested_update = False
        if float(self.aircraft_opacity) > 0.0:
            self._schedule_next_aircraft_refresh()
            self.reproject_aircraft_overlay()
            requested_update = True
        if not requested_update:
            self.request_client_update()

    def _on_aircraft_failed(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.aircraft_state.set_error_banner(banner)
        if float(self.aircraft_opacity) > 0.0:
            self._schedule_next_aircraft_refresh()
        self.request_client_update()

    def start_background_tropical_cyclone_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down:
            return False
        if not self._tropical_cyclone_layer_enabled():
            return False
        controller = self._tropical_cyclone_controller
        if controller is None:
            return False
        return controller.update(reason=reason)

    def _on_tropical_cyclone_started(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.tropical_cyclone_state.banner_text = banner
        self.request_client_update()

    def _on_tropical_cyclone_ready(self, payload: dict) -> None:
        snapshot_collection_payload = payload.get("snapshot_collection")
        if isinstance(snapshot_collection_payload, dict):
            from ..tropical_cyclones.models import TropicalCycloneSnapshotCollection

            parsed_collection = TropicalCycloneSnapshotCollection.from_dict(
                snapshot_collection_payload
            )
        else:
            legacy_snapshot_payload = payload.get("snapshot")
            if isinstance(legacy_snapshot_payload, dict):
                from ..tropical_cyclones.models import TropicalCycloneSnapshotCollection

                parsed_collection = TropicalCycloneSnapshotCollection.from_dict(
                    {"snapshot": legacy_snapshot_payload}
                )
            else:
                parsed_collection = None
        if parsed_collection is None:
            banner = str(payload.get("banner", "")).strip()
            if banner:
                self.tropical_cyclone_state.set_error_banner(banner)
            self.request_client_update()
            return
        cached_at = payload.get("cached_at_utc")
        if not isinstance(cached_at, datetime):
            cached_at = datetime.now(timezone.utc)
        last_checked = payload.get("last_checked_utc")
        if not isinstance(last_checked, datetime):
            last_checked = cached_at
        next_check = payload.get("next_check_utc")
        if not isinstance(next_check, datetime):
            next_check = cached_at + timedelta(minutes=90)
        next_refresh = payload.get("next_refresh_utc")
        if not isinstance(next_refresh, datetime):
            next_refresh = cached_at + timedelta(hours=3)
        banner = str(payload.get("banner", "")).strip()
        self.tropical_cyclone_state.set_result(
            parsed_collection.snapshots,
            cached_at_utc=cached_at,
            last_checked_utc=last_checked,
            next_check_utc=next_check,
            next_refresh_utc=next_refresh,
            banner_text=banner or None,
            source_url=parsed_collection.source_url or None,
            current_source=parsed_collection.service_name or None,
        )
        self.reproject_tropical_cyclone_overlay()
        self.request_client_update()

    def _on_tropical_cyclone_failed(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.tropical_cyclone_state.set_error_banner(banner)
        retry_at = datetime.now(timezone.utc) + timedelta(
            seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS
        )
        self.tropical_cyclone_state.next_check_utc = retry_at
        self.tropical_cyclone_state.next_refresh_utc = datetime.now(
            timezone.utc
        ) + timedelta(seconds=TROPICAL_CYCLONE_CACHE_TTL_SECONDS)
        self.request_client_update()

    def reproject_tropical_cyclone_overlay(
        self,
        *,
        allow_during_viewport_interaction: bool = False,
    ) -> None:
        if not self._tropical_cyclone_layer_enabled():
            return
        if (
            self._viewport_interaction_active()
            and not allow_during_viewport_interaction
        ):
            return
        state = self.tropical_cyclone_state
        snapshots = state.snapshots
        if not snapshots:
            state.projection_next_refresh_utc = None
            self.state.tropical_cyclone_projection_next_refresh_utc = None
            self.request_client_update()
            return
        next_refresh = datetime.now(timezone.utc) + timedelta(
            seconds=AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS
        )
        if state is not None:
            state.projection_next_refresh_utc = next_refresh
        self.state.tropical_cyclone_projection_next_refresh_utc = next_refresh
        self.request_client_update()

    def _on_terrain_horizon_started(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.terrain_horizon_state.banner_text = banner
        self.request_client_update()

    def _on_terrain_horizon_ready(self, payload: dict) -> None:
        self.terrain_horizon_state.set_result(
            payload["profile_altaz"],
            profile_distances_m=payload.get("profile_distances_m"),
            secondary_ridges_altaz_layers=payload.get("secondary_ridges_altaz_layers"),
            secondary_ridges_distances_m_layers=payload.get(
                "secondary_ridges_distances_m_layers"
            ),
            sample_distances_m=payload.get("sample_distances_m"),
            sample_terrain_elevation_m=payload.get("sample_terrain_elevation_m"),
            source=str(payload.get("source", "")).strip(),
        )
        ground_elevation_m = payload.get("ground_elevation_m")
        if isinstance(ground_elevation_m, (int, float)):
            ground_value = float(ground_elevation_m)
            self.terrain_horizon_state.ground_elevation_m = ground_value
            self.viewer_data = replace(
                self.viewer_data, ground_elevation_m=ground_value
            )
        self.state.terrain_horizon_profile = payload["profile_altaz"]
        self.state.terrain_horizon_profile_distances_m = payload.get(
            "profile_distances_m"
        )
        self.state.terrain_secondary_ridges_altaz_layers = payload.get(
            "secondary_ridges_altaz_layers"
        )
        self.state.terrain_secondary_ridges_distances_m_layers = payload.get(
            "secondary_ridges_distances_m_layers"
        )
        self.state.night_light_glow_profile = None
        if (
            _initial_data_load_active(self)
            and float(getattr(self, "night_light_opacity", 0.0)) <= 0.0
            and float(getattr(self, "ridge_glow_opacity", 0.0)) <= 0.0
        ):
            self._startup_initial_night_light_loaded = True
        self._refresh_water_overlay_active_dots()
        startup_initial_load = _initial_data_load_active(self)
        if startup_initial_load:
            self._startup_initial_terrain_loaded = True
        if not self._is_shutting_down:
            if hasattr(self, "start_background_road_night_lights_update"):
                self.start_background_road_night_lights_update(reason="terrain-ready")
            if startup_initial_load:
                if self.water_overlay_opacity > 0.0:
                    self.start_background_water_overlay_update(reason="initial")
            else:
                self.start_background_water_overlay_update(reason="terrain-ready")
        self._sync_water_overlay_action_enabled()
        self._compositor.invalidate()
        self.request_client_update()
        if not self._is_shutting_down:
            self.request_sky_data_update(reason="terrain-ready")
        if startup_initial_load and self.water_overlay_opacity <= 0.0:
            self._continue_initial_data_load()
            return

    def _on_terrain_horizon_failed(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        self.terrain_horizon_state.clear_profile()
        self.state.terrain_horizon_profile = None
        self.state.terrain_horizon_profile_distances_m = None
        self.state.terrain_secondary_ridges_altaz_layers = None
        self.state.terrain_secondary_ridges_distances_m_layers = None
        self.state.night_light_glow_profile = None
        if _initial_data_load_active(self):
            self._startup_initial_night_light_loaded = True
        self.terrain_horizon_state.sample_distances_m = None
        self.terrain_horizon_state.sample_terrain_elevation_m = None
        self._refresh_water_overlay_active_dots()
        if banner:
            self.terrain_horizon_state.set_error_banner(banner)
        if _initial_data_load_active(self):
            self._startup_initial_terrain_loaded = True
        self._sync_water_overlay_action_enabled()
        self._compositor.invalidate()
        self.request_client_update()
        if _initial_data_load_active(self):
            self._continue_initial_data_load()

    def _on_water_overlay_started(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner and self.water_overlay_state.dots is None:
            self.water_overlay_state.banner_text = banner
        self.request_client_update()

    def _on_road_night_lights_started(self) -> None:
        self.road_night_lights_status = "loading"
        self.request_client_update()

    def _on_road_night_lights_ready(self, payload: dict) -> None:
        polylines = payload.get("polylines")
        self.road_night_light_polylines = (
            polylines if isinstance(polylines, list) else None
        )
        self.road_night_lights_status = (
            str(len(polylines))
            if isinstance(polylines, list)
            else "unavailable"
        )
        self._compositor.invalidate()
        self.request_client_update()

    def _on_road_night_lights_failed(self, payload: dict) -> None:
        self.road_night_light_polylines = None
        self.road_night_lights_status = "unavailable"
        logger.warning(
            "%s", str(payload.get("banner", "Road night lights unavailable"))
        )
        self._compositor.invalidate()
        self.request_client_update()

    def _on_precipitation_started(self) -> None:
        self.precipitation_status = "loading"
        self.request_client_update()

    def _on_precipitation_ready(self, payload: dict) -> None:
        self.state.precipitation_columns = list(payload.get("columns") or [])
        self.precipitation_forecast_time_utc = payload.get("forecast_time_utc")
        self.precipitation_interval_seconds = payload.get("interval_seconds")
        self.precipitation_status = "ready"
        self.state.precipitation_next_refresh_utc = datetime.now(
            timezone.utc
        ) + timedelta(seconds=PRECIPITATION_REFRESH_SECONDS)
        self._compositor.invalidate()
        self.request_client_update()

    def _on_precipitation_failed(self) -> None:
        self.state.precipitation_columns = None
        self.precipitation_status = "unavailable"
        self.state.precipitation_next_refresh_utc = datetime.now(
            timezone.utc
        ) + timedelta(seconds=PRECIPITATION_REFRESH_SECONDS)
        self._compositor.invalidate()
        self.request_client_update()

    def _on_water_overlay_ready(self, payload: dict) -> None:
        dots = payload.get("dots")
        sea_dots = payload.get("sea_dots")
        inland_dots = payload.get("inland_dots")
        dem_dots = payload.get("dem_dots")
        water_polylines = payload.get("water_polylines")
        mode = str(payload.get("mode", "")).strip().lower() or "sea"
        source = str(payload.get("source", "")).strip() or "ready"
        count = len(dots) if isinstance(dots, list) else 0
        logger.info(
            "Water surface dots ready: mode=%s source=%s dots=%d sea_dots=%s inland_dots=%s dem_dots=%s",
            mode,
            source,
            count,
            len(sea_dots) if isinstance(sea_dots, list) else "-",
            len(inland_dots) if isinstance(inland_dots, list) else "-",
            len(dem_dots) if isinstance(dem_dots, list) else "-",
        )
        if mode == "dem":
            self.water_overlay_state.set_dem_dots_result(
                dem_dots or dots, source=source
            )
        else:
            self.water_overlay_state.set_sea_level_dots_result(
                sea_dots or dots, source=source
            )
        if isinstance(sea_dots, list):
            self.water_overlay_state.sea_level_dots = sea_dots
        if isinstance(inland_dots, list):
            self.water_overlay_state.inland_dots = inland_dots
        if isinstance(dem_dots, list):
            self.water_overlay_state.dem_dots = dem_dots
        if isinstance(water_polylines, list):
            self.water_overlay_state.set_polylines(water_polylines)
        self._refresh_water_overlay_active_dots()
        if _initial_data_load_active(self):
            self._startup_initial_water_loaded = True
        self._compositor.invalidate()
        self.request_client_update()
        if _initial_data_load_active(self):
            self._continue_initial_data_load()

    def _on_water_overlay_failed(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.water_overlay_state.set_error_banner(banner)
        self._refresh_water_overlay_active_dots()
        if _initial_data_load_active(self):
            self._startup_initial_water_loaded = True
        self._compositor.invalidate()
        self.request_client_update()
        if _initial_data_load_active(self):
            self._continue_initial_data_load()

    def _on_urban_outline_started(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.urban_outline_state.banner_text = banner
        self.request_client_update()

    def _on_urban_outline_ready(self, payload: dict) -> None:
        outlines = payload.get("outlines")
        self.urban_outline_state.set_result(
            outlines,
            source=str(payload.get("source", "")).strip() or "ready",
            base_outline_count=payload.get("base_outline_count"),
            skyscraper_outline_count=payload.get("skyscraper_outline_count"),
        )
        self.state.urban_outlines = outlines
        if _initial_data_load_active(self):
            self._startup_initial_urban_loaded = True
        self._compositor.invalidate()
        self.request_client_update()
        if _initial_data_load_active(self):
            self._continue_initial_data_load()

    def _on_urban_outline_failed(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        self.urban_outline_state.clear_outlines()
        self.state.urban_outlines = None
        if banner:
            self.urban_outline_state.set_error_banner(banner)
        if _initial_data_load_active(self):
            self._startup_initial_urban_loaded = True
        self._compositor.invalidate()
        self.request_client_update()
        if _initial_data_load_active(self):
            self._continue_initial_data_load()
