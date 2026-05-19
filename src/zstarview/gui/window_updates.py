from __future__ import annotations

import logging
import os
import time
from dataclasses import replace
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, Optional

from ..aircraft_constants import AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS
from ..astro import load_ephemeris
from ..clouddisc.providers.select import GOES_SATELLITES
from ..paths import CACHE_PATH, CLOUD_UPDATE_INTERVAL
from ..satellite_constants import SATELLITE_POSITION_REFRESH_INTERVAL_SECONDS
from ..search.jpl import project_jpl_target_altaz_from_state_vector
from ..render import geometry as render_geometry

logger = logging.getLogger(__name__)
_STATUS_CLOUD = "☁"
_STATUS_WATER = "W"
_STATUS_SATELLITE = "🛰"
_STATUS_AIRCRAFT = "✈"
_STATUS_TERRAIN = "△"
_STATUS_URBAN = "🂓"


def _status_segment(icon: str, text: str, *, hidden: bool = False) -> str:
    if hidden:
        return f"{icon} ---"
    clean = str(text).strip()
    return f"{icon} {clean}" if clean else icon


def _strip_status_prefix(text: str, prefix: str) -> str:
    clean = str(text).strip()
    if clean.casefold().startswith(prefix.casefold()):
        clean = clean[len(prefix):].strip()
    return clean


def _cloud_satellite_group(satellite: str) -> str:
    sat = str(satellite).strip()
    if sat in GOES_SATELLITES:
        return "GOES"
    if sat == "HIMAWARI":
        return "HIMAWARI"
    return sat

def _initial_data_load_active(obj: object) -> bool:
    return bool(obj._startup_initial_load_started) and not bool(
        obj._startup_initial_data_loaded
    )


class SkyWindowUpdatesMixin:
    def _viewport_interaction_active(self) -> bool:
        return bool(self.state.viewport_interaction_mode)

    def _background_updates_busy(self) -> bool:
        controllers = (
            getattr(self, "_sky_worker", None),
            getattr(self, "_cloud_controller", None),
            getattr(self, "_satellite_controller", None),
            getattr(self, "_aircraft_controller", None),
            getattr(self, "_jpl_small_body_controller", None),
            getattr(self, "_terrain_horizon_controller", None),
            getattr(self, "_water_overlay_controller", None),
            getattr(self, "_urban_outline_controller", None),
        )
        for controller in controllers:
            has_in_flight_update = getattr(controller, "has_in_flight_update", None)
            if callable(has_in_flight_update) and has_in_flight_update():
                return True
        return False

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

    def _water_overlay_action_enabled(self) -> bool:
        if not bool(self._water_overlay_gui_allowed):
            return False
        return float(self.terrain_horizon_opacity) > 0.0

    def _sync_water_overlay_action_enabled(self) -> None:
        action = self._action_toggle_water_overlay
        if action is not None:
            action.setEnabled(self._water_overlay_action_enabled())

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
        background_updates_busy = self._background_updates_busy()

        sky_next_refresh = self.state.sky_next_refresh_utc
        if (
            not background_updates_busy
            and self.state.sky_update_pending
        ):
            started = self.start_background_sky_data_update(
                star_vmag_limit=self.state.pending_star_vmag_limit,
                reason="scheduler",
                allow_during_viewport_interaction=False,
            )
            if started:
                self.state.sky_update_pending = False
                self.state.pending_star_vmag_limit = None
                self.state.sky_next_refresh_utc = now_utc + timedelta(seconds=self.sky_update_interval)
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
                self.state.sky_next_refresh_utc = now_utc + timedelta(seconds=self.sky_update_interval)
                return

        persistent_next_refresh = self.state.persistent_search_next_refresh_utc
        if (
            not background_updates_busy
            and isinstance(persistent_next_refresh, datetime)
            and now_utc >= persistent_next_refresh
        ):
            started = self._start_persistent_search_refresh(reason="scheduler")
            if started:
                self.state.persistent_search_next_refresh_utc = persistent_next_refresh + timedelta(hours=1)
                return

        cloud_next_refresh = self.state.cloud_next_refresh_utc
        if (
            not background_updates_busy
            and isinstance(cloud_next_refresh, datetime)
            and now_utc >= cloud_next_refresh
            and self._clouddisc
            and self.cloud_disc_alpha > 0.0
        ):
            if self._cloud_controller is None or self._cloud_controller.has_in_flight_update():
                return
            self.start_background_cloud_update(reason="scheduler")
            if self._cloud_controller.has_in_flight_update():
                self.state.cloud_next_refresh_utc = now_utc + timedelta(seconds=CLOUD_UPDATE_INTERVAL)
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

        # Lowest-priority idle work: keep satellite and aircraft positions fresh.
        if self._aircraft_layer_enabled() and self._aircraft_projection_next_refresh_delay_ms() == 0:
            self.reproject_aircraft_overlay()
            return

        if self._satellite_layer_enabled() and self._satellite_projection_next_refresh_delay_ms() == 0:
            self.reproject_satellite_overlay()
            return

    def _resolve_aircraft_debug_snapshot_dir(self) -> Path | None:
        raw = os.getenv("ZSTARVIEW_DEBUG_SAVE_AIRCRAFT_READY_FRAME", "").strip()
        if not raw:
            return None
        lowered = raw.lower()
        if lowered in {"0", "false", "no", "off"}:
            return None
        if lowered in {"1", "true", "yes", "on"}:
            return Path(CACHE_PATH) / "debug" / "aircraft-ready"
        return Path(raw).expanduser()

    def _queue_aircraft_debug_snapshot(self, payload: Dict) -> None:
        output_dir = self._resolve_aircraft_debug_snapshot_dir()
        if output_dir is None:
            return
        source = str(payload.get("source", "")).strip().lower() or "ready"
        if source == "cache-fresh":
            return
        self._aircraft_debug_snapshot_payload = dict(payload)
        self._aircraft_debug_snapshot_save_queued = False

    def _save_aircraft_debug_snapshot_image(
        self,
        image,
        payload: Dict,
    ) -> None:
        current_payload = getattr(self, "_aircraft_debug_snapshot_payload", None)
        if current_payload is not payload:
            return
        try:
            output_dir = self._resolve_aircraft_debug_snapshot_dir()
            if output_dir is None:
                return
            refreshed_at = payload.get("refreshed_at_utc")
            if not isinstance(refreshed_at, datetime):
                refreshed_at = datetime.now(timezone.utc)
            source = str(payload.get("source", "")).strip().lower() or "ready"
            if source == "cache-fresh":
                return
            safe_source = "".join(
                ch if (ch.isascii() and (ch.isalnum() or ch in {"-", "_", "."})) else "-"
                for ch in source
            ).strip("-")
            if not safe_source:
                safe_source = "ready"
            filename = f"aircraft-ready-{refreshed_at.strftime('%Y%m%dT%H%M%SZ')}-{safe_source}.png"
            output_dir.mkdir(parents=True, exist_ok=True)
            output_path = output_dir / filename
            if not image.save(str(output_path), "PNG"):
                logger.warning(
                    "Failed to save aircraft debug snapshot: %s", output_path
                )
                return
            logger.info("Saved aircraft debug snapshot: %s", output_path)
        except Exception as exc:
            logger.warning("Aircraft debug snapshot failed: %s", exc, exc_info=True)
        finally:
            current_payload = getattr(self, "_aircraft_debug_snapshot_payload", None)
            if current_payload is payload:
                self._aircraft_debug_snapshot_payload = None
                self._aircraft_debug_snapshot_save_queued = False

    def _status_line_message(self) -> str:
        vertical_bar = "\u23ae"
        parts: list[str] = []
        cloud_message = self._cloud_status_line()
        if cloud_message:
            parts.append(cloud_message)
        satellite_status_line = self._satellite_status_line
        if callable(satellite_status_line):
            satellite_message = satellite_status_line()
            if satellite_message:
                parts.append(satellite_message)
        aircraft_status_line = self._aircraft_status_line
        if callable(aircraft_status_line):
            aircraft_message = aircraft_status_line()
            if aircraft_message:
                parts.append(aircraft_message)
        jpl_message = self._jpl_small_body_status_line()
        if jpl_message:
            parts.append(jpl_message)
        terrain_message = self._terrain_horizon_status_line()
        if terrain_message:
            parts.append(terrain_message)
        water_message = self._water_overlay_status_line()
        if water_message:
            parts.append(water_message)
        urban_message = self._urban_outline_status_line()
        if urban_message:
            parts.append(urban_message)
        return f"{vertical_bar} %s {vertical_bar}" % f" {vertical_bar} ".join(parts)

    def _safe_request_cloud_repaint(self) -> None:
        """Best-effort repaint request; ignores teardown-time signal errors."""
        if self._is_shutting_down:
            return
        try:
            self.cloud_repaint_requested.emit()
        except RuntimeError:
            logger.debug("Skip cloud repaint emit during shutdown.")

    def _cloud_status_line(self) -> str:
        cloud_disc_alpha = float(self.cloud_disc_alpha)
        if cloud_disc_alpha <= 0.0:
            return _status_segment(_STATUS_CLOUD, "", hidden=True)
        sat = self.cloud_state.current_satellite or self._predicted_cloud_satellite()
        sat_group = _cloud_satellite_group(sat)
        if self.cloud_state.banner_text:
            detail = _strip_status_prefix(self.cloud_state.banner_text, "Clouds:")
            detail_lower = detail.lower()
            if detail_lower.startswith("downloading"):
                return _status_segment(_STATUS_CLOUD, f"{sat_group} downloading")
            if any(token in detail_lower for token in ("timed out", "error", "failed", "failure")):
                return _status_segment(_STATUS_CLOUD, f"{sat_group} failed")
            return _status_segment(_STATUS_CLOUD, f"{sat_group} {detail}")
        meta = self.cloud_state.meta
        if meta is not None:
            try:
                t = meta.time_utc.strftime("%H:%MZ")
                sat_group = _cloud_satellite_group(getattr(meta, "satellite", sat))
                coverage = self.cloud_state.coverage_ratio
                if coverage is not None and coverage < 0.999:
                    pct = int(round(max(0.0, min(1.0, float(coverage))) * 100.0))
                    return _status_segment(_STATUS_CLOUD, f"{sat_group} {pct}% {t}")
                return _status_segment(_STATUS_CLOUD, f"{sat_group} {t}")
            except Exception:
                pass
        return _status_segment(_STATUS_CLOUD, f"{sat_group} idle")

    def _terrain_horizon_status_line(self) -> str:
        if self.terrain_horizon_opacity <= 0.0:
            return _status_segment(_STATUS_TERRAIN, "", hidden=True)
        if self.terrain_horizon_state.banner_text:
            detail = _strip_status_prefix(
                self.terrain_horizon_state.banner_text,
                "Terrain horizon:",
            )
            return _status_segment(_STATUS_TERRAIN, detail)
        if self.terrain_horizon_state.current_source:
            detail = _strip_status_prefix(
                self.terrain_horizon_state.current_source,
                "Dem:",
            )
            return _status_segment(_STATUS_TERRAIN, detail)
        return ""

    def _water_overlay_status_line(self) -> str:
        if self.water_overlay_opacity <= 0.0:
            return _status_segment(_STATUS_WATER, "", hidden=True)
        state = self.water_overlay_state
        if state.banner_text:
            detail = _strip_status_prefix(state.banner_text, "Water:")
            return _status_segment(_STATUS_WATER, detail)
        sea_count = "?" if state.sea_level_dots is None else str(len(state.sea_level_dots))
        inland_count = "?" if state.inland_dots is None else str(len(state.inland_dots))
        return _status_segment(_STATUS_WATER, f"{sea_count}+{inland_count}")

    def _urban_outline_status_line(self) -> str:
        if self.urban_outline_opacity <= 0.0:
            return _status_segment(_STATUS_URBAN, "", hidden=True)
        if self.urban_outline_state.banner_text:
            detail = _strip_status_prefix(
                self.urban_outline_state.banner_text,
                "Urban outline:",
            )
            return _status_segment(_STATUS_URBAN, detail)
        if self.urban_outline_state.current_source:
            detail = _strip_status_prefix(
                self.urban_outline_state.current_source,
                "Urban:",
            )
            return _status_segment(_STATUS_URBAN, detail)
        return ""

    def _aircraft_status_line(self) -> str:
        if float(self.aircraft_opacity) <= 0.0:
            return _status_segment(_STATUS_AIRCRAFT, "", hidden=True)
        aircraft_state = self.aircraft_state
        if aircraft_state is None:
            return ""
        if aircraft_state.banner_text:
            detail = _strip_status_prefix(aircraft_state.banner_text, "Aircraft:")
            return _status_segment(_STATUS_AIRCRAFT, detail)
        if aircraft_state.last_success_utc is None:
            return _status_segment(_STATUS_AIRCRAFT, "idle")
        return _status_segment(_STATUS_AIRCRAFT, aircraft_state.last_success_utc.strftime("%H:%MZ"))

    def _satellite_status_line(self) -> str:
        if float(self.satellite_opacity) <= 0.0:
            return _status_segment(_STATUS_SATELLITE, "", hidden=True)
        satellite_state = self.satellite_state
        if satellite_state is None:
            return ""
        if satellite_state.banner_text:
            detail = _strip_status_prefix(satellite_state.banner_text, "Satellites:")
            return _status_segment(_STATUS_SATELLITE, detail)
        if satellite_state.element_epoch_utc is None:
            return _status_segment(_STATUS_SATELLITE, "idle")
        return _status_segment(_STATUS_SATELLITE, satellite_state.element_epoch_utc.strftime("%H:%MZ"))

    def _jpl_small_body_status_line(self) -> str:
        target = self.state.persistent_search_target
        if target is None:
            return ""
        label = str(target.label).strip()
        if not label:
            return ""
        altaz_suffix = self._target_altaz_suffix(target)
        next_refresh_utc = self.state.persistent_search_next_refresh_utc
        last_error = str(self.state.persistent_search_last_error).strip()
        if last_error.casefold() == "none":
            last_error = ""
        if target.jpl_group != "sb" and not last_error:
            return f"JPL [{label}]: held{altaz_suffix}"
        if isinstance(next_refresh_utc, datetime):
            refresh_part = next_refresh_utc.strftime("%H:%MZ")
        else:
            refresh_part = "pending"
        if last_error:
            return f"JPL [{label}]: retry {refresh_part} ({last_error}){altaz_suffix}"
        if target.jpl_group == "sb":
            return f"JPL [{label}]: retry {refresh_part}{altaz_suffix}"
        return f"JPL [{label}]: refresh {refresh_part}{altaz_suffix}"

    def _target_altaz_suffix(self, target: object) -> str:
        alt = target.alt_deg
        az = target.az_deg
        if alt is None or az is None:
            return ""
        try:
            alt_deg = float(alt)
            az_deg = float(az) % 360.0
        except (TypeError, ValueError):
            return ""
        return f" [alt={alt_deg:.1f} az={az_deg:.1f}]"

    def _on_sky_data_calculated(self, payload: Dict) -> None:
        current_generation = int(self._disc_generation)
        payload_generation = int(payload.get("render_generation", current_generation))
        payload_geometry = payload.get("geometry")
        current_geometry = render_geometry.get_screen_geometry(
            max(2, int(self.client_width())),
            max(2, int(self.client_height())),
            self.viewer_data.view_alt_deg,
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
                self.request_sky_data_update(reason="stale-render")
            return
        if not self.state.viewport_interaction_mode:
            view_center = payload.get("view_center", self.viewer_data.view_center)
            current_view_center = tuple(self.viewer_data.view_center)
            if not (
                isinstance(view_center, (tuple, list))
                and len(view_center) >= 2
                and abs(float(view_center[0]) - float(current_view_center[0])) < 1e-9
                and abs(float(view_center[1]) - float(current_view_center[1])) < 1e-9
            ):
                logger.debug(
                    "Discard stale sky payload view_center=%s current=%s generation=%s",
                    view_center,
                    current_view_center,
                    payload_generation,
                )
                if not self._is_shutting_down:
                    self.request_sky_data_update(reason="stale-view-center")
                return
            if isinstance(view_center, (tuple, list)) and len(view_center) >= 2:
                self.state.render_view_center = (
                    float(view_center[0]),
                    float(view_center[1]),
                )
            else:
                self.state.render_view_center = tuple(self.viewer_data.view_center)
        self.state.celestial_data = payload["celestial"]
        self.state.sky_disc_image = payload["sky_disc"]
        self.state.night_light_glow_profile = payload.get("night_light_glow_profile")

        if self.state.viewport_interaction_release_pending:
            self.state.viewport_interaction_release_pending = False
            self.state.viewport_interaction_mode = False
            self.state.viewport_interaction_stars = None
            if not self._is_shutting_down:
                self.start_background_cloud_update(reason="view-change-release")
                self.start_background_terrain_horizon_update(
                    reason="view-change-release"
                )

        self._compositor.invalidate()
        self.request_client_update()

        if _initial_data_load_active(self):
            self._startup_initial_sky_loaded = True
            self._continue_initial_data_load()
            return

        self.state.sky_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=self.sky_update_interval
        )
        if self._clouddisc and self.cloud_disc_alpha > 0.0:
            self.state.cloud_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
                seconds=CLOUD_UPDATE_INTERVAL
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
        if self.terrain_horizon_opacity > 0.0 and not self._startup_initial_terrain_loaded:
            if self.start_background_terrain_horizon_update(reason="initial"):
                return
            return
        if (
            self.terrain_horizon_opacity > 0.0
            and self.water_overlay_opacity > 0.0
            and self.terrain_horizon_state.profile_altaz is not None
            and not self._startup_initial_water_loaded
        ):
            if self.start_background_water_overlay_update(reason="initial"):
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

    def request_sky_data_update(
        self,
        star_vmag_limit: Optional[float] = None,
        *,
        reason: str = "manual",
        allow_during_viewport_interaction: bool = False,
    ) -> None:
        if self._viewport_interaction_active() and not allow_during_viewport_interaction:
            self.state.sky_update_pending = True
            self.state.pending_star_vmag_limit = star_vmag_limit
            logger.debug(
                "Sky data update deferred during viewport interaction (reason=%s).",
                reason,
            )
            return
        if self.start_background_sky_data_update(
            star_vmag_limit=star_vmag_limit,
            reason=reason,
            allow_during_viewport_interaction=allow_during_viewport_interaction,
        ):
            self.state.sky_update_pending = False
            self.state.pending_star_vmag_limit = None
            return
        self.state.sky_update_pending = True
        self.state.pending_star_vmag_limit = star_vmag_limit
        logger.debug(
            "Sky data update deferred; worker is busy (reason=%s).",
            reason,
        )

    def start_background_sky_data_update(
        self,
        is_initial_load: bool = False,
        star_vmag_limit: Optional[float] = None,
        reason: str = "manual",
        allow_during_viewport_interaction: bool = False,
    ) -> bool:
        if self._viewport_interaction_active() and not allow_during_viewport_interaction:
            return False
        lat, lon = self.viewer_data.location
        use_lod6_catalog = star_vmag_limit is not None and float(star_vmag_limit) <= 6.0
        star_catalog = self.star_catalog_np
        star_subset_indices = (
            self.star_catalog_lod6_indices if use_lod6_catalog else None
        )
        worker_star_vmag_limit = None if use_lod6_catalog else star_vmag_limit
        ephemeris = self._ephemeris
        if ephemeris is None:
            ephemeris = load_ephemeris()
        started = self._sky_worker.update(
            ephemeris=ephemeris,
            viewer_data=self.viewer_data,
            geometry=render_geometry.get_screen_geometry(
                max(2, int(self.client_width())),
                max(2, int(self.client_height())),
                self.viewer_data.view_alt_deg,
            ),
            star_catalog=star_catalog,
            dso_catalog=self.dso_catalog_np,
            star_vmag_limit=worker_star_vmag_limit,
            star_subset_indices=star_subset_indices,
            delta_t=self.delta_t,
            sky_disc_alpha=self.sky_disc_alpha,
            sky_disc_style=self.sky_disc_style,
            theme=self.theme,
            star_catalog_meta=self.star_catalog_meta,
            image_size=(
                max(2, int(self.client_width())),
                max(2, int(self.client_height())),
            ),
            render_generation=int(self._disc_generation),
        )
        if started:
            if is_initial_load:
                logger.info("Calculating initial sky data (reason=%s)...", reason)
            else:
                logger.info("Updating sky data (reason=%s)...", reason)
        return started

    def start_background_cloud_update(self, reason: str = "manual") -> None:
        if self._is_shutting_down:
            return
        if self._viewport_interaction_active():
            return
        if not (self._cloud_controller and self.cloud_disc_alpha > 0.0):
            return
        lat, lon = self.viewer_data.location
        alt, az = self.viewer_data.view_center
        self._cloud_controller.update(
            lat=lat,
            lon=lon,
            alt=alt,
            az=az,
            radius_px=self.state.cloud_base_size,
            content_fov_deg=float(self.content_fov_deg),
            reason=reason,
            render_generation=int(self._disc_generation),
        )

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
            enabled_groups=tuple(
                self._enabled_satellite_groups
            ),
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
        self.state.satellite_projection_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=SATELLITE_POSITION_REFRESH_INTERVAL_SECONDS
        )
        self.request_client_update()

    def refresh_projected_satellite_overlay(self) -> None:
        self.reproject_satellite_overlay()

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

    def _on_cloud_started(self, payload: Dict) -> None:
        sat = str(payload.get("satellite", "")).strip()
        banner = str(payload.get("banner", "")).strip()
        if sat:
            self.cloud_state.current_satellite = sat
        if banner:
            self.cloud_state.set_error_banner(banner)

    def _on_satellite_started(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.satellite_state.set_banner(banner)
            self.request_client_update()

    def _on_satellite_ready(self, payload: Dict) -> None:
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

    def _on_satellite_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.satellite_state.set_error_banner(banner)
        if float(self.satellite_opacity) > 0.0:
            self._schedule_satellite_retry_after_failure()
        self.request_client_update()

    def _on_cloud_ready(self, payload: Dict) -> None:
        current_generation = int(self._disc_generation)
        payload_generation = int(payload.get("render_generation", current_generation))
        if payload_generation != current_generation:
            logger.debug(
                "Discard stale cloud payload generation=%s current=%s",
                payload_generation,
                current_generation,
            )
            if not self._is_shutting_down:
                self.start_background_cloud_update(reason="stale-render")
            return
        self.cloud_state.set_result(
            payload["image"],
            payload.get("meta"),
            az=float(payload["az"]),
            time_utc=payload["time_utc"],
            cloud_amount_field=payload.get("cloud_amount_field"),
            missing_mask=payload.get("missing_mask"),
            coverage_ratio=payload.get("coverage_ratio"),
            source_key=payload.get("source_key"),
            request_id=payload.get("request_id"),
        )
        self.state.cloud_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=CLOUD_UPDATE_INTERVAL
        )
        self._compositor.invalidate()
        if self.state.interaction_mode:
            self.state.cloud_repaint_deferred = True
            return
        self._safe_request_cloud_repaint()

    def _on_cloud_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.cloud_state.set_error_banner(banner)
        self.state.cloud_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=CLOUD_UPDATE_INTERVAL
        )
        if self.state.interaction_mode:
            self.state.cloud_repaint_deferred = True
        if banner:
            self._safe_request_cloud_repaint()

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
        if self._water_overlay_controller is None:
            return False
        return self._water_overlay_controller.update(
            viewer_data=self.viewer_data,
            observer_ground_m=self._water_overlay_ground_elevation_m(),
            use_dem_ground=self._water_overlay_use_dem_ground(),
            reason=reason,
            terrain_horizon_profile_altaz=self.state.terrain_horizon_profile,
            terrain_horizon_profile_distances_m=self.state.terrain_horizon_profile_distances_m,
            terrain_horizon_secondary_profile_altaz_layers=self.state.terrain_horizon_secondary_profile_altaz_layers,
            terrain_horizon_secondary_profile_distances_m_layers=self.state.terrain_horizon_secondary_profile_distances_m_layers,
        )

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
            observer_height_m=self.viewer_data.observer_height_m,
            time_obj=self._current_time_obj(),
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
        self.state.aircraft_projection_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS
        )
        self.request_client_update()

    def refresh_projected_aircraft_overlay(self) -> None:
        self.reproject_aircraft_overlay()

    def _on_aircraft_started(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.aircraft_state.set_banner(banner)
            self.request_client_update()

    def _on_aircraft_ready(self, payload: Dict) -> None:
        refreshed_at = payload.get("refreshed_at_utc")
        if not isinstance(refreshed_at, datetime):
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
        self._queue_aircraft_debug_snapshot(payload)

    def _on_aircraft_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.aircraft_state.set_error_banner(banner)
        if float(self.aircraft_opacity) > 0.0:
            self._schedule_next_aircraft_refresh()
        self.request_client_update()

    def _on_terrain_horizon_started(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.terrain_horizon_state.banner_text = banner
        self.request_client_update()

    def _on_terrain_horizon_ready(self, payload: Dict) -> None:
        self.terrain_horizon_state.set_result(
            payload["profile_altaz"],
            profile_distances_m=payload.get("profile_distances_m"),
            secondary_profile_altaz_layers=payload.get("secondary_profile_altaz_layers"),
            secondary_profile_distances_m_layers=payload.get(
                "secondary_profile_distances_m_layers"
            ),
            source=str(payload.get("source", "")).strip(),
        )
        ground_elevation_m = payload.get("ground_elevation_m")
        if isinstance(ground_elevation_m, (int, float)):
            ground_value = float(ground_elevation_m)
            self.terrain_horizon_state.ground_elevation_m = ground_value
            self.viewer_data = replace(self.viewer_data, ground_elevation_m=ground_value)
        self.state.terrain_horizon_profile = payload["profile_altaz"]
        self.state.terrain_horizon_profile_distances_m = payload.get("profile_distances_m")
        self.state.terrain_horizon_secondary_profile_altaz_layers = payload.get(
            "secondary_profile_altaz_layers"
        )
        self.state.terrain_horizon_secondary_profile_distances_m_layers = payload.get(
            "secondary_profile_distances_m_layers"
        )
        self._refresh_water_overlay_active_dots()
        startup_initial_load = _initial_data_load_active(self)
        if startup_initial_load:
            self._startup_initial_terrain_loaded = True
        if not self._is_shutting_down:
            if startup_initial_load:
                if self.water_overlay_opacity > 0.0:
                    self.start_background_water_overlay_update(reason="initial")
            else:
                self.start_background_water_overlay_update(reason="terrain-ready")
        self._sync_water_overlay_action_enabled()
        self._compositor.invalidate()
        self.request_client_update()
        if startup_initial_load and self.water_overlay_opacity <= 0.0:
            self._continue_initial_data_load()
            return

    def _on_terrain_horizon_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        self.terrain_horizon_state.clear_profile()
        self.state.terrain_horizon_profile = None
        self.state.terrain_horizon_profile_distances_m = None
        self.state.terrain_horizon_secondary_profile_altaz_layers = None
        self.state.terrain_horizon_secondary_profile_distances_m_layers = None
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

    def _on_water_overlay_started(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner and self.water_overlay_state.dots is None:
            self.water_overlay_state.banner_text = banner
        self.request_client_update()

    def _on_water_overlay_ready(self, payload: Dict) -> None:
        dots = payload.get("dots")
        sea_dots = payload.get("sea_dots")
        inland_dots = payload.get("inland_dots")
        dem_dots = payload.get("dem_dots")
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
            self.water_overlay_state.set_dem_dots_result(dem_dots or dots, source=source)
        else:
            self.water_overlay_state.set_sea_level_dots_result(sea_dots or dots, source=source)
        if isinstance(sea_dots, list):
            self.water_overlay_state.sea_level_dots = sea_dots
        if isinstance(inland_dots, list):
            self.water_overlay_state.inland_dots = inland_dots
        if isinstance(dem_dots, list):
            self.water_overlay_state.dem_dots = dem_dots
        self._refresh_water_overlay_active_dots()
        if _initial_data_load_active(self):
            self._startup_initial_water_loaded = True
        self._compositor.invalidate()
        self.request_client_update()
        if _initial_data_load_active(self):
            self._continue_initial_data_load()

    def _on_water_overlay_failed(self, payload: Dict) -> None:
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

    def _on_urban_outline_started(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.urban_outline_state.banner_text = banner
        self.request_client_update()

    def _on_urban_outline_ready(self, payload: Dict) -> None:
        outlines = payload.get("outlines")
        self.urban_outline_state.set_result(
            outlines,
            source=str(payload.get("source", "")).strip() or "ready",
        )
        self.state.urban_outlines = outlines
        if _initial_data_load_active(self):
            self._startup_initial_urban_loaded = True
        self._compositor.invalidate()
        self.request_client_update()
        if _initial_data_load_active(self):
            self._continue_initial_data_load()

    def _on_urban_outline_failed(self, payload: Dict) -> None:
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
