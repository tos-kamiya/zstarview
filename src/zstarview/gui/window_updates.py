from __future__ import annotations

import logging
import os
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional

from ..aircraft import project_aircraft_snapshots
from ..paths import CACHE_PATH
from ..satellites import project_satellite_records
from ..search.jpl import project_jpl_target_altaz_from_state_vector

logger = logging.getLogger(__name__)
_JPL_DEBUG_ENV = "ZSTARVIEW_DEBUG_JPL_SEARCH"
_STATUS_CLOUD = "☁"
_STATUS_SATELLITE = "🛰"
_STATUS_AIRCRAFT = "✈"
_STATUS_TERRAIN = "▲"
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


def _jpl_debug_enabled() -> bool:
    raw = os.getenv(_JPL_DEBUG_ENV, "").strip().casefold()
    return raw in {"1", "true", "yes", "on"}


def _jpl_debug_print(message: str) -> None:
    if not _jpl_debug_enabled():
        return
    print(f"[jpl-debug] {message}", file=sys.stderr, flush=True)


class SkyWindowUpdatesMixin:
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

    def _maybe_save_aircraft_debug_snapshot(self, payload: Dict) -> None:
        output_dir = self._resolve_aircraft_debug_snapshot_dir()
        if output_dir is None:
            return
        render_current_image = self.render_current_image
        if not callable(render_current_image):
            return
        refreshed_at = payload.get("refreshed_at_utc")
        if not isinstance(refreshed_at, datetime):
            refreshed_at = datetime.now(timezone.utc)
        source = str(payload.get("source", "")).strip().lower() or "ready"
        safe_source = "".join(
            ch if (ch.isascii() and (ch.isalnum() or ch in {"-", "_", "."})) else "-"
            for ch in source
        ).strip("-")
        if not safe_source:
            safe_source = "ready"
        filename = f"aircraft-ready-{refreshed_at.strftime('%Y%m%dT%H%M%SZ')}-{safe_source}.png"
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            image = render_current_image(include_hud=True)
            output_path = output_dir / filename
            if not image.save(str(output_path), "PNG"):
                logger.warning(
                    "Failed to save aircraft debug snapshot: %s", output_path
                )
                return
            logger.info("Saved aircraft debug snapshot: %s", output_path)
        except Exception as exc:
            logger.warning("Aircraft debug snapshot failed: %s", exc, exc_info=True)

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
        if self.cloud_state.banner_text:
            detail = _strip_status_prefix(self.cloud_state.banner_text, "Clouds:")
            detail_lower = detail.lower()
            if detail_lower.startswith("downloading"):
                return _status_segment(_STATUS_CLOUD, f"{sat} downloading")
            if any(token in detail_lower for token in ("timed out", "error", "failed", "failure")):
                return _status_segment(_STATUS_CLOUD, f"{sat} failed")
            return _status_segment(_STATUS_CLOUD, f"{sat} {detail}")
        meta = self.cloud_state.meta
        if meta is not None:
            try:
                t = meta.time_utc.strftime("%H:%MZ")
                coverage = self.cloud_state.coverage_ratio
                if coverage is not None and coverage < 0.999:
                    pct = int(round(max(0.0, min(1.0, float(coverage))) * 100.0))
                    return _status_segment(_STATUS_CLOUD, f"{meta.satellite} {pct}% {t}")
                return _status_segment(_STATUS_CLOUD, f"{meta.satellite} {meta.product} {t}")
            except Exception:
                pass
        return _status_segment(_STATUS_CLOUD, f"{sat} idle")

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
            return _status_segment(_STATUS_TERRAIN, str(self.terrain_horizon_state.current_source))
        return ""

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
            return _status_segment(_STATUS_URBAN, str(self.urban_outline_state.current_source))
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
            return ""
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
            return ""
        return _status_segment(_STATUS_SATELLITE, satellite_state.element_epoch_utc.strftime("%H:%MZ"))

    def _jpl_small_body_status_line(self) -> str:
        target = getattr(self.state, "persistent_search_target", None)
        if target is None:
            return ""
        label = str(getattr(target, "label", "")).strip()
        if not label:
            return ""
        altaz_suffix = self._target_altaz_suffix(target)
        next_refresh_utc = getattr(self.state, "persistent_search_next_refresh_utc", None)
        last_error = str(getattr(self.state, "persistent_search_last_error", "")).strip()
        if last_error.casefold() == "none":
            last_error = ""
        if getattr(target, "jpl_group", "") != "sb" and not last_error:
            return f"JPL [{label}]: held{altaz_suffix}"
        if isinstance(next_refresh_utc, datetime):
            refresh_part = next_refresh_utc.strftime("%H:%MZ")
        else:
            refresh_part = "pending"
        if last_error:
            return f"JPL [{label}]: retry {refresh_part} ({last_error}){altaz_suffix}"
        if getattr(target, "jpl_group", "") == "sb":
            return f"JPL [{label}]: retry {refresh_part}{altaz_suffix}"
        return f"JPL [{label}]: refresh {refresh_part}{altaz_suffix}"

    def _target_altaz_suffix(self, target: object) -> str:
        alt = getattr(target, "alt_deg", None)
        az = getattr(target, "az_deg", None)
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
        payload_width = int(
            payload.get("render_width_px", max(2, int(self.client_width())))
        )
        payload_height = int(
            payload.get("render_height_px", max(2, int(self.client_height())))
        )
        current_width = max(2, int(self.client_width()))
        current_height = max(2, int(self.client_height()))
        if (
            payload_generation != current_generation
            or payload_width != current_width
            or payload_height != current_height
        ):
            logger.debug(
                "Discard stale sky payload generation=%s current=%s size=%sx%s current_size=%sx%s",
                payload_generation,
                current_generation,
                payload_width,
                payload_height,
                current_width,
                current_height,
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

        if (
            getattr(self.state, "viewport_interaction_release_pending", False)
        ):
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

        if not self._sky_data_update_timer.isActive():
            self._sky_data_update_timer.start(self.sky_update_interval * 1000)
            if (
                self._clouddisc
                and self.cloud_disc_alpha > 0.0
                and not self._cloud_update_timer.isActive()
            ):
                self._cloud_update_timer.start()
            self.initial_data_loaded.emit()

        if self.state.sky_update_pending and not self._is_shutting_down:
            self.request_sky_data_update(
                self.state.pending_star_vmag_limit,
                reason="pending-follow-up",
            )

        if self.state.cloud_repaint_deferred and not self.state.interaction_mode:
            self.state.cloud_repaint_deferred = False
            self._safe_request_cloud_repaint()

    def request_sky_data_update(
        self,
        star_vmag_limit: Optional[float] = None,
        *,
        reason: str = "manual",
    ) -> None:
        if self.start_background_sky_data_update(
            star_vmag_limit=star_vmag_limit,
            reason=reason,
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
    ) -> bool:
        lat, lon = self.viewer_data.location
        use_lod6_catalog = star_vmag_limit is not None and float(star_vmag_limit) <= 6.0
        star_catalog = self.star_catalog_np
        star_subset_indices = (
            self.star_catalog_lod6_indices if use_lod6_catalog else None
        )
        worker_star_vmag_limit = None if use_lod6_catalog else star_vmag_limit
        started = self._sky_worker.update(
            lat=lat,
            lon=lon,
            observer_height_m=self.viewer_data.observer_height_m,
            view_center=self.viewer_data.view_center,
            star_catalog=star_catalog,
            dso_catalog=self.dso_catalog_np,
            star_vmag_limit=worker_star_vmag_limit,
            star_subset_indices=star_subset_indices,
            delta_t=self.delta_t,
            sky_disc_alpha=self.sky_disc_alpha,
            sky_disc_style=self.sky_disc_style,
            sky_disc_base_size=self.state.sky_disc_base_size,
            edge_fov_deg=float(self.viewer_data.edge_fov_deg),
            content_fov_deg=float(self.content_fov_deg),
            theme=self.theme,
            star_catalog_meta=self.star_catalog_meta,
            render_width_px=max(2, int(self.client_width())),
            render_height_px=max(2, int(self.client_height())),
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

    def refresh_projected_satellite_overlay(self) -> None:
        if float(self.satellite_opacity) <= 0.0:
            return
        validity_remaining_ms = self._satellite_validity_remaining_ms()
        if validity_remaining_ms is not None and validity_remaining_ms <= 0:
            self.state.satellite_overlay_points = None
            self.satellite_state.overlay_points = None
            self.request_client_update()
            self.start_background_satellite_update(reason="time-window-shift")
            return
        records_by_group = self.satellite_state.records_by_group or {}
        if not records_by_group:
            records_by_group = self._load_cached_satellite_records(
                tuple(self._enabled_satellite_groups)
            )
        if not records_by_group:
            self.state.satellite_overlay_points = None
            self.request_client_update()
            return
        lat, lon = self.viewer_data.location
        overlay_points = project_satellite_records(
            records_by_group,
            observer_lat=lat,
            observer_lon=lon,
            observer_height_m=self.viewer_data.observer_height_m,
            time_obj=self._current_time_obj(),
        )
        self.satellite_state.overlay_points = overlay_points
        self.state.satellite_overlay_points = overlay_points
        self.request_client_update()

    def refresh_projected_persistent_search_target(self) -> None:
        target = getattr(self.state, "persistent_search_target", None)
        if target is None:
            return
        if not bool(getattr(target, "persistent_keep_marker", False)):
            return
        if getattr(target, "kind", "") not in {"jpl_small_body", "jpl_body"}:
            return
        projected_altaz = project_jpl_target_altaz_from_state_vector(
            target,
            observer_lat=float(self.viewer_data.location[0]),
            observer_lon=float(self.viewer_data.location[1]),
            observer_height_m=float(self.viewer_data.observer_height_m),
            time_obj=self._current_time_obj(),
        )
        if projected_altaz is None:
            _jpl_debug_print(
                "refresh-project-none "
                f"label={target.label} command={target.command} "
                f"target_time_utc={getattr(target, 'target_time_utc', None)!r} "
                f"epoch={getattr(target, 'horizons_epoch_utc', None)!r} "
                f"pos={getattr(target, 'horizons_position_km', None)!r} "
                f"vel={getattr(target, 'horizons_velocity_km_s', None)!r}"
            )
            return
        alt_deg, az_deg = projected_altaz
        _jpl_debug_print(
            "refresh "
            f"label={target.label} command={target.command} "
            f"target_time_utc={getattr(target, 'target_time_utc', None)!r} "
            f"epoch={getattr(target, 'horizons_epoch_utc', None)!r} "
            f"pos={getattr(target, 'horizons_position_km', None)!r} "
            f"vel={getattr(target, 'horizons_velocity_km_s', None)!r} "
            f"projected_alt={float(alt_deg):.3f} projected_az={float(az_deg) % 360.0:.3f}"
        )
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
            overlay_points=payload.get("overlay_points"),
            element_epoch_utc=element_epoch,
            refreshed_at_utc=payload.get("refreshed_at_utc"),
        )
        if banner:
            self.satellite_state.set_banner(banner)
        if float(self.satellite_opacity) > 0.0:
            self.state.satellite_overlay_points = payload.get("overlay_points")
            self._schedule_next_satellite_refresh()
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
        self._compositor.invalidate()
        if self.state.interaction_mode:
            self.state.cloud_repaint_deferred = True
            return
        self._safe_request_cloud_repaint()

    def _on_cloud_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.cloud_state.set_error_banner(banner)
        if self.state.interaction_mode:
            self.state.cloud_repaint_deferred = True
        if banner:
            self._safe_request_cloud_repaint()

    def start_background_terrain_horizon_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down:
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

    def start_background_urban_outline_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down:
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

    def refresh_projected_aircraft_overlay(self) -> None:
        if float(self.aircraft_opacity) <= 0.0:
            return
        snapshots = self.aircraft_state.snapshots
        if not snapshots:
            self.state.aircraft_overlay_points = None
            self.request_client_update()
            return
        lat, lon = self.viewer_data.location
        overlay_points = project_aircraft_snapshots(
            snapshots,
            observer_lat=lat,
            observer_lon=lon,
            observer_height_m=self.viewer_data.observer_height_m,
            time_obj=self._current_time_obj(),
        )
        self.aircraft_state.overlay_points = overlay_points
        self.state.aircraft_overlay_points = overlay_points
        self.request_client_update()

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
            overlay_points=payload.get("overlay_points"),
            bbox=payload.get("bbox"),
            refreshed_at_utc=refreshed_at,
        )
        if banner:
            self.aircraft_state.set_banner(banner)
        if float(self.aircraft_opacity) > 0.0:
            self.state.aircraft_overlay_points = payload.get("overlay_points")
            self._schedule_next_aircraft_refresh()
        self.request_client_update()
        self._maybe_save_aircraft_debug_snapshot(payload)

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
        self.state.terrain_horizon_profile = payload["profile_altaz"]
        self.state.terrain_horizon_profile_distances_m = payload.get("profile_distances_m")
        self.state.terrain_horizon_secondary_profile_altaz_layers = payload.get(
            "secondary_profile_altaz_layers"
        )
        self.state.terrain_horizon_secondary_profile_distances_m_layers = payload.get(
            "secondary_profile_distances_m_layers"
        )
        self._compositor.invalidate()
        self.request_client_update()

    def _on_terrain_horizon_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        self.terrain_horizon_state.clear_profile()
        self.state.terrain_horizon_profile = None
        self.state.terrain_horizon_profile_distances_m = None
        self.state.terrain_horizon_secondary_profile_altaz_layers = None
        self.state.terrain_horizon_secondary_profile_distances_m_layers = None
        if banner:
            self.terrain_horizon_state.set_error_banner(banner)
        self._compositor.invalidate()
        self.request_client_update()

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
        self._compositor.invalidate()
        self.request_client_update()

    def _on_urban_outline_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        self.urban_outline_state.clear_outlines()
        self.state.urban_outlines = None
        if banner:
            self.urban_outline_state.set_error_banner(banner)
        self._compositor.invalidate()
        self.request_client_update()
