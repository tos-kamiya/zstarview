from __future__ import annotations

import logging
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional

from ..aircraft import project_aircraft_snapshots
from ..paths import CACHE_PATH
from ..satellite_constants import SATELLITE_ISS_CACHE_KEY
from ..satellites import project_satellite_records

logger = logging.getLogger(__name__)


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
        render_current_image = getattr(self, "render_current_image", None)
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
                logger.warning("Failed to save aircraft debug snapshot: %s", output_path)
                return
            logger.info("Saved aircraft debug snapshot: %s", output_path)
        except Exception as exc:
            logger.warning("Aircraft debug snapshot failed: %s", exc, exc_info=True)

    def _status_line_message(self) -> str:
        parts: list[str] = []
        cloud_message = self._cloud_status_line()
        if cloud_message:
            parts.append(cloud_message)
        satellite_status_line = getattr(self, "_satellite_status_line", None)
        if callable(satellite_status_line):
            satellite_message = satellite_status_line()
            if satellite_message:
                parts.append(satellite_message)
        aircraft_status_line = getattr(self, "_aircraft_status_line", None)
        if callable(aircraft_status_line):
            aircraft_message = aircraft_status_line()
            if aircraft_message:
                parts.append(aircraft_message)
        terrain_message = self._terrain_horizon_status_line()
        if terrain_message:
            parts.append(terrain_message)
        urban_message = self._urban_outline_status_line()
        if urban_message:
            parts.append(urban_message)
        return " | ".join(parts)

    def _safe_request_cloud_repaint(self) -> None:
        """Best-effort repaint request; ignores teardown-time signal errors."""
        if self._is_shutting_down:
            return
        try:
            self.cloud_repaint_requested.emit()
        except RuntimeError:
            logger.debug("Skip cloud repaint emit during shutdown.")

    def _cloud_status_line(self) -> str:
        cloud_disc_alpha = float(getattr(self, "cloud_disc_alpha", 1.0))
        if cloud_disc_alpha <= 0.0 and not self.cloud_state.banner_text:
            return ""
        sat = self.cloud_state.current_satellite or self._predicted_cloud_satellite()
        if self.cloud_state.banner_text:
            detail = self.cloud_state.banner_text.removeprefix("Clouds:").strip()
            if "download" in detail.lower():
                return f"Clouds [{sat}]: downloading"
            return f"Clouds [{sat}]: {detail}"
        meta = self.cloud_state.meta
        if meta is not None:
            try:
                t = meta.time_utc.strftime("%H:%MZ")
                coverage = self.cloud_state.coverage_ratio
                if coverage is not None and coverage < 0.999:
                    pct = int(round(max(0.0, min(1.0, float(coverage))) * 100.0))
                    return f"Clouds [{meta.satellite}]: partial {pct}% {t}"
                return f"Clouds [{meta.satellite}]: {meta.product} {t}"
            except Exception:
                pass
        return f"Clouds [{sat}]: idle"

    def _terrain_horizon_status_line(self) -> str:
        if self.terrain_horizon_opacity <= 0.0 and not self.terrain_horizon_state.banner_text:
            return ""
        if self.terrain_horizon_state.banner_text:
            return self.terrain_horizon_state.banner_text
        return ""

    def _urban_outline_status_line(self) -> str:
        if self.urban_outline_opacity <= 0.0 and not self.urban_outline_state.banner_text:
            return ""
        if self.urban_outline_state.banner_text:
            return self.urban_outline_state.banner_text
        if self.urban_outline_state.current_source:
            return f"Urban outline: {self.urban_outline_state.current_source}"
        return ""

    def _aircraft_status_line(self) -> str:
        if float(getattr(self, "aircraft_opacity", 0.0)) <= 0.0:
            return ""
        aircraft_state = getattr(self, "aircraft_state", None)
        if aircraft_state is None:
            return ""
        if aircraft_state.banner_text:
            return aircraft_state.banner_text
        if aircraft_state.last_success_utc is None:
            return ""
        return f"Aircraft: {aircraft_state.last_success_utc.strftime('%H:%MZ')}"

    def _satellite_status_line(self) -> str:
        if float(getattr(self, "satellite_opacity", 0.0)) <= 0.0:
            return ""
        satellite_state = getattr(self, "satellite_state", None)
        if satellite_state is None:
            return ""
        if satellite_state.banner_text:
            return satellite_state.banner_text
        if satellite_state.element_epoch_utc is None:
            return ""
        return f"Satellites: {satellite_state.element_epoch_utc.strftime('%H:%MZ')}"

    def _on_sky_data_calculated(self, payload: Dict) -> None:
        current_generation = int(getattr(self, "_disc_generation", 0))
        payload_generation = int(payload.get("render_generation", current_generation))
        payload_width = int(payload.get("render_width_px", max(2, int(self.width()))))
        payload_height = int(payload.get("render_height_px", max(2, int(self.height()))))
        current_width = max(2, int(self.width()))
        current_height = max(2, int(self.height()))
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
                self.request_sky_data_update()
            return
        if not self.state.viewport_interaction_mode:
            self.state.render_view_center = tuple(
                payload.get("view_center", self.viewer_data.view_center)
            )
        self.state.celestial_data = payload["celestial"]
        self.state.sky_disc_image = payload["sky_disc"]

        self._compositor.invalidate()
        self.update()

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
            )

        if self.state.cloud_repaint_deferred and not self.state.interaction_mode:
            self.state.cloud_repaint_deferred = False
            self._safe_request_cloud_repaint()

    def request_sky_data_update(
        self,
        star_vmag_limit: Optional[float] = None,
    ) -> None:
        if self.start_background_sky_data_update(
            star_vmag_limit=star_vmag_limit,
        ):
            self.state.sky_update_pending = False
            self.state.pending_star_vmag_limit = None
            return
        self.state.sky_update_pending = True
        self.state.pending_star_vmag_limit = star_vmag_limit
        logger.debug("Sky data update deferred; worker is busy.")

    def start_background_sky_data_update(
        self,
        is_initial_load: bool = False,
        star_vmag_limit: Optional[float] = None,
    ) -> bool:
        lat, lon = self.viewer_data.location
        use_lod6_catalog = star_vmag_limit is not None and float(star_vmag_limit) <= 6.0
        star_catalog = self.star_catalog_np
        star_subset_indices = self.star_catalog_lod6_indices if use_lod6_catalog else None
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
            sky_disc_base_size=self.state.sky_disc_base_size,
            content_fov_deg=float(self.content_fov_deg),
            visual_preset=self.visual_preset,
            star_catalog_meta=self.star_catalog_meta,
            render_width_px=max(2, int(self.width())),
            render_height_px=max(2, int(self.height())),
            render_generation=int(getattr(self, "_disc_generation", 0)),
        )
        if started:
            if is_initial_load:
                logger.info("Calculating initial sky data...")
            else:
                logger.info("Updating sky data...")
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
            render_generation=int(getattr(self, "_disc_generation", 0)),
        )

    def start_background_satellite_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down:
            return False
        if float(getattr(self, "satellite_opacity", 0.0)) <= 0.0:
            return False
        if self._satellite_controller is None:
            return False
        lat, lon = self.viewer_data.location
        return self._satellite_controller.update(
            observer_lat=lat,
            observer_lon=lon,
            observer_height_m=self.viewer_data.observer_height_m,
            time_obj=self._current_time_obj(),
            enabled_groups=tuple(getattr(self, "_enabled_satellite_groups", (SATELLITE_ISS_CACHE_KEY,))),
            reason=reason,
        )

    def refresh_projected_satellite_overlay(self) -> None:
        if float(getattr(self, "satellite_opacity", 0.0)) <= 0.0:
            return
        validity_remaining_ms = getattr(self, "_satellite_validity_remaining_ms", lambda: None)()
        if validity_remaining_ms is not None and validity_remaining_ms <= 0:
            self.state.satellite_overlay_points = None
            self.satellite_state.overlay_points = None
            self.update()
            self.start_background_satellite_update(reason="time-window-shift")
            return
        records_by_group = getattr(self.satellite_state, "records_by_group", None) or {}
        if not records_by_group:
            load_cached_records = getattr(self, "_load_cached_satellite_records", None)
            if callable(load_cached_records):
                records_by_group = load_cached_records(
                    tuple(getattr(self, "_enabled_satellite_groups", (SATELLITE_ISS_CACHE_KEY,)))
                )
        if not records_by_group:
            self.state.satellite_overlay_points = None
            self.update()
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
        self.update()

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
            self.update()

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
        if float(getattr(self, "satellite_opacity", 0.0)) > 0.0:
            self.state.satellite_overlay_points = payload.get("overlay_points")
            schedule_next = getattr(self, "_schedule_next_satellite_refresh", None)
            if callable(schedule_next):
                schedule_next()
        self.update()

    def _on_satellite_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.satellite_state.set_error_banner(banner)
        schedule_retry = getattr(self, "_schedule_satellite_retry_after_failure", None)
        if callable(schedule_retry) and float(getattr(self, "satellite_opacity", 0.0)) > 0.0:
            schedule_retry()
        self.update()

    def _on_cloud_ready(self, payload: Dict) -> None:
        current_generation = int(getattr(self, "_disc_generation", 0))
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
        self.cloud_state.image = None
        self.cloud_state.missing_mask = None
        self.cloud_state.cloud_amount_field = None
        self._compositor.invalidate()
        if banner:
            self.cloud_state.set_error_banner(banner)
        if self.state.interaction_mode:
            self.state.cloud_repaint_deferred = True
            return
        self._safe_request_cloud_repaint()

    def start_background_terrain_horizon_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down:
            return False
        if self.terrain_horizon_opacity <= 0.0 or self._terrain_horizon_controller is None:
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
        if float(getattr(self, "aircraft_opacity", 0.0)) <= 0.0:
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
        if float(getattr(self, "aircraft_opacity", 0.0)) <= 0.0:
            return
        snapshots = self.aircraft_state.snapshots
        if not snapshots:
            self.state.aircraft_overlay_points = None
            self.update()
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
        self.update()

    def _on_aircraft_started(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.aircraft_state.set_banner(banner)
            self.update()

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
        if float(getattr(self, "aircraft_opacity", 0.0)) > 0.0:
            self.state.aircraft_overlay_points = payload.get("overlay_points")
            schedule_next = getattr(self, "_schedule_next_aircraft_refresh", None)
            if callable(schedule_next):
                schedule_next()
        self.update()
        self._maybe_save_aircraft_debug_snapshot(payload)

    def _on_aircraft_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.aircraft_state.set_error_banner(banner)
        schedule_next = getattr(self, "_schedule_next_aircraft_refresh", None)
        if callable(schedule_next) and float(getattr(self, "aircraft_opacity", 0.0)) > 0.0:
            schedule_next()
        self.update()

    def _on_terrain_horizon_started(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.terrain_horizon_state.banner_text = banner
        self.update()

    def _on_terrain_horizon_ready(self, payload: Dict) -> None:
        self.terrain_horizon_state.set_result(
            payload["profile_altaz"],
            source=str(payload.get("source", "")).strip(),
        )
        self.state.terrain_horizon_profile = payload["profile_altaz"]
        self._compositor.invalidate()
        self.update()

    def _on_terrain_horizon_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        self.terrain_horizon_state.clear_profile()
        self.state.terrain_horizon_profile = None
        if banner:
            self.terrain_horizon_state.set_error_banner(banner)
        self._compositor.invalidate()
        self.update()

    def _on_urban_outline_started(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.urban_outline_state.banner_text = banner
        self.update()

    def _on_urban_outline_ready(self, payload: Dict) -> None:
        outlines = payload.get("outlines")
        self.urban_outline_state.set_result(
            outlines,
            source=str(payload.get("source", "")).strip() or "ready",
        )
        self.state.urban_outlines = outlines
        self._compositor.invalidate()
        self.update()

    def _on_urban_outline_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        self.urban_outline_state.clear_outlines()
        self.state.urban_outlines = None
        if banner:
            self.urban_outline_state.set_error_banner(banner)
        self._compositor.invalidate()
        self.update()
