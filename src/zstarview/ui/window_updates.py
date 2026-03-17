from __future__ import annotations

import logging
from datetime import datetime, timezone
from typing import Dict, Optional

logger = logging.getLogger(__name__)


class SkyWindowUpdatesMixin:
    def _status_line_message(self) -> str:
        parts: list[str] = []
        cloud_message = self._cloud_status_line()
        if cloud_message:
            parts.append(cloud_message)
        terrain_message = self._terrain_horizon_status_line()
        if terrain_message:
            parts.append(terrain_message)
        urban_message = self._urban_outline_status_line()
        if urban_message:
            parts.append(urban_message)
        aircraft_message = self._aircraft_status_line()
        if aircraft_message:
            parts.append(aircraft_message)
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
        aircraft_state = getattr(self, "aircraft_state", None)
        if aircraft_state is None:
            return ""
        if aircraft_state.banner_text:
            return aircraft_state.banner_text
        if aircraft_state.last_success_utc is None:
            return ""
        return f"Aircraft: {aircraft_state.last_success_utc.strftime('%H:%MZ')}"

    def _on_sky_data_calculated(self, payload: Dict) -> None:
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
        star_catalog = self.star_catalog_lod6_np if use_lod6_catalog else self.star_catalog_full_np
        worker_star_vmag_limit = None if use_lod6_catalog else star_vmag_limit
        started = self._sky_worker.update(
            lat=lat,
            lon=lon,
            observer_height_m=self.viewer_data.observer_height_m,
            view_center=self.viewer_data.view_center,
            star_catalog=star_catalog,
            dso_catalog=self.dso_catalog_np,
            star_vmag_limit=worker_star_vmag_limit,
            delta_t=self.delta_t,
            sky_disc_alpha=self.sky_disc_alpha,
            sky_disc_base_size=self.state.sky_disc_base_size,
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
            reason=reason,
        )

    def _on_cloud_started(self, payload: Dict) -> None:
        sat = str(payload.get("satellite", "")).strip()
        banner = str(payload.get("banner", "")).strip()
        if sat:
            self.cloud_state.current_satellite = sat
        if banner:
            self.cloud_state.set_error_banner(banner)

    def _on_cloud_ready(self, payload: Dict) -> None:
        self.cloud_state.set_result(
            payload["image"],
            payload.get("meta"),
            az=float(payload["az"]),
            time_utc=payload["time_utc"],
            stripe_density=payload.get("stripe_density"),
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
        self.cloud_state.stripe_density = None
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

    def _on_aircraft_started(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.aircraft_state.set_banner(banner)
            self.update()

    def _on_aircraft_ready(self, payload: Dict) -> None:
        refreshed_at = payload.get("refreshed_at_utc")
        if not isinstance(refreshed_at, datetime):
            refreshed_at = datetime.now(timezone.utc)
        self.aircraft_state.set_result(
            payload.get("snapshots", []),
            overlay_points=payload.get("overlay_points"),
            bbox=payload.get("bbox"),
            refreshed_at_utc=refreshed_at,
        )
        self.state.aircraft_overlay_points = payload.get("overlay_points")
        self.update()

    def _on_aircraft_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.aircraft_state.set_error_banner(banner)
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
