from __future__ import annotations

import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)


def _get_state_value(obj, name: str, *, legacy_name: str, default=None):
    state = getattr(obj, "state", None)
    if state is not None and hasattr(state, name):
        return getattr(state, name)
    return getattr(obj, legacy_name, default)


def _set_state_value(obj, name: str, value, *, legacy_name: str) -> None:
    state = getattr(obj, "state", None)
    if state is not None and hasattr(state, name):
        setattr(state, name, value)
    setattr(obj, legacy_name, value)


class SkyWindowUpdatesMixin:
    def _safe_request_cloud_repaint(self) -> None:
        """Best-effort repaint request; ignores teardown-time signal errors."""
        if self._is_shutting_down:
            return
        try:
            self.cloud_repaint_requested.emit()
        except RuntimeError:
            logger.debug("Skip cloud repaint emit during shutdown.")

    def _cloud_status_line(self) -> str:
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

    def _on_sky_data_calculated(self, payload: Dict) -> None:
        _set_state_value(
            self,
            "render_view_center",
            tuple(payload.get("view_center", self.viewer_data.view_center)),
            legacy_name="_render_view_center",
        )
        _set_state_value(
            self,
            "celestial_data",
            payload["celestial"],
            legacy_name="celestial_data",
        )
        _set_state_value(
            self,
            "sky_disc_image",
            payload["sky_disc"],
            legacy_name="_sky_disc_image",
        )

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

        if _get_state_value(
            self,
            "sky_update_pending",
            legacy_name="_sky_update_pending",
            default=False,
        ) and not self._is_shutting_down:
            self.request_sky_data_update(
                _get_state_value(
                    self,
                    "pending_star_vmag_limit",
                    legacy_name="_pending_star_vmag_limit",
                ),
            )

        if _get_state_value(
            self,
            "cloud_repaint_deferred",
            legacy_name="_cloud_repaint_deferred",
            default=False,
        ) and not _get_state_value(
            self,
            "interaction_mode",
            legacy_name="_interaction_mode",
            default=False,
        ):
            _set_state_value(
                self,
                "cloud_repaint_deferred",
                False,
                legacy_name="_cloud_repaint_deferred",
            )
            self._safe_request_cloud_repaint()

    def request_sky_data_update(
        self,
        star_vmag_limit: Optional[float] = None,
    ) -> None:
        if self.start_background_sky_data_update(
            star_vmag_limit=star_vmag_limit,
        ):
            _set_state_value(
                self,
                "sky_update_pending",
                False,
                legacy_name="_sky_update_pending",
            )
            _set_state_value(
                self,
                "pending_star_vmag_limit",
                None,
                legacy_name="_pending_star_vmag_limit",
            )
            return
        _set_state_value(
            self,
            "sky_update_pending",
            True,
            legacy_name="_sky_update_pending",
        )
        _set_state_value(
            self,
            "pending_star_vmag_limit",
            star_vmag_limit,
            legacy_name="_pending_star_vmag_limit",
        )
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
            view_center=self.viewer_data.view_center,
            star_catalog=star_catalog,
            dso_catalog=self.dso_catalog_np,
            star_vmag_limit=worker_star_vmag_limit,
            delta_t=self.delta_t,
            sky_disc_alpha=self.sky_disc_alpha,
            sky_disc_base_size=_get_state_value(
                self,
                "sky_disc_base_size",
                legacy_name="_sky_disc_base_size",
                default=1024,
            ),
            debug_eclipses=self._debug_eclipses,
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
            radius_px=_get_state_value(
                self,
                "cloud_base_size",
                legacy_name="_cloud_base_size",
                default=256,
            ),
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
        if _get_state_value(
            self,
            "interaction_mode",
            legacy_name="_interaction_mode",
            default=False,
        ):
            _set_state_value(
                self,
                "cloud_repaint_deferred",
                True,
                legacy_name="_cloud_repaint_deferred",
            )
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
        if _get_state_value(
            self,
            "interaction_mode",
            legacy_name="_interaction_mode",
            default=False,
        ):
            _set_state_value(
                self,
                "cloud_repaint_deferred",
                True,
                legacy_name="_cloud_repaint_deferred",
            )
            return
        self._safe_request_cloud_repaint()
