from __future__ import annotations

import logging
from datetime import datetime, timedelta, timezone
from typing import cast

from ..geosatellite.pipeline import is_within_europe_band
from ..overlay_time import overlay_availability_for_delta
from ..paths import CLOUD_UPDATE_INTERVAL
from .window_update_common import (
    _CloudProjectionUpdateOwner,
    _request_cloud_projection_update,
)

logger = logging.getLogger(__name__)


class SkyWindowCloudUpdatesMixin:
    def request_cloud_projection_update(self, *, reason: str) -> None:
        _request_cloud_projection_update(
            cast(_CloudProjectionUpdateOwner, self), reason=reason
        )

    def _geo_satellite_mode_active(self) -> bool:
        return bool(
            self._geo_satellite_enabled
            and self._geosatellite_controller is not None
            and float(self.cloud_disc_alpha) > 0.0
            and overlay_availability_for_delta(self.delta_t).cloud
            and is_within_europe_band(
                float(self.viewer_data.lat_deg), float(self.viewer_data.lon_deg)
            )
        )

    def _active_cloud_state(self):
        if self._geo_satellite_mode_active():
            return self.geosatellite_state
        return self.cloud_state

    def _cloud_projection_next_refresh_delay_ms(self) -> int | None:
        next_refresh_utc = self.state.cloud_projection_next_refresh_utc
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

    def _safe_request_cloud_repaint(self) -> None:
        """Best-effort repaint request; ignores teardown-time signal errors."""
        if self._is_shutting_down:
            return
        try:
            self.cloud_repaint_requested.emit()
        except RuntimeError:
            logger.debug("Skip cloud repaint emit during shutdown.")

    def _cloud_layer_enabled(self) -> bool:
        if float(self.cloud_disc_alpha) <= 0.0:
            return False
        if self._geo_satellite_mode_active():
            return self._geosatellite_controller is not None
        return self._clouddisc is not None and self._cloud_controller is not None

    def start_background_cloud_update(self, reason: str = "manual") -> None:
        if self._is_shutting_down:
            return
        if self._viewport_interaction_active():
            return
        if self._geo_satellite_mode_active():
            self.start_background_geo_satellite_update(reason=reason)
            return
        if not (self._cloud_controller and self._cloud_layer_enabled()):
            return
        lat, lon = self.viewer_data.location
        self._cloud_controller.update_source(
            lat=lat,
            lon=lon,
            reason=reason,
        )

    def start_background_geo_satellite_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down:
            return False
        if self._viewport_interaction_active():
            return False
        if not self._geo_satellite_mode_active():
            return False
        if self._geosatellite_controller is None:
            return False
        lat, lon = self.viewer_data.location
        return self._geosatellite_controller.update(
            observer_lat=lat,
            observer_lon=lon,
            alt=float(self.viewer_data.view_alt_deg),
            az=float(self.viewer_data.view_az_deg),
            fov_deg=float(self.viewer_data.edge_fov_deg) + 2.0,
            render_generation=int(self._disc_generation),
            reason=reason,
        )

    def reproject_geo_satellite_overlay(self, reason: str = "manual") -> None:
        self.start_background_geo_satellite_update(reason=reason)

    def _start_cloud_projection_update(self, reason: str = "manual") -> bool:
        if self._is_shutting_down:
            return False
        if self._viewport_interaction_active():
            return False
        if self._geo_satellite_mode_active():
            return self.start_background_geo_satellite_update(reason=reason)
        if not (self._cloud_controller and self._cloud_layer_enabled()):
            return False
        if not self._cloud_controller.has_source_data():
            return False
        lat, lon = self.viewer_data.location
        alt, az = self.viewer_data.view_center
        return self._cloud_controller.update_render(
            lat=lat,
            lon=lon,
            alt=alt,
            az=az,
            radius_px=self.state.cloud_base_size,
            content_fov_deg=float(self.content_fov_deg),
            reason=reason,
            render_generation=int(self._disc_generation),
        )

    def reproject_cloud_overlay(self, reason: str = "manual") -> None:
        if self._is_shutting_down:
            return
        if self._viewport_interaction_active():
            return
        if self._geo_satellite_mode_active():
            self.start_background_geo_satellite_update(reason=reason)
            return
        if not (self._cloud_controller and self._cloud_layer_enabled()):
            return
        if not self._cloud_controller.has_source_data():
            self.state.cloud_projection_next_refresh_utc = None
            self.start_background_cloud_update(reason=reason)
            return
        self.state.cloud_projection_next_refresh_utc = datetime.now(timezone.utc)
        self._on_scheduler_tick()

    def _on_geosatellite_started(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.geosatellite_state.set_banner(banner)
        self.request_client_update()

    def _on_geosatellite_source_ready(self, payload: dict) -> None:
        refreshed_at = payload.get("refreshed_at_utc")
        if not isinstance(refreshed_at, datetime):
            refreshed_at = datetime.now(timezone.utc)
        banner = str(payload.get("banner", "")).strip()
        self.geosatellite_state.set_source_ready(
            refreshed_at_utc=refreshed_at,
            banner_text=banner or "Geo-sat + Projecting",
            current_source="Geo-sat",
        )
        self.state.cloud_projection_next_refresh_utc = None

    def _on_geosatellite_ready(self, payload: dict) -> None:
        current_generation = int(self._disc_generation)
        payload_generation = int(payload.get("render_generation", current_generation))
        if payload_generation != current_generation and not self._is_shutting_down:
            logger.debug(
                "Discard stale Geo-satellite payload generation=%s current=%s",
                payload_generation,
                current_generation,
            )
            return
        captured_at_utc = payload.get("captured_at_utc")
        if not isinstance(captured_at_utc, datetime):
            captured_at_utc = payload.get("time_utc")
        if not isinstance(captured_at_utc, datetime):
            captured_at_utc = datetime.now(timezone.utc)
        fetched_at_utc = payload.get("fetched_at_utc")
        if not isinstance(fetched_at_utc, datetime):
            fetched_at_utc = captured_at_utc
        self.geosatellite_state.set_result(
            payload["image"],
            payload.get("meta"),
            az=float(payload["az"]),
            time_utc=captured_at_utc,
            cloud_amount_field=payload.get("cloud_amount_field"),
            altaz_grid=payload.get("altaz_grid"),
            missing_mask=payload.get("missing_mask"),
            coverage_ratio=payload.get("coverage_ratio"),
            source_key=payload.get("source_key"),
            request_id=payload.get("request_id"),
            current_source="Geo-sat",
            captured_at_utc=captured_at_utc,
            fetched_at_utc=fetched_at_utc,
        )
        refreshed_at_utc = payload.get("refreshed_at_utc")
        if not isinstance(refreshed_at_utc, datetime):
            refreshed_at_utc = datetime.now(timezone.utc)
        self.state.cloud_next_refresh_utc = refreshed_at_utc + timedelta(
            seconds=CLOUD_UPDATE_INTERVAL
        )
        self.state.cloud_projection_next_refresh_utc = None
        self._compositor.invalidate()
        if self.state.interaction_mode:
            self.state.cloud_repaint_deferred = True
            return
        self._safe_request_cloud_repaint()

    def _on_geosatellite_failed(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.geosatellite_state.set_error_banner(banner)
        self.state.cloud_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=CLOUD_UPDATE_INTERVAL
        )
        self.state.cloud_projection_next_refresh_utc = None
        if self.state.interaction_mode:
            self.state.cloud_repaint_deferred = True
        if banner:
            self._safe_request_cloud_repaint()

    def _on_cloud_started(self, payload: dict) -> None:
        sat = str(payload.get("satellite", "")).strip()
        banner = str(payload.get("banner", "")).strip()
        if sat:
            self.cloud_state.current_satellite = sat
        if banner:
            self.cloud_state.set_error_banner(banner)
        if sat or banner:
            self.request_client_update()

    def _on_cloud_source_ready(self, payload: dict) -> None:
        sat = str(payload.get("satellite", "")).strip()
        source_key = payload.get("source_key")
        refreshed_at = payload.get("refreshed_at_utc")
        if not isinstance(refreshed_at, datetime):
            refreshed_at = datetime.now(timezone.utc)
        banner = str(payload.get("banner", "")).strip()
        self.cloud_state.set_source_ready(
            refreshed_at_utc=refreshed_at,
            satellite=sat or None,
            source_key=source_key,
            banner_text=banner or "Clouds: projecting...",
            altaz_grid=payload.get("altaz_grid"),
        )
        self.state.cloud_next_refresh_utc = refreshed_at + timedelta(
            seconds=CLOUD_UPDATE_INTERVAL
        )
        self.state.cloud_projection_next_refresh_utc = datetime.now(timezone.utc)
        if self.state.interaction_mode:
            self.state.cloud_repaint_deferred = True
        else:
            _request_cloud_projection_update(
                cast(_CloudProjectionUpdateOwner, self), reason="source-ready"
            )

    def _on_cloud_ready(self, payload: dict) -> None:
        current_generation = int(self._disc_generation)
        payload_generation = int(payload.get("render_generation", current_generation))
        if payload_generation != current_generation:
            logger.debug(
                "Discard stale cloud payload generation=%s current=%s",
                payload_generation,
                current_generation,
            )
            if not self._is_shutting_down:
                _request_cloud_projection_update(
                    cast(_CloudProjectionUpdateOwner, self), reason="stale-render"
                )
            return
        self.cloud_state.set_result(
            payload["image"],
            payload.get("meta"),
            az=float(payload["az"]),
            time_utc=payload["time_utc"],
            cloud_amount_field=payload.get("cloud_amount_field"),
            altaz_grid=payload.get("altaz_grid"),
            missing_mask=payload.get("missing_mask"),
            coverage_ratio=payload.get("coverage_ratio"),
            source_key=payload.get("source_key"),
            request_id=payload.get("request_id"),
            source_expected_count=payload.get("source_expected_count"),
            source_available_count=payload.get("source_available_count"),
            source_completeness_ratio=payload.get("source_completeness_ratio"),
        )
        self.state.cloud_projection_next_refresh_utc = None
        self._compositor.invalidate()
        if self.state.interaction_mode:
            self.state.cloud_repaint_deferred = True
            return
        self._safe_request_cloud_repaint()

    def _on_cloud_failed(self, payload: dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        if banner:
            self.cloud_state.set_error_banner(banner)
        self.state.cloud_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=CLOUD_UPDATE_INTERVAL
        )
        self.state.cloud_projection_next_refresh_utc = None
        if self.state.interaction_mode:
            self.state.cloud_repaint_deferred = True
        if banner:
            self._safe_request_cloud_repaint()
