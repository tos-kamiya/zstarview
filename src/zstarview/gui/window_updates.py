from __future__ import annotations

import logging
import time
from dataclasses import replace
from datetime import datetime, timedelta, timezone
from typing import Protocol, cast

from ..aircraft_constants import AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS
from ..astro import load_ephemeris
from ..clouddisc.providers.select import GOES_SATELLITES
from ..geosatellite.pipeline import is_within_europe_band
from ..night_lights import is_night_light_enabled
from ..overlay_time import overlay_availability_for_delta
from ..paths import CLOUD_UPDATE_INTERVAL
from ..render import geometry as render_geometry
from ..satellite_constants import SATELLITE_POSITION_REFRESH_INTERVAL_SECONDS
from ..search.jpl import project_jpl_target_altaz_from_state_vector

logger = logging.getLogger(__name__)
_STATUS_CLOUD = "☁"
_STATUS_WATER = "W"
_STATUS_SATELLITE = "🛰"
_STATUS_AIRCRAFT = "✈"
_STATUS_TROPICAL_CYCLONE = "TC"
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
        clean = clean[len(prefix) :].strip()
    return clean


def _urban_outline_source_name(source_value: object) -> str:
    source = str(source_value or "").strip()
    normalized = source.casefold()
    if normalized == "urban: plateau":
        return "PLATEAU"
    if normalized in {"urban: cache", "urban: overture maps"}:
        return "Overture Maps"
    return _strip_status_prefix(source, "Urban:")


def _cloud_satellite_group(satellite: str) -> str:
    sat = str(satellite).strip()
    if sat in GOES_SATELLITES:
        return "GOES"
    if sat == "HIMAWARI":
        return "HIMAWARI"
    return sat


def _cloud_source_label(satellite: str) -> str:
    sat = str(satellite).strip()
    group = _cloud_satellite_group(sat)
    if group == sat:
        return sat
    return f"{group} {sat}"


class _CloudProjectionUpdateOwner(Protocol):
    def _geo_satellite_mode_active(self) -> bool: ...

    def reproject_geo_satellite_overlay(self, reason: str = "manual") -> None: ...

    def reproject_cloud_overlay(self, reason: str = "manual") -> None: ...


class _ViewportInteractionState(Protocol):
    viewport_interaction_release_pending: bool
    viewport_interaction_completion_reason: str | None
    viewport_interaction_mode: bool
    viewport_interaction_stars: object | None


class _ViewportInteractionWaitOwner(Protocol):
    state: _ViewportInteractionState

    def _sync_viewport_interaction_chrome_visibility(self) -> None: ...

    def request_client_update(self) -> None: ...


def _request_cloud_projection_update(
    obj: _CloudProjectionUpdateOwner, *, reason: str
) -> None:
    if obj._geo_satellite_mode_active():
        obj.reproject_geo_satellite_overlay(reason=reason)
        return
    obj.reproject_cloud_overlay(reason=reason)


def _initial_data_load_active(obj: object) -> bool:
    return bool(obj._startup_initial_load_started) and not bool(
        obj._startup_initial_data_loaded
    )


def _extract_sun_altitude_deg(celestial_data: object) -> float | None:
    planets = getattr(celestial_data, "planets", None)
    if not isinstance(planets, (list, tuple)):
        return None
    for body in planets:
        if getattr(body, "name", "").strip().casefold() == "sun":
            alt = getattr(body, "alt", None)
            if isinstance(alt, (int, float)):
                return float(alt)
            try:
                return float(alt)
            except (TypeError, ValueError):
                return None
    return None


def _startup_night_light_requires_warmup(obj: object, payload: dict) -> bool:
    if float(getattr(obj, "terrain_horizon_opacity", 0.0)) <= 0.0:
        return False
    if (
        float(getattr(obj, "night_light_opacity", 0.0)) <= 0.0
        and float(getattr(obj, "ridge_glow_opacity", 0.0)) <= 0.0
    ):
        return False
    celestial_data = payload.get("celestial")
    sun_alt_deg = _extract_sun_altitude_deg(celestial_data)
    return sun_alt_deg is not None and is_night_light_enabled(float(sun_alt_deg))


def _clear_viewport_interaction_wait(obj: _ViewportInteractionWaitOwner) -> None:
    state = obj.state
    state.viewport_interaction_release_pending = False
    state.viewport_interaction_completion_reason = None
    state.viewport_interaction_mode = False
    state.viewport_interaction_stars = None
    obj._sync_viewport_interaction_chrome_visibility()
    obj.request_client_update()


class SkyWindowUpdatesMixin:
    def _viewport_interaction_active(self) -> bool:
        return bool(self.state.viewport_interaction_mode)

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

    def _background_updates_busy(self) -> bool:
        controllers = (
            self._sky_worker,
            self._cloud_controller,
            self._geosatellite_controller,
            self._satellite_controller,
            self._aircraft_controller,
            self._tropical_cyclone_controller,
            self._jpl_small_body_controller,
            self._terrain_horizon_controller,
            self._water_overlay_controller,
            self._urban_outline_controller,
        )
        for controller in controllers:
            if controller is not None and controller.has_in_flight_update():
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
        self._request_dynamic_planet_update()

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

        # Lowest-priority idle work: keep satellite, aircraft, and cyclone overlays fresh.
        if (
            self._aircraft_layer_enabled()
            and self._aircraft_projection_next_refresh_delay_ms() == 0
        ):
            self.reproject_aircraft_overlay()
            return

        if (
            self._tropical_cyclone_layer_enabled()
            and self._tropical_cyclone_projection_next_refresh_delay_ms() == 0
        ):
            self.reproject_tropical_cyclone_overlay()
            return

        if (
            self._satellite_layer_enabled()
            and self._satellite_projection_next_refresh_delay_ms() == 0
        ):
            self.reproject_satellite_overlay()
            return

        if (
            self._cloud_layer_enabled()
            and self._cloud_projection_next_refresh_delay_ms() == 0
        ):
            self._start_cloud_projection_update(reason="scheduler")
            return

    def _status_line_message(self) -> str:
        simplified_view_mode = self._effective_simplified_view_mode()
        if simplified_view_mode == "labels":
            return "Simplified: with labels [Space]"
        if simplified_view_mode == "nolabels":
            return "Simplified: no labels [Space]"
        vertical_bar = "\u23ae"
        parts: list[str] = []
        cloud_message = self._cloud_status_line()
        if cloud_message:
            parts.append(cloud_message)
        cyclone_message = self._tropical_cyclone_status_line()
        if cyclone_message:
            parts.append(cyclone_message)
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
        if self._geo_satellite_mode_active():
            state = self.geosatellite_state
            if state.banner_text:
                detail = _strip_status_prefix(state.banner_text, "Geo-sat")
                detail_lower = detail.lower()
                if detail_lower.startswith("+"):
                    detail = detail[1:].strip()
                    detail_lower = detail.lower()
                if detail_lower.startswith("downloading"):
                    return _status_segment(_STATUS_CLOUD, "Geo-sat + Downloading")
                if detail_lower.startswith("projecting"):
                    return _status_segment(_STATUS_CLOUD, "Geo-sat + Projecting")
                if any(
                    token in detail_lower
                    for token in ("timed out", "error", "failed", "failure")
                ):
                    return _status_segment(_STATUS_CLOUD, "Geo-sat + error")
                return _status_segment(_STATUS_CLOUD, f"Geo-sat + {detail}")
            captured_at_utc = (
                state.captured_at_utc
                or (state.meta.time_utc if state.meta is not None else None)
                or state.fetched_at_utc
                or state.last_time_utc
            )
            if captured_at_utc is not None:
                try:
                    return _status_segment(
                        _STATUS_CLOUD,
                        f"Geo-sat + {captured_at_utc.astimezone(timezone.utc).strftime('%H:%MZ')}",
                    )
                except Exception:
                    pass
            return _status_segment(_STATUS_CLOUD, "Geo-sat + idle")
        sat = self.cloud_state.current_satellite or self._predicted_cloud_satellite()
        sat_label = _cloud_source_label(sat)
        if self.cloud_state.banner_text:
            detail = _strip_status_prefix(self.cloud_state.banner_text, "Clouds:")
            detail_lower = detail.lower()
            if detail_lower.startswith("downloading"):
                return _status_segment(_STATUS_CLOUD, f"{sat_label} downloading")
            if any(
                token in detail_lower
                for token in ("timed out", "error", "failed", "failure")
            ):
                return _status_segment(_STATUS_CLOUD, f"{sat_label} failed")
            return _status_segment(_STATUS_CLOUD, f"{sat_label} {detail}")
        meta = self.cloud_state.meta
        if meta is not None:
            try:
                t = meta.time_utc.strftime("%H:%MZ")
                sat_label = _cloud_source_label(meta.satellite or sat)
                source_ratio = self.cloud_state.source_completeness_ratio
                if source_ratio is None:
                    expected = self.cloud_state.source_expected_count
                    available = self.cloud_state.source_available_count
                    if (
                        expected is not None
                        and available is not None
                        and int(expected) > 0
                    ):
                        source_ratio = float(available) / float(expected)
                if source_ratio is not None and float(source_ratio) < 0.999:
                    return _status_segment(_STATUS_CLOUD, f"{sat_label} ? {t}")
                coverage = self.cloud_state.coverage_ratio
                if coverage is not None and coverage < 0.999:
                    pct = int(round(max(0.0, min(1.0, float(coverage))) * 100.0))
                    return _status_segment(_STATUS_CLOUD, f"{sat_label} {pct}% {t}")
                return _status_segment(_STATUS_CLOUD, f"{sat_label} {t}")
            except Exception:
                pass
        return _status_segment(_STATUS_CLOUD, f"{sat_label} idle")

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
            if detail.casefold().startswith("cache"):
                detail = "cache"
            return _status_segment(_STATUS_TERRAIN, detail)
        return ""

    def _water_overlay_status_line(self) -> str:
        if self.water_overlay_opacity <= 0.0:
            return _status_segment(_STATUS_WATER, "", hidden=True)
        state = self.water_overlay_state
        if state.banner_text:
            detail = _strip_status_prefix(state.banner_text, "Water:")
            if detail.casefold() in {"http 504", "offline", "network error"}:
                detail = "unavailable"
            return _status_segment(_STATUS_WATER, detail)
        sea_count = (
            "?" if state.sea_level_dots is None else str(len(state.sea_level_dots))
        )
        inland_count = "?" if state.inland_dots is None else str(len(state.inland_dots))
        detail = f"{sea_count}+{inland_count}"
        if str(state.current_source or "").strip().casefold().startswith("water: cache"):
            detail = "cache"
        return _status_segment(_STATUS_WATER, detail)

    def _urban_outline_status_line(self) -> str:
        if self.urban_outline_opacity <= 0.0:
            return _status_segment(_STATUS_URBAN, "", hidden=True)
        if self.urban_outline_state.banner_text:
            detail = _strip_status_prefix(
                self.urban_outline_state.banner_text,
                "Urban outline:",
            )
            return _status_segment(_STATUS_URBAN, detail)
        base_count = self.urban_outline_state.base_outline_count
        skyscraper_count = self.urban_outline_state.skyscraper_outline_count
        source_name = _urban_outline_source_name(
            self.urban_outline_state.current_source
        )

        def with_source(detail: str) -> str:
            return f"{source_name} {detail}" if source_name else detail

        if base_count is not None or skyscraper_count is not None:
            if base_count is not None and skyscraper_count is not None:
                return _status_segment(
                    _STATUS_URBAN, with_source(f"{base_count}+{skyscraper_count}")
                )
            if base_count is not None:
                return _status_segment(_STATUS_URBAN, with_source(str(base_count)))
            return _status_segment(_STATUS_URBAN, with_source(str(skyscraper_count)))
        if self.urban_outline_state.outlines is not None:
            return _status_segment(
                _STATUS_URBAN,
                with_source(str(len(self.urban_outline_state.outlines))),
            )
        if self.urban_outline_state.current_source:
            detail = _urban_outline_source_name(self.urban_outline_state.current_source)
            return _status_segment(_STATUS_URBAN, detail)
        return ""

    def _tropical_cyclone_status_line(self) -> str:
        if not self.show_tropical_cyclone_overlay:
            return _status_segment(_STATUS_TROPICAL_CYCLONE, "", hidden=True)
        state = self.tropical_cyclone_state
        if state is None:
            return _status_segment(_STATUS_TROPICAL_CYCLONE, "idle")
        if state.banner_text:
            detail = _strip_status_prefix(state.banner_text, "Typhoon:")
            return _status_segment(_STATUS_TROPICAL_CYCLONE, detail)
        snapshots = state.snapshots
        if not snapshots:
            if state.snapshot_collection is not None:
                return _status_segment(_STATUS_TROPICAL_CYCLONE, "none")
            return _status_segment(_STATUS_TROPICAL_CYCLONE, "idle")
        collection = state.snapshot_collection
        if hasattr(collection, "summary_text"):
            return _status_segment(_STATUS_TROPICAL_CYCLONE, collection.summary_text())
        if len(snapshots) == 1:
            snapshot = snapshots[0]
            advdate = snapshot.advdate_utc
            if advdate is not None:
                try:
                    advdate_text = advdate.astimezone(timezone.utc).strftime(
                        "%m-%d %H:%MZ"
                    )
                except Exception:
                    advdate_text = "?"
            else:
                advdate_text = "?"
            return _status_segment(
                _STATUS_TROPICAL_CYCLONE,
                f"{snapshot.storm_name} {advdate_text}",
            )
        preview_names = ", ".join(
            snapshot.storm_name for snapshot in snapshots[:3] if snapshot.storm_name
        )
        if len(snapshots) > 3:
            preview_names = (
                f"{preview_names}, +{len(snapshots) - 3}"
                if preview_names
                else f"+{len(snapshots) - 3}"
            )
        detail = f"{len(snapshots)} storms"
        if preview_names:
            detail = f"{detail}: {preview_names}"
        return _status_segment(_STATUS_TROPICAL_CYCLONE, detail)

    def _aircraft_status_line(self) -> str:
        if float(self.aircraft_opacity) <= 0.0:
            return _status_segment(_STATUS_AIRCRAFT, "", hidden=True)
        aircraft_state = self.aircraft_state
        if aircraft_state is None:
            return ""
        if aircraft_state.banner_text:
            detail = _strip_status_prefix(aircraft_state.banner_text, "Aircraft:")
            if detail.startswith("cache"):
                detail = "cache"
            return _status_segment(_STATUS_AIRCRAFT, detail)
        if aircraft_state.last_success_utc is None:
            return _status_segment(_STATUS_AIRCRAFT, "idle")
        return _status_segment(
            _STATUS_AIRCRAFT, aircraft_state.last_success_utc.strftime("%H:%MZ")
        )

    def _satellite_status_line(self) -> str:
        if float(self.satellite_opacity) <= 0.0:
            return _status_segment(_STATUS_SATELLITE, "", hidden=True)
        satellite_state = self.satellite_state
        if satellite_state is None:
            return ""
        if satellite_state.banner_text:
            detail = _strip_status_prefix(satellite_state.banner_text, "Satellites:")
            if detail.startswith("cache"):
                detail = "cache"
            return _status_segment(_STATUS_SATELLITE, detail)
        if satellite_state.element_epoch_utc is None:
            return _status_segment(_STATUS_SATELLITE, "idle")
        return _status_segment(
            _STATUS_SATELLITE, satellite_state.element_epoch_utc.strftime("%H:%MZ")
        )

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
        self.state.dynamic_planets = None
        self.state.dynamic_planet_bucket = None
        self.state.sky_disc_image = payload["sky_disc"]
        self.state.night_light_glow_profile = payload.get("night_light_glow_profile")

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
            if not _startup_night_light_requires_warmup(self, payload) or self.state.night_light_glow_profile is not None:
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

    def _request_dynamic_planet_update(self) -> None:
        if str(getattr(self, "presentation_id", "scenic")).strip().lower() != "scenic":
            return
        if self.state.celestial_data is None or self._viewport_interaction_active():
            return
        current_time = self._current_time_obj()
        bucket = int(float(current_time.unix) // 2.0)
        if self.state.dynamic_planet_bucket == bucket:
            return
        ephemeris = self._ephemeris
        if ephemeris is None:
            ephemeris = load_ephemeris()
        if self._sky_worker.update_planets(
            ephemeris=ephemeris,
            viewer_data=self._viewer_data_for_render(),
            time_obj=current_time,
        ):
            self.state.dynamic_planet_bucket = bucket

    def _on_planet_data_calculated(self, payload: dict) -> None:
        if self._is_shutting_down:
            return
        planets = payload.get("planets")
        if not isinstance(planets, list):
            return
        self.state.dynamic_planets = planets
        self.request_client_update()

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
            ephemeris = load_ephemeris()
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
            star_data_policy=getattr(
                self,
                "star_data_policy",
                "scenic_view_scoped",
            ),
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
            terrain_secondary_ridges_altaz_layers=self.state.terrain_secondary_ridges_altaz_layers,
            terrain_secondary_ridges_distances_m_layers=self.state.terrain_secondary_ridges_distances_m_layers,
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
        self._queue_aircraft_debug_snapshot(payload)

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
