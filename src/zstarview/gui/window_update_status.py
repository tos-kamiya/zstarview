from __future__ import annotations

import math
from datetime import datetime, timezone

from ..clouddisc.providers.select import GOES_SATELLITES

_STATUS_CLOUD = "☁"
_STATUS_WATER = "W"
_STATUS_SATELLITE = "🛰"
_STATUS_AIRCRAFT = "✈"
_STATUS_METEOR = "M"
_STATUS_TROPICAL_CYCLONE = "TC"
_STATUS_TERRAIN = "△"
_STATUS_URBAN = "🂓"
_STATUS_ROAD = "R"
_STATUS_PRECIPITATION = "☂"


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


class SkyWindowStatusUpdatesMixin:
    def _status_line_message(self) -> str:
        vertical_bar = "\u23ae"
        dynamic_parts: list[str] = []
        fixed_parts: list[str] = []
        cloud_message = self._cloud_status_line()
        if cloud_message:
            dynamic_parts.append(cloud_message)
        precipitation_status_line = getattr(self, "_precipitation_status_line", None)
        if callable(precipitation_status_line):
            precipitation_message = precipitation_status_line()
            if precipitation_message:
                dynamic_parts.append(precipitation_message)
        cyclone_message = self._tropical_cyclone_status_line()
        if cyclone_message:
            dynamic_parts.append(cyclone_message)
        satellite_status_line = self._satellite_status_line
        if callable(satellite_status_line):
            satellite_message = satellite_status_line()
            if satellite_message:
                dynamic_parts.append(satellite_message)
        aircraft_status_line = self._aircraft_status_line
        if callable(aircraft_status_line):
            aircraft_message = aircraft_status_line()
            if aircraft_message:
                dynamic_parts.append(aircraft_message)
        meteor_status_line = getattr(self, "_meteor_status_line", None)
        meteor_message = meteor_status_line() if callable(meteor_status_line) else ""
        if meteor_message:
            dynamic_parts.append(meteor_message)
        jpl_message = self._jpl_small_body_status_line()
        if jpl_message:
            dynamic_parts.append(jpl_message)
        terrain_message = self._terrain_horizon_status_line()
        if terrain_message:
            fixed_parts.append(terrain_message)
        water_message = self._water_overlay_status_line()
        if water_message:
            fixed_parts.append(water_message)
        road_status_line = getattr(self, "_road_night_lights_status_line", None)
        if callable(road_status_line):
            road_message = road_status_line()
            if road_message:
                fixed_parts.append(road_message)
        urban_message = self._urban_outline_status_line()
        if urban_message:
            fixed_parts.append(urban_message)

        def format_line(parts: list[str]) -> str:
            return (
                f"{vertical_bar} %s {vertical_bar}"
                % f" {vertical_bar} ".join(parts)
                if parts
                else ""
            )

        return "\n".join(
            line for line in (format_line(fixed_parts), format_line(dynamic_parts)) if line
        )

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
        if (
            str(state.current_source or "")
            .strip()
            .casefold()
            .startswith("water: cache")
        ):
            detail = "cache"
        return _status_segment(_STATUS_WATER, detail)

    def _road_night_lights_status_line(self) -> str:
        if float(getattr(self, "road_night_lights_opacity", 0.0)) <= 0.0:
            return _status_segment(_STATUS_ROAD, "", hidden=True)
        status = str(getattr(self, "road_night_lights_status", "") or "").strip()
        return _status_segment(_STATUS_ROAD, status) if status else ""

    def _precipitation_status_line(self) -> str:
        if float(getattr(self, "precipitation_opacity", 0.0)) <= 0.0:
            return ""
        status = str(getattr(self, "precipitation_status", "") or "").strip()
        if status != "ready":
            return _status_segment(
                _STATUS_PRECIPITATION, f"Open Meteo {status or 'loading'}"
            )
        forecast_time = getattr(self, "precipitation_forecast_time_utc", None)
        if isinstance(forecast_time, datetime):
            return _status_segment(
                _STATUS_PRECIPITATION, f"Open Meteo {forecast_time:%H:%MZ}"
            )
        return _status_segment(_STATUS_PRECIPITATION, "Open Meteo")

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

    def _meteor_status_line(self) -> str:
        if float(getattr(self, "meteor_opacity", 0.0)) <= 0.0:
            return _status_segment(_STATUS_METEOR, "", hidden=True)
        state = getattr(self, "meteor_state", None)
        if state is None:
            return _status_segment(_STATUS_METEOR, "idle")
        if state.banner_text:
            return _status_segment(
                _STATUS_METEOR,
                _strip_status_prefix(state.banner_text, "GMN meteors:"),
            )
        if state.result is None:
            return _status_segment(_STATUS_METEOR, "idle")
        oldest_hours = max(
            0,
            math.ceil(
                (state.result.display_time_utc - state.result.window_start_utc)
                .total_seconds()
                / 3600.0
            ),
        )
        newest_hours = max(
            0,
            math.floor(
                (state.result.display_time_utc - state.result.window_end_utc)
                .total_seconds()
                / 3600.0
            ),
        )
        detail = (
            f"{len(state.result.trails)}, "
            f"{oldest_hours}-{newest_hours}h ago"
        )
        if state.result.used_stale_index or state.result.used_stale_files:
            detail += " cache"
        return _status_segment(_STATUS_METEOR, detail)

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
