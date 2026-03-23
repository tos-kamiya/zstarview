from __future__ import annotations

from collections.abc import Mapping, Sequence
import re

import astropy.time
from skyfield.api import Topos
import skyfield.api

from .fetch import build_earth_satellites
from .types import SatelliteOmmRecord, SatelliteOverlayPoint

_DEFAULT_GROUP_ORDER = ("station", "starlink")
_MAX_MARKERS_BY_GROUP = {
    "station": 8,
    "starlink": 20,
}
_MARKER_SCALE_BY_GROUP = {
    "station": 0.3,
    "starlink": 0.156,
}


def project_satellite_records(
    records_by_group: Mapping[str, Sequence[SatelliteOmmRecord]],
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    time_obj: astropy.time.Time,
) -> list[SatelliteOverlayPoint]:
    projected_points = compute_satellite_altaz_points(
        records_by_group,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
        observer_height_m=observer_height_m,
        time_obj=time_obj,
    )

    overlay_points: list[SatelliteOverlayPoint] = []
    for group_key in _iter_group_order(records_by_group):
        group_points = [point for point in projected_points if point.group_key == str(group_key)]
        group_points.sort(key=lambda point: float(point.alt_deg), reverse=True)
        overlay_points.extend(group_points[: _MAX_MARKERS_BY_GROUP.get(group_key, len(group_points))])
    return overlay_points


def find_satellite_altaz(
    records_by_group: Mapping[str, Sequence[SatelliteOmmRecord]],
    *,
    object_key: str,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    time_obj: astropy.time.Time,
) -> tuple[float, float] | None:
    query = str(object_key).strip().casefold()
    if not query:
        return None
    best_score = 0
    best_altaz: tuple[float, float] | None = None
    for point in compute_satellite_altaz_points(
        records_by_group,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
        observer_height_m=observer_height_m,
        time_obj=time_obj,
    ):
        score = _satellite_match_score(query, str(point.satellite_name))
        if score > best_score:
            best_score = score
            best_altaz = (float(point.alt_deg), float(point.az_deg))
    return best_altaz


def compute_satellite_altaz_points(
    records_by_group: Mapping[str, Sequence[SatelliteOmmRecord]],
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    time_obj: astropy.time.Time,
) -> list[SatelliteOverlayPoint]:
    ts = skyfield.api.load.timescale()
    t = ts.from_astropy(time_obj)
    observer = Topos(
        latitude_degrees=float(observer_lat),
        longitude_degrees=float(observer_lon),
        elevation_m=float(observer_height_m),
    )

    points: list[SatelliteOverlayPoint] = []
    for group_key in _iter_group_order(records_by_group):
        records = records_by_group.get(group_key) or ()
        satellites = build_earth_satellites(records, ts=ts)
        for satellite in satellites:
            topocentric = (satellite - observer).at(t)
            alt, az, _distance = topocentric.altaz()
            points.append(
                SatelliteOverlayPoint(
                    group_key=str(group_key),
                    satellite_name=str(getattr(satellite, "name", "") or group_key.upper()),
                    alt_deg=float(alt.degrees),
                    az_deg=float(az.degrees),
                    marker_scale=float(_MARKER_SCALE_BY_GROUP.get(group_key, 0.13)),
                    show_label=group_key == "station",
                )
            )
    return points


def _iter_group_order(records_by_group: Mapping[str, Sequence[SatelliteOmmRecord]]) -> list[str]:
    ordered = [group_key for group_key in _DEFAULT_GROUP_ORDER if group_key in records_by_group]
    for group_key in records_by_group:
        if group_key not in ordered:
            ordered.append(str(group_key))
    return ordered


def _satellite_match_score(query: str, satellite_name: str) -> int:
    aliases = _satellite_name_aliases(satellite_name)
    if query in aliases:
        return 3
    if any(alias.startswith(f"{query} ") for alias in aliases):
        return 2
    normalized_name = _normalize_satellite_name(satellite_name)
    if query and query in normalized_name:
        return 1
    return 0


def _satellite_name_aliases(satellite_name: str) -> set[str]:
    normalized_name = _normalize_satellite_name(satellite_name)
    if not normalized_name:
        return set()
    aliases = {normalized_name}
    first_token = normalized_name.split(" ", 1)[0]
    if first_token:
        aliases.add(first_token)
    return aliases


def _normalize_satellite_name(value: str) -> str:
    parts = [part for part in re.split(r"[^0-9a-z]+", str(value).casefold()) if part]
    return " ".join(parts)
