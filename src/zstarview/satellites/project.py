from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone

import astropy.time
import skyfield.api
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    AltAz,
    CartesianRepresentation,
    EarthLocation,
    SkyCoord,
)

from ..satellite_constants import SATELLITE_HORIZONS_CACHE_KEY, SATELLITE_ISS_CACHE_KEY
from .fetch import build_earth_satellites
from .types import (
    SatelliteOmmRecord,
    SatelliteOverlayPoint,
    satellite_altaz_from_record,
)

_DEFAULT_GROUP_ORDER = (SATELLITE_ISS_CACHE_KEY, SATELLITE_HORIZONS_CACHE_KEY)
_MAX_MARKERS_BY_GROUP = {
    SATELLITE_ISS_CACHE_KEY: 1,
    SATELLITE_HORIZONS_CACHE_KEY: 4,
}
_MARKER_SCALE_BY_GROUP = {
    SATELLITE_ISS_CACHE_KEY: 0.42,
    SATELLITE_HORIZONS_CACHE_KEY: 0.30,
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
        group_points = [
            point for point in projected_points if point.group_key == str(group_key)
        ]
        group_points.sort(key=lambda point: float(point.alt_deg), reverse=True)
        overlay_points.extend(
            group_points[: _MAX_MARKERS_BY_GROUP.get(group_key, len(group_points))]
        )
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
    points: list[SatelliteOverlayPoint] = []
    for group_key in _iter_group_order(records_by_group):
        records = records_by_group.get(group_key) or ()
        if group_key == SATELLITE_HORIZONS_CACHE_KEY:
            for record in records:
                alt_deg, az_deg = _project_horizons_record_to_altaz(
                    record,
                    observer_lat=observer_lat,
                    observer_lon=observer_lon,
                    observer_height_m=observer_height_m,
                    time_obj=time_obj,
                )
                if alt_deg is None or az_deg is None:
                    continue
                points.append(
                    SatelliteOverlayPoint(
                        group_key=str(group_key),
                        satellite_name=str(
                            record.get("OBJECT_NAME")
                            or record.get("HORIZONS_TARGET_NAME")
                            or group_key.upper()
                        ),
                        alt_deg=float(alt_deg),
                        az_deg=float(az_deg),
                        marker_scale=float(_MARKER_SCALE_BY_GROUP.get(group_key, 0.13)),
                    )
                )
            continue
        ts = skyfield.api.load.timescale()
        t = ts.from_astropy(time_obj)
        observer = skyfield.api.Topos(
            latitude_degrees=float(observer_lat),
            longitude_degrees=float(observer_lon),
            elevation_m=float(observer_height_m),
        )
        satellites = build_earth_satellites(records, ts=ts)
        for satellite in satellites:
            topocentric = (satellite - observer).at(t)
            alt, az, _distance = topocentric.altaz()
            points.append(
                SatelliteOverlayPoint(
                    group_key=str(group_key),
                    satellite_name=str(
                        getattr(satellite, "name", "") or group_key.upper()
                    ),
                    alt_deg=float(alt.degrees),
                    az_deg=float(az.degrees),
                    marker_scale=float(_MARKER_SCALE_BY_GROUP.get(group_key, 0.13)),
                )
            )
    return points


def _project_horizons_record_to_altaz(
    record: SatelliteOmmRecord,
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    time_obj: astropy.time.Time,
) -> tuple[float | None, float | None]:
    state_vector = _horizons_state_vector_from_record(record)
    if state_vector is None:
        return satellite_altaz_from_record(record)
    epoch_utc = _horizons_epoch_utc_from_record(record)
    if epoch_utc is None:
        return satellite_altaz_from_record(record)
    position_km, velocity_km_s = state_vector
    try:
        current_utc = time_obj.to_datetime(timezone.utc)
    except Exception:
        current_utc = datetime.now(timezone.utc)
    delta_seconds = (current_utc - epoch_utc).total_seconds()
    x_km = position_km[0] + velocity_km_s[0] * delta_seconds
    y_km = position_km[1] + velocity_km_s[1] * delta_seconds
    z_km = position_km[2] + velocity_km_s[2] * delta_seconds
    location = EarthLocation(
        lat=float(observer_lat) * u.deg,
        lon=float(observer_lon) * u.deg,
        height=float(observer_height_m) * u.m,
    )
    altaz_frame = AltAz(obstime=time_obj, location=location)
    coords = SkyCoord(
        CartesianRepresentation(
            x=float(x_km) * u.km,
            y=float(y_km) * u.km,
            z=float(z_km) * u.km,
        ),
        frame=GCRS(obstime=time_obj),
    )
    altaz = coords.transform_to(altaz_frame)
    return float(altaz.alt.deg), float(altaz.az.deg) % 360.0


def _horizons_state_vector_from_record(
    record: SatelliteOmmRecord,
) -> tuple[tuple[float, float, float], tuple[float, float, float]] | None:
    position_keys = (
        "HORIZONS_X_KM",
        "HORIZONS_Y_KM",
        "HORIZONS_Z_KM",
    )
    velocity_keys = (
        "HORIZONS_VX_KM_S",
        "HORIZONS_VY_KM_S",
        "HORIZONS_VZ_KM_S",
    )
    try:
        position = tuple(float(record[key]) for key in position_keys)
        velocity = tuple(float(record[key]) for key in velocity_keys)
    except (KeyError, TypeError, ValueError):
        return None
    return position, velocity


def _horizons_epoch_utc_from_record(record: SatelliteOmmRecord) -> datetime | None:
    raw_epoch = str(record.get("EPOCH", "")).strip()
    if not raw_epoch:
        return None
    try:
        parsed = datetime.fromisoformat(raw_epoch.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        return parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _iter_group_order(
    records_by_group: Mapping[str, Sequence[SatelliteOmmRecord]],
) -> list[str]:
    ordered = [
        group_key for group_key in _DEFAULT_GROUP_ORDER if group_key in records_by_group
    ]
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
