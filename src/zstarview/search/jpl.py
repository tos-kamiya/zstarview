from __future__ import annotations

from datetime import datetime, timezone
import logging
from typing import Mapping
from collections.abc import Sequence
from collections.abc import Callable

from astropy import units as u
from astropy.coordinates import AltAz, CartesianRepresentation, EarthLocation, GCRS, SkyCoord

from ..satellites import (
    fetch_horizons_lookup,
    fetch_horizons_observer_csv,
    fetch_horizons_vector_csv,
)
from .constants import SOLAR_SYSTEM_BODY_QUERIES
from .models import SearchJumpTarget
from .query import parse_search_query

logger = logging.getLogger(__name__)

JPL_BYPASS_QUERIES = SOLAR_SYSTEM_BODY_QUERIES
JPL_SEARCH_RESULT_LIMIT = 500


def build_horizons_command(target: Mapping[str, object], *, group: str) -> str:
    spkid = str(target.get("spkid", "")).strip()
    if spkid:
        if group == "sb":
            return f"DES={spkid};"
        return spkid
    pdes = str(target.get("pdes", "")).strip()
    if pdes:
        if group == "sb":
            return f"DES={pdes};"
        return pdes
    name = str(target.get("name", "")).strip()
    if name:
        if group == "sb":
            return name if name.endswith(";") else f"{name};"
        return name
    return ""


def extract_horizons_altaz(rows: list[list[str]]) -> tuple[float, float] | None:
    for row in rows:
        numeric_values: list[float] = []
        for value in row:
            try:
                numeric_values.append(float(str(value).strip()))
            except (TypeError, ValueError):
                continue
        if len(numeric_values) >= 2:
            # Horizons observer CSV reports azimuth first and elevation second.
            return numeric_values[1], numeric_values[0]
    return None


def extract_horizons_state_vector(
    rows: list[list[str]],
) -> tuple[tuple[float, float, float], tuple[float, float, float]] | None:
    for row in rows:
        numeric_values: list[float] = []
        for value in row:
            try:
                numeric_values.append(float(str(value).strip()))
            except (TypeError, ValueError):
                continue
        if len(numeric_values) >= 7 and _looks_like_julian_date(numeric_values[0]):
            numeric_values = numeric_values[1:]
        if len(numeric_values) >= 6:
            return (
                (numeric_values[0], numeric_values[1], numeric_values[2]),
                (numeric_values[3], numeric_values[4], numeric_values[5]),
            )
    return None


def _looks_like_julian_date(value: float) -> bool:
    return 2_000_000.0 <= float(value) <= 3_000_000.0


def search_jpl_targets(
    query: str,
    *,
    target_time_utc: datetime | None = None,
    groups: Sequence[str] = ("mb", "sb"),
    timeout_s: float | None = None,
    lookup_base_url: str | None = None,
    lookup_fetch: Callable[..., dict[str, object]] | None = None,
) -> list[SearchJumpTarget]:
    spec = parse_search_query(query)
    search_text = spec.value or spec.raw
    if not search_text or search_text.casefold() in JPL_BYPASS_QUERIES:
        return []
    fetch_kwargs: dict[str, object] = {}
    if timeout_s is not None:
        fetch_kwargs["timeout_s"] = float(timeout_s)
    if lookup_base_url is not None:
        fetch_kwargs["base_url"] = lookup_base_url
    lookup_impl = lookup_fetch or fetch_horizons_lookup
    resolved_time_utc = target_time_utc or datetime.now(timezone.utc)
    targets: list[SearchJumpTarget] = []
    for group in tuple(str(group).strip() for group in groups if str(group).strip()):
        group_label = "major body" if group == "mb" else "small body" if group == "sb" else group
        try:
            lookup_payload = lookup_impl(search_text, group=group, **fetch_kwargs)
        except Exception as exc:
            logger.info("JPL %s lookup failed for %s: %s", group_label, query, exc)
            continue
        result = lookup_payload.get("result")
        if not isinstance(result, list) or not result:
            continue
        for item in result[:JPL_SEARCH_RESULT_LIMIT]:
            if not isinstance(item, dict):
                continue
            name = str(item.get("name", "")).strip()
            if not name:
                continue
            command = build_horizons_command(item, group=group)
            if not command:
                continue
            type_name = str(item.get("type", "")).strip()
            pdes = str(item.get("pdes", "")).strip()
            subtitle_parts = [part for part in (group_label, type_name, pdes) if part]
            targets.append(
                SearchJumpTarget(
                    label=name,
                    kind="jpl_body",
                    sort_key=(0.0 if group == "mb" else 1.0, name.casefold()),
                    subtitle=" / ".join(subtitle_parts) if subtitle_parts else "JPL Body",
                    object_key=str(item.get("spkid", "")).strip(),
                    command=command,
                    target_time_utc=resolved_time_utc,
                    jpl_group=group,
                )
            )
    return targets


def resolve_jpl_target_altaz(
    target: SearchJumpTarget,
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    target_time_utc: datetime | None = None,
    observer_fetch: Callable[..., list[list[str]]] | None = None,
    timeout_s: float | None = None,
    horizons_base_url: str | None = None,
) -> tuple[float, float] | None:
    command = str(target.command).strip()
    if not command:
        return None
    effective_target_time_utc = target_time_utc or target.target_time_utc or datetime.now(timezone.utc)
    fetch_kwargs: dict[str, object] = {}
    if timeout_s is not None:
        fetch_kwargs["timeout_s"] = float(timeout_s)
    if horizons_base_url is not None:
        fetch_kwargs["base_url"] = horizons_base_url
    observer_impl = observer_fetch or fetch_horizons_observer_csv
    logger.info(
        "Resolving JPL target alt/az: label=%s group=%s command=%s target_time_utc=%s observer=(lat=%s lon=%s height_m=%s)",
        str(target.label).strip() or "<unnamed>",
        str(target.jpl_group).strip() or "<none>",
        command,
        effective_target_time_utc.astimezone(timezone.utc).isoformat(),
        observer_lat,
        observer_lon,
        observer_height_m,
    )
    rows = observer_impl(
        command,
        target_time_utc=effective_target_time_utc,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
        observer_height_m=observer_height_m,
        **fetch_kwargs,
    )
    altaz = extract_horizons_altaz(rows)
    if altaz is None:
        logger.info(
            "Resolved JPL target alt/az: label=%s command=%s result=<none>",
            str(target.label).strip() or "<unnamed>",
            command,
        )
        return None
    logger.info(
        "Resolved JPL target alt/az: label=%s command=%s alt=%.1f az=%.1f",
        str(target.label).strip() or "<unnamed>",
        command,
        float(altaz[0]),
        float(altaz[1]) % 360.0,
    )
    return altaz


def resolve_jpl_target_state_vector(
    target: SearchJumpTarget,
    *,
    target_time_utc: datetime | None = None,
    vector_fetch: Callable[..., list[list[str]]] | None = None,
    timeout_s: float | None = None,
    horizons_base_url: str | None = None,
) -> tuple[datetime, tuple[float, float, float], tuple[float, float, float]] | None:
    command = str(target.command).strip()
    if not command:
        return None
    effective_target_time_utc = target_time_utc or target.target_time_utc or datetime.now(timezone.utc)
    fetch_kwargs: dict[str, object] = {}
    if timeout_s is not None:
        fetch_kwargs["timeout_s"] = float(timeout_s)
    if horizons_base_url is not None:
        fetch_kwargs["base_url"] = horizons_base_url
    vector_impl = vector_fetch or fetch_horizons_vector_csv
    logger.info(
        "Resolving JPL target state vector: label=%s group=%s command=%s target_time_utc=%s",
        str(target.label).strip() or "<unnamed>",
        str(target.jpl_group).strip() or "<none>",
        command,
        effective_target_time_utc.astimezone(timezone.utc).isoformat(),
    )
    rows = vector_impl(
        command,
        target_time_utc=effective_target_time_utc,
        **fetch_kwargs,
    )
    state_vector = extract_horizons_state_vector(rows)
    if state_vector is None:
        logger.info(
            "Resolved JPL target state vector: label=%s command=%s result=<none>",
            str(target.label).strip() or "<unnamed>",
            command,
        )
        return None
    position_km, velocity_km_s = state_vector
    logger.info(
        "Resolved JPL target state vector: label=%s command=%s x=%.3f y=%.3f z=%.3f",
        str(target.label).strip() or "<unnamed>",
        command,
        float(position_km[0]),
        float(position_km[1]),
        float(position_km[2]),
    )
    return effective_target_time_utc.astimezone(timezone.utc), position_km, velocity_km_s


def project_jpl_target_altaz_from_state_vector(
    target: SearchJumpTarget,
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    time_obj,
) -> tuple[float, float] | None:
    epoch_utc = target.horizons_epoch_utc or target.target_time_utc
    position_km = target.horizons_position_km
    velocity_km_s = target.horizons_velocity_km_s
    if epoch_utc is None or position_km is None or velocity_km_s is None:
        return None
    try:
        current_utc = time_obj.to_datetime(timezone.utc)
    except Exception:
        current_utc = datetime.now(timezone.utc)
    delta_seconds = (current_utc - epoch_utc.astimezone(timezone.utc)).total_seconds()
    x_km = float(position_km[0]) + float(velocity_km_s[0]) * delta_seconds
    y_km = float(position_km[1]) + float(velocity_km_s[1]) * delta_seconds
    z_km = float(position_km[2]) + float(velocity_km_s[2]) * delta_seconds

    location = EarthLocation(
        lat=float(observer_lat) * u.deg,
        lon=float(observer_lon) * u.deg,
        height=float(observer_height_m) * u.m,
    )
    altaz_frame = AltAz(obstime=time_obj, location=location)
    coords = SkyCoord(
        CartesianRepresentation(
            x=x_km * u.km,
            y=y_km * u.km,
            z=z_km * u.km,
        ),
        frame=GCRS(obstime=time_obj),
    )
    altaz = coords.transform_to(altaz_frame)
    return float(altaz.alt.deg), float(altaz.az.deg) % 360.0
