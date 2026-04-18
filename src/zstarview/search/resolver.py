from __future__ import annotations

from collections.abc import Callable, Sequence
from datetime import datetime, timezone

from ..location_resolver import project_place_target_to_altaz
from ..astro import radec_to_altaz
from .models import SearchJumpTarget, SearchResolution
from .query import parse_search_query, search_target_matches_query


def resolve_search_targets(
    query: str,
    local_targets: Sequence[SearchJumpTarget],
    *,
    satellite_search_callback: Callable[[str], Sequence[SearchJumpTarget]] | None = None,
    jpl_search_callback: Callable[[str], Sequence[SearchJumpTarget]] | None = None,
) -> SearchResolution:
    spec = parse_search_query(query)
    if not spec.normalized:
        return SearchResolution(query=spec.raw, candidates=tuple(), status="")

    local_matches = tuple(
        target for target in local_targets if search_target_matches_query(target, spec)
    )
    if local_matches:
        return SearchResolution(
            query=spec.raw,
            candidates=local_matches,
            selected_target=local_matches[0] if len(local_matches) == 1 else None,
            status=f"Found {len(local_matches)} local result(s)",
        )

    if satellite_search_callback is not None:
        satellite_query = spec.value or spec.raw
        satellite_candidates = tuple(satellite_search_callback(satellite_query))
        if satellite_candidates:
            return SearchResolution(
                query=spec.raw,
                candidates=satellite_candidates,
                selected_target=(
                    satellite_candidates[0] if len(satellite_candidates) == 1 else None
                ),
                status=(
                    f"Found {len(satellite_candidates)} satellite result(s)"
                ),
            )

    if jpl_search_callback is None:
        return SearchResolution(query=spec.raw, candidates=tuple(), status="")

    jpl_query = spec.value or spec.raw
    candidates = tuple(jpl_search_callback(jpl_query))
    return SearchResolution(
        query=spec.raw,
        candidates=candidates,
        selected_target=candidates[0] if len(candidates) == 1 else None,
        status=(
            f"Found {len(candidates)} JPL result(s)"
            if candidates
            else "No JPL results found"
        ),
    )


def compute_search_target_altaz(
    target: SearchJumpTarget,
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    target_time_utc: datetime | None = None,
    satellite_altaz_resolver: Callable[[SearchJumpTarget], tuple[float, float] | None] | None = None,
) -> tuple[float, float] | None:
    if target.kind in {"jpl_small_body", "jpl_body"}:
        if target.alt_deg is None or target.az_deg is None:
            return None
        return float(target.alt_deg), float(target.az_deg) % 360.0
    if target.kind == "satellite":
        if target.alt_deg is not None and target.az_deg is not None:
            return float(target.alt_deg), float(target.az_deg) % 360.0
        if satellite_altaz_resolver is not None:
            altaz = satellite_altaz_resolver(target)
            if altaz is not None:
                return float(altaz[0]), float(altaz[1]) % 360.0
        return None
    if target.kind == "place":
        if target.latitude_deg is None or target.longitude_deg is None:
            return None
        projection = project_place_target_to_altaz(
            observer_latitude_deg=float(observer_lat),
            observer_longitude_deg=float(observer_lon),
            observer_height_m=float(observer_height_m),
            target_latitude_deg=float(target.latitude_deg),
            target_longitude_deg=float(target.longitude_deg),
        )
        return float(projection.alt_deg), float(projection.az_deg) % 360.0
    if target_time_utc is None:
        target_time_utc = datetime.now(timezone.utc)
    alt, az = radec_to_altaz(
        float(target.ra_hours),
        float(target.dec_deg),
        float(observer_lat),
        float(observer_lon),
        float(observer_height_m),
        target_time_utc,
    )
    return float(alt), float(az) % 360.0
