from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
import re
from datetime import datetime, timezone
import logging

from ..satellites.cache import resolve_satellite_elements_for_time
from ..satellites.project import find_satellite_altaz
from ..satellites.types import SatelliteOmmRecord
from .models import SearchJumpTarget
from .query import parse_search_query

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SatelliteTargetSpec:
    label: str
    aliases: tuple[str, ...]


SATELLITE_TARGETS: tuple[SatelliteTargetSpec, ...] = (
    SatelliteTargetSpec(
        label="ISS",
        aliases=(
            "ISS",
            "ISS (ZARYA)",
            "ZARYA",
            "INTERNATIONAL SPACE STATION",
        ),
    ),
    SatelliteTargetSpec(
        label="JWST",
        aliases=(
            "JWST",
            "JAMES WEBB SPACE TELESCOPE",
            "JAMES WEBB",
        ),
    ),
    SatelliteTargetSpec(
        label="Voyager 1",
        aliases=(
            "VOYAGER 1",
            "VOYAGER-1",
        ),
    ),
    SatelliteTargetSpec(
        label="Voyager 2",
        aliases=(
            "VOYAGER 2",
            "VOYAGER-2",
        ),
    ),
    SatelliteTargetSpec(
        label="Parker",
        aliases=(
            "PARKER",
            "PARKER SOLAR PROBE",
            "SOLAR PROBE PLUS",
        ),
    ),
)


def fetch_current_satellite_records(
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    enabled_groups: Sequence[str] = ("iss", "horizons"),
    timeout_s: float | None = None,
    target_time_utc: datetime | None = None,
    force_refresh: bool = False,
) -> dict[str, list[SatelliteOmmRecord]]:
    now_utc = target_time_utc or datetime.now(timezone.utc)
    records_by_group: dict[str, list[SatelliteOmmRecord]] = {}
    for group_key in enabled_groups:
        try:
            fetched = resolve_satellite_elements_for_time(
                group_key,
                target_time_utc=now_utc,
                time_mode="present",
                timeout_s=timeout_s if timeout_s is not None else 60.0,
                validity_seconds=-1 if force_refresh else None,
                observer_lat=observer_lat,
                observer_lon=observer_lon,
                observer_height_m=observer_height_m,
            )
        except Exception as exc:
            logger.info("Satellite fetch unavailable for %s: %s", group_key, exc)
            continue
        records_by_group[str(group_key)] = list(fetched.records)
    return records_by_group


def search_satellite_targets(
    query: str,
    *,
    target_time_utc: datetime | None = None,
) -> list[SearchJumpTarget]:
    spec = parse_search_query(query)
    search_text = spec.value or spec.raw
    if not search_text:
        return []

    normalized_query = _normalize_satellite_name(search_text)
    if not normalized_query:
        return []

    matching_specs = [
        target_spec
        for target_spec in SATELLITE_TARGETS
        if normalized_query in {_normalize_satellite_name(alias) for alias in target_spec.aliases}
    ]
    if not matching_specs:
        return []

    search_time_utc = target_time_utc or datetime.now(timezone.utc)
    return [
        SearchJumpTarget(
            label=target_spec.label,
            kind="satellite",
            sort_key=(0.0, target_spec.label.casefold()),
            object_key=target_spec.label,
            subtitle="Artificial satellite",
            target_time_utc=search_time_utc,
        )
        for target_spec in matching_specs
    ]


def resolve_satellite_target_altaz(
    target: SearchJumpTarget,
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    target_time_utc: datetime | None = None,
    records_by_group: Mapping[str, Sequence[SatelliteOmmRecord]] | None = None,
    fetch_records_by_group: Callable[[], Mapping[str, Sequence[SatelliteOmmRecord]]] | None = None,
    enabled_groups: Sequence[str] = ("iss", "horizons"),
) -> tuple[float, float] | None:
    object_key = target.object_key or target.label
    if not object_key:
        return None
    records = dict(records_by_group or {})
    if not records:
        if fetch_records_by_group is not None:
            records = dict(fetch_records_by_group())
        else:
            records = dict(
                fetch_current_satellite_records(
                    observer_lat=observer_lat,
                    observer_lon=observer_lon,
                    observer_height_m=observer_height_m,
                    enabled_groups=enabled_groups,
                    target_time_utc=target_time_utc or target.target_time_utc or datetime.now(timezone.utc),
                    force_refresh=True,
                )
            )
    if not records:
        return None
    time_obj = _datetime_to_astropy_time(
        target_time_utc or target.target_time_utc or datetime.now(timezone.utc)
    )
    return find_satellite_altaz(
        records,
        object_key=object_key,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
        observer_height_m=observer_height_m,
        time_obj=time_obj,
    )


def _normalize_satellite_name(value: str) -> str:
    parts = [part for part in re.split(r"[^0-9a-z]+", str(value).casefold()) if part]
    return " ".join(parts)


def _datetime_to_astropy_time(value: datetime):
    import astropy.time

    dt = value if value.tzinfo is not None else value.replace(tzinfo=timezone.utc)
    return astropy.time.Time(dt)
