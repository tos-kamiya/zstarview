from __future__ import annotations

from datetime import datetime
import logging
from typing import Mapping
from collections.abc import Sequence
from collections.abc import Callable

from ..satellites import fetch_horizons_lookup, fetch_horizons_observer_csv
from .models import SearchJumpTarget
from .query import parse_search_query

logger = logging.getLogger(__name__)

JPL_BYPASS_QUERIES = {"sun", "moon"}
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
            return numeric_values[-1], numeric_values[-2]
    return None


def search_jpl_targets(
    query: str,
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    target_time_utc: datetime,
    groups: Sequence[str] = ("mb", "sb"),
    timeout_s: float | None = None,
    lookup_base_url: str | None = None,
    horizons_base_url: str | None = None,
    lookup_fetch: Callable[..., dict[str, object]] | None = None,
    observer_fetch: Callable[..., list[list[str]]] | None = None,
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
    observer_fetch_kwargs: dict[str, object] = {}
    if timeout_s is not None:
        observer_fetch_kwargs["timeout_s"] = float(timeout_s)
    if horizons_base_url is not None:
        observer_fetch_kwargs["base_url"] = horizons_base_url
    lookup_impl = lookup_fetch or fetch_horizons_lookup
    observer_impl = observer_fetch or fetch_horizons_observer_csv
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
            try:
                rows = observer_impl(
                    command,
                    target_time_utc=target_time_utc,
                    observer_lat=observer_lat,
                    observer_lon=observer_lon,
                    observer_height_m=observer_height_m,
                    **observer_fetch_kwargs,
                )
            except Exception as exc:
                logger.info(
                    "JPL %s observer fetch failed for %s: %s",
                    group_label,
                    name,
                    exc,
                )
                continue
            alt_az = extract_horizons_altaz(rows)
            if alt_az is None:
                continue
            alt_deg, az_deg = alt_az
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
                    alt_deg=float(alt_deg),
                    az_deg=float(az_deg) % 360.0,
                    target_time_utc=target_time_utc,
                    jpl_group=group,
                )
            )
    targets.sort(key=lambda target: target.sort_key)
    return targets
