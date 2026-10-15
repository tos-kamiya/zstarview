"""Cache coverage reports without rendering or live acquisition."""

from __future__ import annotations

from datetime import datetime
from typing import Any

from .inventory import Inventory
from .models import MINUTE, UTC, ceil_minute, iso
from .selection import HistoryIndex, coverage, intersect


def serialize_ranges(ranges: list[tuple[datetime, datetime]]) -> list[dict[str, str]]:
    return [{"start": iso(start), "end": iso(end)} for start, end in ranges]


def gaps(
    ranges: list[tuple[datetime, datetime]], start: datetime, end: datetime
) -> list[tuple[datetime, datetime]]:
    result = []
    cursor = start
    for begin, finish in ranges:
        if cursor < begin:
            result.append((cursor, begin))
        cursor = finish
    if cursor < end:
        result.append((cursor, end))
    return result


def build_report(
    inventory: Inventory,
    enabled: set[str],
    *,
    start: datetime | None = None,
    end: datetime | None = None,
    ephemeris_ranges: list[tuple[datetime, datetime]] | None = None,
    dependencies: list[dict[str, str]] | None = None,
) -> dict[str, Any]:
    records = [r for values in inventory.records.values() for r in values]
    first = start or (min(r.starts_at for r in records) if records else None)
    last = end or (max(r.expires_at for r in records) + MINUTE if records else None)
    bounds = (
        (ceil_minute(first), ceil_minute(last))
        if first and last and first < last
        else None
    )
    layers = {}
    common: list[tuple[datetime, datetime]] | None = None
    for name, values in inventory.records.items():
        ranges = coverage(HistoryIndex(values), *bounds) if bounds else []
        if name in enabled:
            common = ranges if common is None else intersect(common, ranges)
        source_times = [r.source_time for r in values]
        fetched_times = [r.fetched_at for r in values]
        layers[name] = {
            "enabled": name in enabled,
            "record_count": len(values),
            "source_time_range": {
                "min": iso(min(source_times)),
                "max": iso(max(source_times)),
            }
            if values
            else None,
            "fetched_at_range": {
                "min": iso(min(fetched_times)),
                "max": iso(max(fetched_times)),
            }
            if values
            else None,
            "ranges": serialize_ranges(ranges),
            "gaps": serialize_ranges(gaps(ranges, *bounds)) if bounds else [],
            "diagnostics": dict(inventory.diagnostics[name]),
        }
    dynamic_common = common or []
    celestial_common = (
        intersect(dynamic_common, ephemeris_ranges)
        if ephemeris_ranges is not None
        else None
    )
    dependencies = dependencies or []
    return {
        "schema_version": 1,
        "verification_level": "metadata-only",
        "payloads_verified": False,
        "cloud_selection_policy": {
            "basis": "acquisition-time",
            "between_records": "hold-until-next-successful-acquisition",
            "native_tail_seconds": 1800,
            "geo_tail_seconds": 3600,
        },
        "checked_at_utc": iso(datetime.now(UTC)),
        "cache_root": str(inventory.root),
        "scan_range": {"start": iso(bounds[0]), "end": iso(bounds[1])}
        if bounds
        else None,
        "layers": layers,
        "common_status": "not-applicable"
        if not enabled
        else "available"
        if dynamic_common
        else "no-overlap"
        if all(inventory.records[n] for n in enabled)
        else "missing-history",
        "common_ranges": serialize_ranges(dynamic_common),
        "ephemeris_ranges": serialize_ranges(ephemeris_ranges)
        if ephemeris_ranges is not None
        else None,
        "celestial_common_ranges": serialize_ranges(celestial_common)
        if celestial_common is not None
        else None,
        "dependencies": dependencies,
        "complete_output_ready": False,
        "image_export_implemented": False,
    }


def format_report(report: dict[str, Any]) -> str:
    location = report["location"]
    lines = [
        f"Location: {location['display_name']}",
        f"Coordinates: {location['lat']}, {location['lon']}",
        f"Cache root: {report['cache_root']}",
        "Metadata only: payload integrity and rendering readiness are not verified.",
        "Clouds: hold from acquisition until the next acquisition; final tail: native 30 min, geo 60 min.",
    ]
    for name, layer in report["layers"].items():
        lines.append(
            f"{name}: records={layer['record_count']} enabled={str(layer['enabled']).lower()}"
        )
        for field in ("source_time_range", "fetched_at_range"):
            value = layer[field]
            lines.append(
                f"  {field}: {value['min']} .. {value['max']}"
                if value
                else f"  {field}: none"
            )
        for kind in ("ranges", "gaps"):
            lines.append(
                f"  {kind}: "
                + (
                    ", ".join(f"[{r['start']}, {r['end']})" for r in layer[kind])
                    or "none"
                )
            )
        if layer["diagnostics"]:
            lines.append(f"  diagnostics: {layer['diagnostics']}")
    lines.append(f"Common status: {report['common_status']}")
    for key in ("common_ranges", "celestial_common_ranges"):
        if report[key] is None:
            lines.append(f"{key}: unverified")
            continue
        lines.append(
            key
            + ": "
            + (", ".join(f"[{r['start']}, {r['end']})" for r in report[key]) or "none")
        )
    for item in report["dependencies"]:
        lines.append(f"Dependency {item['name']}: {item['status']} ({item['reason']})")
    return "\n".join(lines).encode("ascii", "backslashreplace").decode("ascii")
