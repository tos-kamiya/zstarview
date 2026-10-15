"""Cache coverage reports without rendering or live acquisition."""

from __future__ import annotations

from datetime import datetime, timedelta
from typing import Any

from ..satellite_constants import SATELLITE_GROUP_VALIDITY_SECONDS
from .inventory import Inventory
from .models import MINUTE, UTC, HistoryRecord, ceil_minute, iso
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


def _satellite_object_reports(
    values: list[HistoryRecord],
    bounds: tuple[datetime, datetime] | None,
) -> list[dict[str, Any]]:
    snapshots_by_time: dict[datetime, HistoryRecord] = {}
    for record in sorted(values, key=lambda item: (item.fetched_at, str(item.path))):
        snapshots_by_time[record.fetched_at] = record
    snapshots = sorted(snapshots_by_time.items())
    accumulated: dict[str, dict[str, Any]] = {}
    for index, (fetched_at, snapshot) in enumerate(snapshots):
        tail_seconds = SATELLITE_GROUP_VALIDITY_SECONDS[
            snapshot.layer.removeprefix("satellite-")
        ]
        if index + 1 < len(snapshots):
            raw_end = snapshots[index + 1][0]
            end = ceil_minute(raw_end)
        else:
            raw_end = fetched_at + timedelta(seconds=tail_seconds)
            end = ceil_minute(raw_end)
            if end == raw_end:
                end += MINUTE
        start = ceil_minute(fetched_at)
        if end <= start:
            continue
        for observation in snapshot.satellite_objects:
            item = accumulated.setdefault(
                observation.object_id,
                {
                    "object_id": observation.object_id,
                    "name": observation.name,
                    "sources": set(),
                    "source_times": [],
                    "fetched_times": [],
                    "ranges": [],
                    "acquisitions": set(),
                },
            )
            item["sources"].add(observation.source)
            item["source_times"].append(observation.source_time)
            item["fetched_times"].append(observation.fetched_at)
            item["ranges"].append((start, end))
            item["acquisitions"].add(
                (observation.source, observation.source_time, observation.fetched_at)
            )

    reports = []
    for item in sorted(
        accumulated.values(), key=lambda value: (value["name"].casefold(), value["object_id"])
    ):
        intervals = sorted(item["ranges"])
        merged: list[tuple[datetime, datetime]] = []
        for begin, end in intervals:
            if bounds is not None:
                begin = max(begin, bounds[0])
                end = min(end, bounds[1])
            if begin >= end:
                continue
            if merged and begin <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((begin, end))
        source_times = item["source_times"]
        fetched_times = item["fetched_times"]
        acquisitions = sorted(item["acquisitions"], key=lambda value: value[2])
        reports.append(
            {
                "object_id": item["object_id"],
                "name": item["name"],
                "sources": sorted(item["sources"]),
                "source_time_range": {
                    "min": iso(min(source_times)),
                    "max": iso(max(source_times)),
                },
                "fetched_at_range": {
                    "min": iso(min(fetched_times)),
                    "max": iso(max(fetched_times)),
                },
                "ranges": serialize_ranges(merged),
                "acquisitions": [
                    {
                        "source": source,
                        "source_time": iso(source_time),
                        "fetched_at": iso(fetched_at),
                    }
                    for source, source_time, fetched_at in acquisitions
                ],
            }
        )
    return reports


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
        if name in {"satellite-iss", "satellite-horizons"}:
            layers[name]["objects"] = _satellite_object_reports(values, bounds)
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
        "satellite_selection_policy": {
            "basis": "acquisition-time",
            "between_records": "hold-until-next-successful-acquisition",
            "final_tail_seconds_by_group": SATELLITE_GROUP_VALIDITY_SECONDS,
            "source_epoch_controls_expiry": False,
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
        "Satellites: hold from acquisition until the next acquisition; final tail: 24 hours; source epoch does not set expiry.",
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
        for item in layer.get("objects", []):
            lines.append(
                f"  object: {item['name']} id={item['object_id']} "
                f"sources={','.join(item['sources']) or 'unknown'}"
            )
            for field in ("source_time_range", "fetched_at_range"):
                value = item[field]
                lines.append(f"    {field}: {value['min']} .. {value['max']}")
            lines.append(
                "    ranges: "
                + (
                    ", ".join(
                        f"[{r['start']}, {r['end']})" for r in item["ranges"]
                    )
                    or "none"
                )
            )
            for acquisition in item["acquisitions"]:
                lines.append(
                    "    acquisition: "
                    f"source={acquisition['source']} "
                    f"source_time={acquisition['source_time']} "
                    f"fetched_at={acquisition['fetched_at']}"
                )
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
