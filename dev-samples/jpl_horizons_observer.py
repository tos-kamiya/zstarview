#!/usr/bin/env python3
"""
jpl_horizons_observer.py
Version: 0.1.0

A small standalone command-line tool that resolves a JPL Horizons object
and fetches one observer ephemeris sample for a specific location and time.

This script is the second step in the small-body prototype sequence:

1. search the SBDB
2. resolve the Horizons target
3. fetch a single observer alt/az sample

It does not yet build a time series or a trajectory polyline.

------------------------------------------------------------
OVERVIEW
------------------------------------------------------------

Given a search string, the script:

1. calls the Horizons lookup API to resolve the target
2. calls the Horizons observer API for the specified observer position
3. prints the resulting altitude/azimuth sample

The lookup step is used to find the correct Horizons command/SPK ID before
querying the observer ephemeris.

------------------------------------------------------------
API USED
------------------------------------------------------------

Horizons Lookup API:
    https://ssd.jpl.nasa.gov/api/horizons_lookup.api

Horizons API:
    https://ssd.jpl.nasa.gov/api/horizons.api

------------------------------------------------------------
OPTIONS
------------------------------------------------------------

--group GROUP
    Restrict the lookup to a Horizons object group.
    Default: sb

--lat LATITUDE
    Observer latitude in degrees.

--lon LONGITUDE
    Observer longitude in degrees.

--height-m HEIGHT
    Observer height above sea level in meters.

--time UTC_ISO
    Observation time in UTC ISO-8601 format.
    Default: current UTC time rounded to the nearest minute.

--json
    Emit a structured JSON payload.

--compact
    Reduce output formatting to terse one-line summaries.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.satellites.fetch import (  # noqa: E402
    fetch_horizons_lookup,
    fetch_horizons_observer_csv,
)


VERSION = "0.1.0"
ALLOWED_GROUPS = {"ast", "com", "pln", "sat", "sct", "mb", "sb"}
DEFAULT_GROUP = "sb"


def _normalize_text(value: str) -> str:
    return " ".join(part for part in value.casefold().replace("(", " ").replace(")", " ").split() if part)


def _parse_utc_time(value: str | None) -> datetime:
    if not value:
        now = datetime.now(timezone.utc)
        return now.replace(second=0, microsecond=0)
    text = value.strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    parsed = datetime.fromisoformat(text)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _pick_exact_match(query: str, results: list[object]) -> dict[str, Any] | None:
    normalized_query = _normalize_text(query)
    if not normalized_query:
        return None
    for item in results:
        if not isinstance(item, dict):
            continue
        name = _normalize_text(str(item.get("name", "")))
        pdes = _normalize_text(str(item.get("pdes", "")))
        aliases = {
            _normalize_text(alias)
            for alias in (item.get("alias") or [])
            if isinstance(alias, str)
        }
        if normalized_query == name or normalized_query == pdes or normalized_query in aliases:
            return item
    return None


def _format_candidate(entry: dict[str, Any]) -> str:
    name = str(entry.get("name", "")).strip()
    pdes = str(entry.get("pdes", "")).strip()
    spkid = str(entry.get("spkid", "")).strip()
    type_name = str(entry.get("type", "")).strip()
    parts = [part for part in (name, pdes, f"spkid={spkid}" if spkid else "", type_name) if part]
    return " | ".join(parts)


def _extract_altaz(rows: list[list[str]]) -> tuple[float | None, float | None]:
    for row in rows:
        numeric_values: list[float] = []
        for value in row:
            try:
                numeric_values.append(float(str(value).strip()))
            except (TypeError, ValueError):
                continue
        if len(numeric_values) >= 2:
            return numeric_values[-1], numeric_values[-2]
    return None, None


def _build_command(target: dict[str, Any]) -> str:
    spkid = str(target.get("spkid", "")).strip()
    if spkid:
        return f"DES={spkid};"
    pdes = str(target.get("pdes", "")).strip()
    if pdes:
        return f"DES={pdes};"
    name = str(target.get("name", "")).strip()
    if name:
        return name
    return ""


def _fetch_lookup(query: str, group: str | None) -> dict[str, Any]:
    payload = fetch_horizons_lookup(query, group=group or DEFAULT_GROUP)
    if not isinstance(payload, dict):
        raise RuntimeError("Horizons lookup returned an invalid payload")
    return payload


def _resolve_target(query: str, group: str | None) -> dict[str, Any] | None:
    payload = _fetch_lookup(query, group)
    result = payload.get("result")
    if not isinstance(result, list) or not result:
        return None
    exact = _pick_exact_match(query, result)
    if exact is not None:
        return exact
    first = result[0]
    if isinstance(first, dict):
        return first
    return None


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="jpl_horizons_observer.py",
        description="Resolve a Horizons target and fetch one observer ephemeris sample.",
    )
    parser.add_argument("query", help="Search string passed to the Horizons lookup API.")
    parser.add_argument(
        "--group",
        choices=sorted(ALLOWED_GROUPS),
        default=DEFAULT_GROUP,
        help="Restrict the lookup to a Horizons object group.",
    )
    parser.add_argument("--lat", type=float, required=True, help="Observer latitude in degrees.")
    parser.add_argument("--lon", type=float, required=True, help="Observer longitude in degrees.")
    parser.add_argument("--height-m", type=float, default=0.0, help="Observer height above sea level in meters.")
    parser.add_argument(
        "--time",
        dest="time_utc",
        help="Observation time in UTC ISO-8601 format. Defaults to the current minute.",
    )
    parser.add_argument("--json", action="store_true", help="Emit the structured payload as JSON.")
    parser.add_argument(
        "--compact",
        action="store_true",
        help="Use terse one-line output for the resolved target and sample.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    query = str(args.query).strip()
    if not query:
        print("Query must not be empty.", file=sys.stderr)
        return 2

    try:
        target_time_utc = _parse_utc_time(args.time_utc)
    except ValueError as exc:
        print(f"Invalid --time value: {exc}", file=sys.stderr)
        return 2

    try:
        lookup_payload = _fetch_lookup(query, str(args.group).strip() if args.group else None)
    except Exception as exc:
        print(f"Horizons lookup request failed: {exc}", file=sys.stderr)
        return 1

    result = lookup_payload.get("result")
    if not isinstance(result, list) or not result:
        print("No matching Horizons object found.")
        return 0

    target = _pick_exact_match(query, result)
    if target is None and isinstance(result[0], dict):
        target = result[0]
    if target is None:
        print("No usable Horizons target was returned.")
        return 1

    target_spkid = str(target.get("spkid", "")).strip()
    command = _build_command(target)
    if not command:
        print("Horizons lookup did not provide a usable target command.")
        return 1

    try:
        rows = fetch_horizons_observer_csv(
            command,
            target_time_utc=target_time_utc,
            observer_lat=float(args.lat),
            observer_lon=float(args.lon),
            observer_height_m=float(args.height_m),
        )
    except Exception as exc:
        print(f"Horizons observer request failed: {exc}", file=sys.stderr)
        return 1

    alt_deg, az_deg = _extract_altaz(rows)
    if alt_deg is None or az_deg is None:
        print("Horizons observer table did not contain an alt/az sample.")
        return 1

    resolved_name = str(target.get("name", "")).strip() or query
    if args.json:
        payload = {
            "query": query,
            "group": args.group,
            "time_utc": target_time_utc.astimezone(timezone.utc).isoformat(),
            "observer": {
                "lat_deg": float(args.lat),
                "lon_deg": float(args.lon),
                "height_m": float(args.height_m),
            },
            "lookup": lookup_payload,
            "target": target,
            "observer_sample": {
                "alt_deg": float(alt_deg),
                "az_deg": float(az_deg),
                "rows": rows,
            },
        }
        print(json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True))
        return 0

    if args.compact:
        print(
            f"{resolved_name} | spkid={target_spkid} | cmd={command} | "
            f"alt={float(alt_deg):.3f} | az={float(az_deg):.3f} | "
            f"utc={target_time_utc.astimezone(timezone.utc).isoformat()}"
        )
        return 0

    print(resolved_name)
    print(f"spkid: {target_spkid}")
    print(f"command: {command}")
    print(f"type: {str(target.get('type', '')).strip()}")
    print(f"observer: {float(args.lat):.6f}, {float(args.lon):.6f}, {float(args.height_m):.1f} m")
    print(f"time_utc: {target_time_utc.astimezone(timezone.utc).isoformat()}")
    print(f"alt_deg: {float(alt_deg):.6f}")
    print(f"az_deg: {float(az_deg):.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
