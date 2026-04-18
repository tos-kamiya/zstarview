#!/usr/bin/env python3
"""
jpl_horizons_trajectory.py
Version: 0.1.0

A small standalone command-line tool that resolves a JPL Horizons object
and fetches a short observer trajectory for a specific location.

This script is the third step in the small-body prototype sequence:

1. search the SBDB
2. resolve the Horizons target
3. fetch a 1-hour observer trajectory

It uses repeated single-time observer queries so the prototype stays close to
the already working observer sample.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.satellites.fetch import fetch_horizons_lookup, fetch_horizons_observer_csv  # noqa: E402


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


def _resolve_target(query: str, group: str | None) -> tuple[dict[str, Any], dict[str, Any]]:
    lookup_payload = fetch_horizons_lookup(query, group=group or DEFAULT_GROUP)
    if not isinstance(lookup_payload, dict):
        raise RuntimeError("Horizons lookup returned an invalid payload")
    result = lookup_payload.get("result")
    if not isinstance(result, list) or not result:
        raise RuntimeError("No matching Horizons object found")
    target = _pick_exact_match(query, result)
    if target is None and isinstance(result[0], dict):
        target = result[0]
    if target is None:
        raise RuntimeError("No usable Horizons target was returned")
    return lookup_payload, target


def _sample_trajectory(
    command: str,
    *,
    start_time_utc: datetime,
    duration_minutes: int,
    step_minutes: int,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
) -> list[dict[str, Any]]:
    samples: list[dict[str, Any]] = []
    offset = 0
    while offset <= duration_minutes:
        sample_time = start_time_utc + timedelta(minutes=offset)
        rows = fetch_horizons_observer_csv(
            command,
            target_time_utc=sample_time,
            observer_lat=observer_lat,
            observer_lon=observer_lon,
            observer_height_m=observer_height_m,
        )
        alt_deg, az_deg = _extract_altaz(rows)
        if alt_deg is None or az_deg is None:
            raise RuntimeError(f"Missing alt/az sample at {sample_time.isoformat()}")
        samples.append(
            {
                "time_utc": sample_time.astimezone(timezone.utc).isoformat(),
                "alt_deg": float(alt_deg),
                "az_deg": float(az_deg),
                "rows": rows,
            }
        )
        offset += step_minutes
    return samples


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="jpl_horizons_trajectory.py",
        description="Resolve a Horizons target and fetch a short observer trajectory.",
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
        help="Start time in UTC ISO-8601 format. Defaults to the current minute.",
    )
    parser.add_argument(
        "--duration-minutes",
        type=int,
        default=60,
        help="How many minutes of trajectory to sample, inclusive of the start time.",
    )
    parser.add_argument(
        "--step-minutes",
        type=int,
        default=15,
        help="Sampling step in minutes.",
    )
    parser.add_argument("--json", action="store_true", help="Emit the structured payload as JSON.")
    parser.add_argument(
        "--compact",
        action="store_true",
        help="Use terse one-line output for each trajectory sample.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    query = str(args.query).strip()
    if not query:
        print("Query must not be empty.", file=sys.stderr)
        return 2
    if args.duration_minutes < 0:
        print("--duration-minutes must be non-negative.", file=sys.stderr)
        return 2
    if args.step_minutes <= 0:
        print("--step-minutes must be positive.", file=sys.stderr)
        return 2

    try:
        start_time_utc = _parse_utc_time(args.time_utc)
    except ValueError as exc:
        print(f"Invalid --time value: {exc}", file=sys.stderr)
        return 2

    try:
        lookup_payload, target = _resolve_target(query, str(args.group).strip() if args.group else None)
    except Exception as exc:
        print(f"Horizons lookup request failed: {exc}", file=sys.stderr)
        return 1

    target_spkid = str(target.get("spkid", "")).strip()
    command = _build_command(target)
    if not command:
        print("Horizons lookup did not provide a usable target command.", file=sys.stderr)
        return 1

    try:
        samples = _sample_trajectory(
            command,
            start_time_utc=start_time_utc,
            duration_minutes=int(args.duration_minutes),
            step_minutes=int(args.step_minutes),
            observer_lat=float(args.lat),
            observer_lon=float(args.lon),
            observer_height_m=float(args.height_m),
        )
    except Exception as exc:
        print(f"Horizons trajectory request failed: {exc}", file=sys.stderr)
        return 1

    resolved_name = str(target.get("name", "")).strip() or query
    if args.json:
        payload = {
            "query": query,
            "group": args.group,
            "start_time_utc": start_time_utc.astimezone(timezone.utc).isoformat(),
            "duration_minutes": int(args.duration_minutes),
            "step_minutes": int(args.step_minutes),
            "observer": {
                "lat_deg": float(args.lat),
                "lon_deg": float(args.lon),
                "height_m": float(args.height_m),
            },
            "lookup": lookup_payload,
            "target": target,
            "trajectory": samples,
        }
        print(json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True))
        return 0

    if args.compact:
        for sample in samples:
            print(
                f"{resolved_name} | spkid={target_spkid} | cmd={command} | "
                f"utc={sample['time_utc']} | alt={sample['alt_deg']:.3f} | az={sample['az_deg']:.3f}"
            )
        return 0

    print(resolved_name)
    print(f"spkid: {target_spkid}")
    print(f"command: {command}")
    print(f"observer: {float(args.lat):.6f}, {float(args.lon):.6f}, {float(args.height_m):.1f} m")
    print(f"start_time_utc: {start_time_utc.astimezone(timezone.utc).isoformat()}")
    print(f"duration_minutes: {int(args.duration_minutes)}")
    print(f"step_minutes: {int(args.step_minutes)}")
    for sample in samples:
        print(
            f"- {sample['time_utc']}: "
            f"alt={sample['alt_deg']:.6f}, az={sample['az_deg']:.6f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
