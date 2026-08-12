#!/usr/bin/env python3
"""Inspect GMN meteor loading and celestial-coordinate round trips.

This development helper exercises the phase-one GMN pipeline without starting
the GUI. It writes machine-readable JSONL and can create a simple red polar SVG
of the event-time visible paths.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TextIO

import astropy.time
import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.meteors.constants import GMN_CANDIDATE_RADIUS_KM, GMN_WINDOW  # noqa: E402
from zstarview.meteors.projection import (  # noqa: E402
    _clip_segment_to_geometric_horizon,
    _earth_location_xyz_m,
    _enu_basis,
    _geodetic_xyz_m,
    _within_candidate_radius,
    project_meteor_observations_to_celestial,
)
from zstarview.meteors.repository import GmnMeteorRepository  # noqa: E402
from zstarview.meteors.types import (  # noqa: E402
    CelestialMeteorTrail,
    GmnLoadResult,
    MeteorObservation,
)

SCHEMA = "zstarview.gmn-meteor-probe.v1"
DEFAULT_SVG_SIZE = 900
TRAIL_COLOR = "#f04848"


@dataclass(frozen=True, slots=True)
class ProbeRun:
    display_time_utc: datetime
    window_start_utc: datetime
    observer_lat: float
    observer_lon: float
    observer_height_m: float
    loaded_count: int
    candidate_count: int
    visible_count: int
    clipped_count: int
    source_files: tuple[str, ...]
    unavailable_files: tuple[str, ...]
    used_stale_index: bool
    used_stale_files: bool
    records: tuple[dict[str, Any], ...]


def parse_utc(value: str) -> datetime:
    text = value.strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    parsed = datetime.fromisoformat(text)
    if parsed.tzinfo is None:
        raise ValueError("--time must include a UTC offset or Z suffix")
    return parsed.astimezone(timezone.utc)


def build_probe_run(
    loaded: GmnLoadResult,
    *,
    display_time_utc: datetime,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    max_records: int | None = None,
) -> ProbeRun:
    display_time = _normalize_utc(display_time_utc)
    candidate_count = 0
    visible_count = 0
    clipped_count = 0
    records: list[dict[str, Any]] = []
    for observation in loaded.observations:
        candidate = _within_candidate_radius(
            observation,
            observer_lat=observer_lat,
            observer_lon=observer_lon,
            radius_km=GMN_CANDIDATE_RADIUS_KM,
        )
        if not candidate:
            continue
        candidate_count += 1
        trails = project_meteor_observations_to_celestial(
            (observation,),
            observer_lat=observer_lat,
            observer_lon=observer_lon,
            observer_height_m=observer_height_m,
        )
        if not trails:
            continue
        visible_count += 1
        record = build_trail_diagnostic(
            observation,
            trails[0],
            observer_lat=observer_lat,
            observer_lon=observer_lon,
            observer_height_m=observer_height_m,
        )
        if bool(record["horizon_clipped"]):
            clipped_count += 1
        if max_records is None or len(records) < max_records:
            records.append(record)
    return ProbeRun(
        display_time_utc=display_time,
        window_start_utc=display_time - GMN_WINDOW,
        observer_lat=float(observer_lat),
        observer_lon=float(observer_lon),
        observer_height_m=float(observer_height_m),
        loaded_count=len(loaded.observations),
        candidate_count=candidate_count,
        visible_count=visible_count,
        clipped_count=clipped_count,
        source_files=loaded.source_files,
        unavailable_files=loaded.unavailable_files,
        used_stale_index=loaded.used_stale_index,
        used_stale_files=loaded.used_stale_files,
        records=tuple(records),
    )


def build_trail_diagnostic(
    observation: MeteorObservation,
    trail: CelestialMeteorTrail,
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
) -> dict[str, Any]:
    observer = EarthLocation(
        lat=float(observer_lat) * u.deg,
        lon=float(observer_lon) * u.deg,
        height=float(observer_height_m) * u.m,
    )
    observer_xyz = _earth_location_xyz_m(observer)
    _, _, up = _enu_basis(
        observation_lat_deg=observer_lat,
        observation_lon_deg=observer_lon,
    )
    begin_xyz = _geodetic_xyz_m(
        observation.begin_lat_deg,
        observation.begin_lon_deg,
        observation.begin_height_km,
    )
    end_xyz = _geodetic_xyz_m(
        observation.end_lat_deg,
        observation.end_lon_deg,
        observation.end_height_km,
    )
    original_begin = _vector_to_altaz(
        begin_xyz - observer_xyz,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
    )
    original_end = _vector_to_altaz(
        end_xyz - observer_xyz,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
    )
    clipped = _clip_segment_to_geometric_horizon(
        begin_xyz,
        end_xyz,
        observer_xyz=observer_xyz,
        up=up,
    )
    if clipped is None:
        raise ValueError("visible trail unexpectedly has no clipped segment")
    clipped_begin, clipped_end = clipped
    expected_begin = _vector_to_altaz(
        clipped_begin - observer_xyz,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
    )
    expected_end = _vector_to_altaz(
        clipped_end - observer_xyz,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
    )
    restored_begin = _icrs_to_event_altaz(
        trail.begin_ra_deg,
        trail.begin_dec_deg,
        observation=observation,
        observer=observer,
    )
    restored_end = _icrs_to_event_altaz(
        trail.end_ra_deg,
        trail.end_dec_deg,
        observation=observation,
        observer=observer,
    )
    begin_error = _angular_error_deg(expected_begin, restored_begin)
    end_error = _angular_error_deg(expected_end, restored_end)
    horizon_clipped = original_begin[0] < 0.0 or original_end[0] < 0.0
    return {
        "schema": SCHEMA,
        "type": "trail",
        "trajectory_id": observation.trajectory_id,
        "beginning_utc": observation.beginning_utc.isoformat(),
        "horizon_clipped": horizon_clipped,
        "source_geodetic": {
            "begin": _geodetic_dict(
                observation.begin_lat_deg,
                observation.begin_lon_deg,
                observation.begin_height_km,
            ),
            "end": _geodetic_dict(
                observation.end_lat_deg,
                observation.end_lon_deg,
                observation.end_height_km,
            ),
        },
        "original_event_altaz_deg": {
            "begin": _altaz_dict(original_begin),
            "end": _altaz_dict(original_end),
        },
        "clipped_event_altaz_deg": {
            "begin": _altaz_dict(expected_begin),
            "end": _altaz_dict(expected_end),
        },
        "fixed_icrs_deg": {
            "begin": {"ra": trail.begin_ra_deg, "dec": trail.begin_dec_deg},
            "end": {"ra": trail.end_ra_deg, "dec": trail.end_dec_deg},
        },
        "roundtrip_event_altaz_deg": {
            "begin": _altaz_dict(restored_begin),
            "end": _altaz_dict(restored_end),
        },
        "roundtrip_error_deg": {
            "begin": begin_error,
            "end": end_error,
            "maximum": max(begin_error, end_error),
        },
        "duration_s": observation.duration_s,
        "peak_abs_magnitude": observation.peak_abs_magnitude,
        "shower_code": observation.shower_code,
    }


def build_summary_record(run: ProbeRun) -> dict[str, Any]:
    maximum_error = max(
        (
            float(record["roundtrip_error_deg"]["maximum"])
            for record in run.records
        ),
        default=0.0,
    )
    return {
        "schema": SCHEMA,
        "type": "summary",
        "display_time_utc": run.display_time_utc.isoformat(),
        "window_start_utc": run.window_start_utc.isoformat(),
        "observer": {
            "lat": run.observer_lat,
            "lon": run.observer_lon,
            "height_m": run.observer_height_m,
        },
        "counts": {
            "loaded": run.loaded_count,
            "within_1000_km": run.candidate_count,
            "visible": run.visible_count,
            "horizon_clipped": run.clipped_count,
            "records_written": len(run.records),
        },
        "source_files": list(run.source_files),
        "unavailable_files": list(run.unavailable_files),
        "used_stale_index": run.used_stale_index,
        "used_stale_files": run.used_stale_files,
        "maximum_roundtrip_error_deg": maximum_error,
    }


def build_svg(run: ProbeRun, *, size: int = DEFAULT_SVG_SIZE) -> str:
    width = max(320, int(size))
    center = width / 2.0
    radius = width * 0.42
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width}" height="{width}" viewBox="0 0 {width} {width}">',
        f'<rect width="{width}" height="{width}" fill="#101014"/>',
        f'<circle cx="{center:.2f}" cy="{center:.2f}" r="{radius:.2f}" '
        'fill="none" stroke="#707078" stroke-width="1.5"/>',
    ]
    for altitude in (30.0, 60.0):
        ring_radius = radius * (90.0 - altitude) / 90.0
        lines.append(
            f'<circle cx="{center:.2f}" cy="{center:.2f}" r="{ring_radius:.2f}" '
            'fill="none" stroke="#383840" stroke-width="1"/>'
        )
    for azimuth, label in ((0.0, "N"), (90.0, "E"), (180.0, "S"), (270.0, "W")):
        x, y = polar_svg_xy(0.0, azimuth, center=center, radius=radius)
        lines.append(
            f'<text x="{x:.2f}" y="{y:.2f}" fill="#c8c8ce" '
            f'font-family="sans-serif" font-size="14">{label}</text>'
        )
    for record in run.records:
        event_altaz = record["roundtrip_event_altaz_deg"]
        begin = event_altaz["begin"]
        end = event_altaz["end"]
        x1, y1 = polar_svg_xy(begin["alt"], begin["az"], center=center, radius=radius)
        x2, y2 = polar_svg_xy(end["alt"], end["az"], center=center, radius=radius)
        title = html.escape(str(record["trajectory_id"]), quote=True)
        lines.append(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{TRAIL_COLOR}" stroke-opacity="0.72" stroke-width="1.6" '
            'stroke-linecap="round">'
            f"<title>{title}</title></line>"
        )
        # Endpoint dots keep very short, near-zenith trails visible in the
        # diagnostic plot without changing their projected geometry.
        for x, y in ((x1, y1), (x2, y2)):
            lines.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.2" '
                f'fill="{TRAIL_COLOR}" fill-opacity="0.88">'
                f"<title>{title}</title></circle>"
            )
    summary = (
        f"loaded={run.loaded_count} candidate={run.candidate_count} "
        f"visible={run.visible_count} clipped={run.clipped_count}"
    )
    lines.append(
        f'<text x="20" y="28" fill="#e0e0e6" font-family="monospace" '
        f'font-size="14">{html.escape(summary)}</text>'
    )
    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def polar_svg_xy(
    alt_deg: float,
    az_deg: float,
    *,
    center: float,
    radius: float,
) -> tuple[float, float]:
    radial = radius * (90.0 - max(0.0, min(90.0, float(alt_deg)))) / 90.0
    azimuth = math.radians(float(az_deg))
    return center + radial * math.sin(azimuth), center - radial * math.cos(azimuth)


def _vector_to_altaz(
    vector_m: np.ndarray,
    *,
    observer_lat: float,
    observer_lon: float,
) -> tuple[float, float]:
    east, north, up = _enu_basis(
        observation_lat_deg=observer_lat,
        observation_lon_deg=observer_lon,
    )
    east_m = float(np.dot(vector_m, east))
    north_m = float(np.dot(vector_m, north))
    up_m = float(np.dot(vector_m, up))
    return (
        math.degrees(math.atan2(up_m, math.hypot(east_m, north_m))),
        math.degrees(math.atan2(east_m, north_m)) % 360.0,
    )


def _icrs_to_event_altaz(
    ra_deg: float,
    dec_deg: float,
    *,
    observation: MeteorObservation,
    observer: EarthLocation,
) -> tuple[float, float]:
    frame = AltAz(
        obstime=astropy.time.Time(observation.beginning_utc),
        location=observer,
    )
    altaz = SkyCoord(
        ra=float(ra_deg) * u.deg,
        dec=float(dec_deg) * u.deg,
        frame="icrs",
    ).transform_to(frame)
    return float(altaz.alt.deg), float(altaz.az.deg) % 360.0


def _angular_error_deg(
    expected: tuple[float, float],
    actual: tuple[float, float],
) -> float:
    expected_alt, expected_az = map(math.radians, expected)
    actual_alt, actual_az = map(math.radians, actual)
    cosine = (
        math.sin(expected_alt) * math.sin(actual_alt)
        + math.cos(expected_alt)
        * math.cos(actual_alt)
        * math.cos(expected_az - actual_az)
    )
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _geodetic_dict(lat: float, lon: float, height_km: float) -> dict[str, float]:
    return {"lat": float(lat), "lon": float(lon), "height_km": float(height_km)}


def _altaz_dict(value: tuple[float, float]) -> dict[str, float]:
    return {"alt": float(value[0]), "az": float(value[1])}


def _normalize_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _offline_fetcher(url: str, *, timeout_s: float) -> str:
    raise OSError(f"offline mode has no cached GMN data for {url}")


def _open_jsonl_output(path: Path | None) -> tuple[TextIO, bool]:
    if path is None:
        return sys.stdout, False
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.open("w", encoding="utf-8"), True


def _write_jsonl(run: ProbeRun, stream: TextIO) -> None:
    payloads = (build_summary_record(run), *run.records)
    for payload in payloads:
        stream.write(json.dumps(payload, ensure_ascii=True, sort_keys=True) + "\n")


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--time",
        required=True,
        help="Display time as ISO 8601 with UTC offset or Z suffix.",
    )
    parser.add_argument("--lat", required=True, type=float, help="Observer latitude in degrees.")
    parser.add_argument("--lon", required=True, type=float, help="Observer longitude in degrees.")
    parser.add_argument("--height-m", type=float, default=0.0, help="Observer WGS84 height in metres.")
    parser.add_argument("--cache-dir", type=Path, help="Optional isolated GMN cache directory.")
    parser.add_argument("--offline", action="store_true", help="Use cache only; do not access GMN.")
    parser.add_argument("--output-jsonl", type=Path, help="Write JSONL here instead of stdout.")
    parser.add_argument("--output-svg", type=Path, help="Write an event-time polar SVG preview.")
    parser.add_argument("--svg-size", type=int, default=DEFAULT_SVG_SIZE)
    parser.add_argument(
        "--max-records",
        type=int,
        default=500,
        help="Maximum detail records and SVG trails; counts still cover all observations.",
    )
    parser.add_argument("--timeout", type=float, default=30.0, help="HTTP timeout in seconds.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_argument_parser().parse_args(argv)
    try:
        display_time = parse_utc(args.time)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if not -90.0 <= args.lat <= 90.0 or not -180.0 <= args.lon <= 180.0:
        print("error: observer latitude or longitude is out of range", file=sys.stderr)
        return 2
    if args.max_records < 0:
        print("error: --max-records must be non-negative", file=sys.stderr)
        return 2
    repository_kwargs: dict[str, Any] = {"timeout_s": args.timeout}
    if args.cache_dir is not None:
        repository_kwargs["cache_root"] = args.cache_dir
    if args.offline:
        repository_kwargs["fetcher"] = _offline_fetcher
    repository = GmnMeteorRepository(**repository_kwargs)
    try:
        loaded = repository.load_window(display_time - GMN_WINDOW, display_time)
        run = build_probe_run(
            loaded,
            display_time_utc=display_time,
            observer_lat=args.lat,
            observer_lon=args.lon,
            observer_height_m=args.height_m,
            max_records=args.max_records,
        )
    except Exception as exc:
        print(f"error: GMN probe failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    stream, should_close = _open_jsonl_output(args.output_jsonl)
    try:
        _write_jsonl(run, stream)
    finally:
        if should_close:
            stream.close()
    if args.output_svg is not None:
        args.output_svg.parent.mkdir(parents=True, exist_ok=True)
        args.output_svg.write_text(build_svg(run, size=args.svg_size), encoding="utf-8")
    summary = build_summary_record(run)["counts"]
    print(
        "GMN probe complete: "
        f"loaded={summary['loaded']} candidate={summary['within_1000_km']} "
        f"visible={summary['visible']} clipped={summary['horizon_clipped']}",
        file=sys.stderr,
    )
    if args.output_jsonl is not None:
        print(f"JSONL: {args.output_jsonl}", file=sys.stderr)
    if args.output_svg is not None:
        print(f"SVG: {args.output_svg}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
