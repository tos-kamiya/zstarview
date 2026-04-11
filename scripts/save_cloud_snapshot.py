#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Save a realtime cloud disc image for a specified city or coordinates.

Examples:
  python scripts/save_cloud_snapshot.py Tokyo
  python scripts/save_cloud_snapshot.py JP/Tokyo --output tokyo_cloud.png
  python scripts/save_cloud_snapshot.py "35.68;139.76" --radius-px 320
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.clouddisc import CloudDisc, CloudDiscConfig, CloudDiscError  # noqa: E402
from zstarview.paths import (  # noqa: E402
    CACHE_PATH,
    CITY_ADMIN1_CODES_FILE,
    CITY_COORD_FILE,
    CLOUD_SHELL_KM,
)
from zstarview.utils.resolve_city import (  # noqa: E402
    CityRec,
    load_admin1_names,
    resolve_city,
    resolve_city_by_geonameid,
    resolve_city_by_name,
)


logger = logging.getLogger("save_cloud_snapshot")


@dataclass(frozen=True)
class ResolvedLocation:
    label: str
    lat: float
    lon: float
    timezone_name: str
    geonameid: int


def _parse_coord_token(token: str, allowed_dirs: str) -> float:
    text = token.strip().upper()
    found = {ch for ch in text if ch in "NSEW"}
    allowed = set(allowed_dirs)
    if found and not found.issubset(allowed):
        raise ValueError(f"Invalid direction in {token!r}. Expected one of {sorted(allowed)}")
    sign = -1.0 if ("S" in found or "W" in found) else 1.0
    num = re.sub(r"[^0-9.\-]", "", text)
    if not num:
        raise ValueError(f"No numeric value found in {token!r}")
    val = float(num)
    return val if val < 0.0 else val * sign


def _resolve_location(query: str) -> ResolvedLocation:
    q = query.strip()
    if not q:
        raise ValueError("Location must not be empty")

    # Direct coordinates: "lat;lon"
    if ";" in q:
        parts = [p.strip() for p in q.split(";")]
        if len(parts) != 2:
            raise ValueError("Coordinate format must be 'lat;lon'")
        lat = _parse_coord_token(parts[0], "NS")
        lon = _parse_coord_token(parts[1], "EW")
        if not (-90.0 <= lat <= 90.0):
            raise ValueError(f"Latitude out of range: {lat}")
        if not (-180.0 <= lon <= 180.0):
            raise ValueError(f"Longitude out of range: {lon}")
        return ResolvedLocation(
            label=f"Lat{lat:.4f}_Lon{lon:.4f}",
            lat=lat,
            lon=lon,
            timezone_name="UTC",
            geonameid=0,
        )

    admin1_map = load_admin1_names(CITY_ADMIN1_CODES_FILE)
    rec: Optional[CityRec] = None
    if re.fullmatch(r"\d+", q):
        rec = resolve_city_by_geonameid(int(q), CITY_COORD_FILE, admin1_map)
    elif "/" in q:
        recs = resolve_city(q, CITY_COORD_FILE, admin1_map)
        if recs:
            rec = recs[0]
    else:
        recs = resolve_city_by_name(q, CITY_COORD_FILE, admin1_map)
        if recs:
            rec = recs[0]

    if rec is None:
        raise ValueError(f"Could not resolve location: {q!r}")

    return ResolvedLocation(
        label=f"{rec.cc}_{rec.name}".replace(" ", "_"),
        lat=rec.lat,
        lon=rec.lon,
        timezone_name=rec.tz,
        geonameid=rec.geonameid,
    )


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Save a realtime cloud disc image for a city or coordinates."
    )
    p.add_argument(
        "location",
        help=(
            "City query (e.g., Tokyo, JP/Tokyo, geonameid) "
            "or coordinates 'lat;lon' (e.g., 35.68;139.76)."
        ),
    )
    p.add_argument("--output", "-o", default="", help="Output PNG path. Auto-generated if omitted.")
    p.add_argument("--meta", default="", help="Optional path to write JSON metadata.")
    p.add_argument("--alt", type=float, default=90.0, help="View altitude in degrees (default: 90).")
    p.add_argument("--az", type=float, default=180.0, help="View azimuth in degrees (default: 180).")
    p.add_argument("--radius-px", type=int, default=256, help="Cloud disc radius in pixels (default: 256).")
    p.add_argument(
        "--sat-priority",
        default="AUTO",
        help="Satellite priority CSV (e.g., AUTO or G16,G18,HIMAWARI).",
    )
    p.add_argument(
        "--search-back-minutes",
        type=int,
        default=120,
        help="Lookback window for satellite data (default: 120).",
    )
    p.add_argument("--cache-dir", default=CACHE_PATH, help=f"Cloud cache directory (default: {CACHE_PATH}).")
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        help="Logging level (default: INFO).",
    )
    return p


def main() -> int:
    args = _build_parser().parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    try:
        loc = _resolve_location(args.location)
    except Exception as e:
        logger.error("%s", e)
        return 2

    sat_priority = tuple(
        t.strip().upper() for t in args.sat_priority.split(",") if t.strip()
    ) or ("AUTO",)
    cfg = CloudDiscConfig(
        cache_dir=Path(args.cache_dir),
        sat_priority=sat_priority,
        search_back_minutes=max(10, int(args.search_back_minutes)),
    )
    clouddisc = CloudDisc(cfg)

    alt = max(-5.0, min(90.0, float(args.alt)))
    az = float(args.az) % 360.0
    radius_px = max(64, int(args.radius_px))

    logger.info(
        "Rendering cloud image for %s (lat=%.4f lon=%.4f alt=%.1f az=%.1f)",
        loc.label,
        loc.lat,
        loc.lon,
        alt,
        az,
    )
    try:
        img, meta = clouddisc.render_now(
            lat=loc.lat,
            lon=loc.lon,
            alt=alt,
            az=az,
            radius_px=radius_px,
            cloud_shell_km=CLOUD_SHELL_KM,
        )
    except CloudDiscError as e:
        logger.error("CloudDisc failed: %s", e)
        return 1
    except Exception as e:
        logger.exception("Unexpected error while rendering: %s", e)
        return 1

    if args.output:
        out_png = Path(args.output)
    else:
        stamp = meta.time_utc.strftime("%Y%m%dT%H%M%SZ")
        out_png = Path(f"cloud_{loc.label}_{stamp}.png")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    img.save(out_png)

    info = {
        "location_query": args.location,
        "label": loc.label,
        "geonameid": loc.geonameid,
        "lat": loc.lat,
        "lon": loc.lon,
        "timezone": loc.timezone_name,
        "alt": alt,
        "az": az,
        "radius_px": radius_px,
        "satellite": meta.satellite,
        "product": meta.product,
        "source_time_utc": meta.time_utc.isoformat(),
        "src_paths": [str(p) for p in meta.src_paths],
        "output_png": str(out_png),
    }

    if args.meta:
        out_meta = Path(args.meta)
        out_meta.parent.mkdir(parents=True, exist_ok=True)
        out_meta.write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps(info, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
