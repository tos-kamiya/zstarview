#!/usr/bin/env python3
"""Render water footprints as outer-ring outlines in a simple SVG plan view.

This script is intended for local debugging of cache snapshots or compact
footprint JSON files. It draws only the outer rings of each footprint as
stroked lines and skips polygon fills entirely.
"""

from __future__ import annotations

import argparse
import html
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = REPO_ROOT / "dev-samples" / "water_footprints_outline_preview.svg"

BACKGROUND = "#f4eadb"
TEXT_FILL = "#0f172a"
GRID_STROKE = "#cbd5e1"
ELE_STROKE = "#0f4c5c"
NO_ELE_STROKE = "#8c5a00"


@dataclass(frozen=True)
class Footprint:
    index: int
    water_id: str
    kind: str
    outer_rings_lonlat: tuple[tuple[tuple[float, float], ...], ...]
    inner_rings_lonlat: tuple[tuple[tuple[float, float], ...], ...]
    source: str
    tags: dict[str, str]

    @property
    def has_explicit_height(self) -> bool:
        return "ele" in self.tags or "water_level" in self.tags


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render water footprints from a zstarview cache JSON into an outline SVG.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Input JSON path containing a water-overlay cache snapshot or footprint list.",
    )
    parser.add_argument(
        "--input-cache",
        type=Path,
        help="Input zstarview cache JSON path. This is equivalent to --input for cache payloads.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"SVG output path. Default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument("--width", type=int, default=1200, help="SVG width in pixels.")
    parser.add_argument("--height", type=int, default=1200, help="SVG height in pixels.")
    parser.add_argument(
        "--padding-ratio",
        type=float,
        default=0.06,
        help="Extra margin around the footprint bounds, as a fraction of the data span.",
    )
    parser.add_argument(
        "--show-legend",
        action="store_true",
        help="Draw a small legend explaining the outline colors.",
    )
    parser.add_argument(
        "--show-labels",
        action="store_true",
        help="Draw footprint indices and water IDs at approximate centroids.",
    )
    parser.add_argument(
        "--ele-only",
        action="store_true",
        help="Skip footprints that do not carry explicit ele or water_level tags.",
    )
    parser.add_argument(
        "--background",
        default=BACKGROUND,
        help="SVG background color.",
    )
    args = parser.parse_args()
    if args.input is None and args.input_cache is None:
        parser.error("one of --input or --input-cache is required")
    return args


def _load_payload(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("Input JSON must be an object")
    return payload


def _load_rings(value: object) -> tuple[tuple[tuple[float, float], ...], ...]:
    if not isinstance(value, list):
        return ()
    rings: list[tuple[tuple[float, float], ...]] = []
    for ring in value:
        if not isinstance(ring, list):
            continue
        points: list[tuple[float, float]] = []
        for point in ring:
            if not isinstance(point, list) or len(point) != 2:
                continue
            try:
                lon = float(point[0])
                lat = float(point[1])
            except (TypeError, ValueError):
                continue
            points.append((lon, lat))
        if len(points) >= 4:
            rings.append(tuple(points))
    return tuple(rings)


def _load_tags(value: object) -> dict[str, str]:
    if not isinstance(value, dict):
        return {}
    tags: dict[str, str] = {}
    for key, item in value.items():
        if isinstance(key, str) and isinstance(item, str):
            tags[key] = item
    return tags


def _is_coastline_footprint(item: dict[str, Any]) -> bool:
    kind = str(item.get("kind", ""))
    source = str(item.get("source", ""))
    if kind == "coastline":
        return True
    return source == "coastline"


def load_footprints(path: Path) -> tuple[Footprint, ...]:
    payload = _load_payload(path)
    footprints: list[Footprint] = []
    items = payload.get("footprints")
    if not isinstance(items, list):
        items = payload.get("water_polygons")
    if not isinstance(items, list):
        raise ValueError(f"No footprints were found in {path}")
    for index, item in enumerate(items):
        if not isinstance(item, dict):
            continue
        if not _is_coastline_footprint(item):
            continue
        outer = _load_rings(item.get("outer_rings_lonlat") or item.get("outer_rings"))
        inner = _load_rings(item.get("inner_rings_lonlat") or item.get("inner_rings"))
        if not outer:
            continue
        footprints.append(
            Footprint(
                index=index,
                water_id=str(item.get("water_id") or item.get("osm_id") or f"footprint/{index}"),
                kind=str(item.get("kind", "water_polygon")),
                outer_rings_lonlat=outer,
                inner_rings_lonlat=inner,
                source=str(item.get("source", "")),
                tags=_load_tags(item.get("tags")),
            )
        )
    if not footprints:
        raise ValueError(f"No footprints were found in {path}")
    return tuple(footprints)


def _ring_bounds(ring: tuple[tuple[float, float], ...]) -> tuple[float, float, float, float]:
    lons = [lon for lon, _lat in ring]
    lats = [lat for _lon, lat in ring]
    return min(lons), min(lats), max(lons), max(lats)


def _footprint_bounds(footprint: Footprint) -> tuple[float, float, float, float]:
    bounds: list[tuple[float, float, float, float]] = []
    for ring in footprint.outer_rings_lonlat:
        bounds.append(_ring_bounds(ring))
    if not bounds:
        return 0.0, 0.0, 0.0, 0.0
    min_lon = min(item[0] for item in bounds)
    min_lat = min(item[1] for item in bounds)
    max_lon = max(item[2] for item in bounds)
    max_lat = max(item[3] for item in bounds)
    return min_lon, min_lat, max_lon, max_lat


def _all_bounds(footprints: Sequence[Footprint]) -> tuple[float, float, float, float]:
    min_lon = min_lat = float("inf")
    max_lon = max_lat = float("-inf")
    for footprint in footprints:
        fmin_lon, fmin_lat, fmax_lon, fmax_lat = _footprint_bounds(footprint)
        min_lon = min(min_lon, fmin_lon)
        min_lat = min(min_lat, fmin_lat)
        max_lon = max(max_lon, fmax_lon)
        max_lat = max(max_lat, fmax_lat)
    return min_lon, min_lat, max_lon, max_lat


def _project(
    lon_deg: float,
    lat_deg: float,
    *,
    min_lon: float,
    max_lon: float,
    min_lat: float,
    max_lat: float,
    width: int,
    height: int,
) -> tuple[float, float]:
    lon_span = max(1.0e-9, max_lon - min_lon)
    lat_span = max(1.0e-9, max_lat - min_lat)
    x = ((lon_deg - min_lon) / lon_span) * width
    y = height - (((lat_deg - min_lat) / lat_span) * height)
    return x, y


def _ring_to_path(
    ring: tuple[tuple[float, float], ...],
    *,
    min_lon: float,
    max_lon: float,
    min_lat: float,
    max_lat: float,
    width: int,
    height: int,
) -> str:
    commands: list[str] = []
    for index, (lon, lat) in enumerate(ring):
        x, y = _project(
            lon,
            lat,
            min_lon=min_lon,
            max_lon=max_lon,
            min_lat=min_lat,
            max_lat=max_lat,
            width=width,
            height=height,
        )
        verb = "M" if index == 0 else "L"
        commands.append(f"{verb}{x:.2f},{y:.2f}")
    commands.append("Z")
    return " ".join(commands)


def _polygon_centroid(
    ring: tuple[tuple[float, float], ...],
) -> tuple[float, float]:
    if len(ring) < 3:
        return 0.0, 0.0
    area2 = 0.0
    cx = 0.0
    cy = 0.0
    for index, (x1, y1) in enumerate(ring):
        x2, y2 = ring[(index + 1) % len(ring)]
        cross = (x1 * y2) - (x2 * y1)
        area2 += cross
        cx += (x1 + x2) * cross
        cy += (y1 + y2) * cross
    if abs(area2) < 1.0e-12:
        lons = [lon for lon, _lat in ring]
        lats = [lat for _lon, lat in ring]
        return sum(lons) / len(lons), sum(lats) / len(lats)
    factor = 1.0 / (3.0 * area2)
    return cx * factor, cy * factor


def _footprint_label(footprint: Footprint) -> str:
    if footprint.has_explicit_height:
        return f"#{footprint.index} ele"
    return f"#{footprint.index} no-ele"


def build_svg(
    footprints: Sequence[Footprint],
    *,
    width: int,
    height: int,
    padding_ratio: float,
    background: str,
    show_labels: bool,
    show_legend: bool,
    ele_only: bool,
) -> str:
    if ele_only:
        footprints = tuple(footprint for footprint in footprints if footprint.has_explicit_height)
    if not footprints:
        raise ValueError("No footprints remain after applying filters")

    min_lon, min_lat, max_lon, max_lat = _all_bounds(footprints)
    lon_span = max_lon - min_lon
    lat_span = max_lat - min_lat
    pad_lon = lon_span * max(0.0, padding_ratio)
    pad_lat = lat_span * max(0.0, padding_ratio)
    min_lon -= pad_lon
    max_lon += pad_lon
    min_lat -= pad_lat
    max_lat += pad_lat

    body: list[str] = []
    body.append(
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="{html.escape(background)}"/>'
    )

    for fraction in (0.25, 0.5, 0.75):
        x = width * fraction
        y = height * fraction
        body.append(
            f'<line x1="{x:.2f}" y1="0" x2="{x:.2f}" y2="{height}" stroke="{GRID_STROKE}" stroke-width="1" opacity="0.35"/>'
        )
        body.append(
            f'<line x1="0" y1="{y:.2f}" x2="{width}" y2="{y:.2f}" stroke="{GRID_STROKE}" stroke-width="1" opacity="0.35"/>'
        )

    for footprint in footprints:
        stroke = ELE_STROKE if footprint.has_explicit_height else NO_ELE_STROKE
        for ring in footprint.outer_rings_lonlat:
            path_d = _ring_to_path(
                ring,
                min_lon=min_lon,
                max_lon=max_lon,
                min_lat=min_lat,
                max_lat=max_lat,
                width=width,
                height=height,
            )
            body.append(
                f'<path d="{path_d}" fill="none" stroke="{stroke}" stroke-width="2.5" vector-effect="non-scaling-stroke"/>'
            )

        if show_labels:
            label_ring = footprint.outer_rings_lonlat[0]
            cx_lon, cy_lat = _polygon_centroid(label_ring)
            cx, cy = _project(
                cx_lon,
                cy_lat,
                min_lon=min_lon,
                max_lon=max_lon,
                min_lat=min_lat,
                max_lat=max_lat,
                width=width,
                height=height,
            )
            body.append(
                f'<text x="{cx:.2f}" y="{cy:.2f}" fill="{TEXT_FILL}" font-size="22" font-family="monospace" text-anchor="middle" dominant-baseline="middle">{html.escape(_footprint_label(footprint))}</text>'
            )
            body.append(
                f'<text x="{cx:.2f}" y="{cy + 24:.2f}" fill="{TEXT_FILL}" font-size="12" font-family="monospace" text-anchor="middle" opacity="0.85">{html.escape(footprint.water_id)}</text>'
            )

    if show_legend:
        legend_x = 24
        legend_y = 28
        body.append(
            f'<rect x="{legend_x}" y="{legend_y - 18}" rx="8" ry="8" width="360" height="92" fill="#ffffff" opacity="0.78" stroke="#cbd5e1"/>'
        )
        body.append(
            f'<line x1="{legend_x + 14}" y1="{legend_y + 10}" x2="{legend_x + 38}" y2="{legend_y + 10}" stroke="{ELE_STROKE}" stroke-width="3" vector-effect="non-scaling-stroke"/>'
        )
        body.append(
            f'<text x="{legend_x + 52}" y="{legend_y + 16}" fill="{TEXT_FILL}" font-size="16" font-family="sans-serif">ele / water_level present</text>'
        )
        body.append(
            f'<line x1="{legend_x + 14}" y1="{legend_y + 42}" x2="{legend_x + 38}" y2="{legend_y + 42}" stroke="{NO_ELE_STROKE}" stroke-width="3" vector-effect="non-scaling-stroke"/>'
        )
        body.append(
            f'<text x="{legend_x + 52}" y="{legend_y + 48}" fill="{TEXT_FILL}" font-size="16" font-family="sans-serif">no explicit height</text>'
        )
        body.append(
            f'<text x="{legend_x + 14}" y="{legend_y + 76}" fill="{TEXT_FILL}" font-size="12" font-family="sans-serif" opacity="0.75">Only outer rings are drawn; inner rings are skipped.</text>'
        )

    return "\n".join(
        [
            '<?xml version="1.0" encoding="UTF-8"?>',
            (
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
                f'viewBox="0 0 {width} {height}">'
            ),
            *body,
            "</svg>",
        ]
    )


def main() -> int:
    args = parse_args()
    input_path = (args.input_cache or args.input).expanduser()
    output_path = args.output.expanduser()
    footprints = load_footprints(input_path)
    svg = build_svg(
        footprints,
        width=int(args.width),
        height=int(args.height),
        padding_ratio=float(args.padding_ratio),
        background=str(args.background),
        show_labels=bool(args.show_labels),
        show_legend=bool(args.show_legend),
        ele_only=bool(args.ele_only),
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(svg, encoding="utf-8")
    print(f"Wrote {output_path} from {input_path} ({len(footprints)} footprints)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
