#!/usr/bin/env python3
"""Render coastline inner rings as filled sea-surface shapes in SVG.

This script is a diagnostic companion to the outline renderer. It keeps only
coastline footprints, then draws their inner rings as filled polygons so the
result approximates the water surface rather than the bbox frame.
"""

from __future__ import annotations

import argparse
import html
import sys
from pathlib import Path
from typing import Sequence


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from render_water_footprints_svg_outline import (  # type: ignore[import-not-found]  # noqa: E402
    BACKGROUND,
    GRID_STROKE,
    TEXT_FILL,
    Footprint,
    _all_bounds,
    _footprint_label,
    _polygon_centroid,
    _project,
    _ring_to_path,
    load_footprints,
)


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = REPO_ROOT / "dev-samples" / "coastline_inner_rings_preview.svg"

SEA_FILL = "#4c8cff"
SEA_STROKE = "#2457b3"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render coastline inner rings from a zstarview cache JSON into SVG.",
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
        help="Draw a small legend explaining the fill color.",
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


def _inner_rings(footprint: Footprint) -> tuple[tuple[tuple[float, float], ...], ...]:
    return footprint.inner_rings_lonlat


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

    footprints = tuple(footprint for footprint in footprints if footprint.inner_rings_lonlat)
    if not footprints:
        raise ValueError("No coastline inner rings remain after applying filters")

    min_lon, min_lat, max_lon, max_lat = _all_bounds(
        tuple(
            Footprint(
                index=footprint.index,
                water_id=footprint.water_id,
                kind=footprint.kind,
                outer_rings_lonlat=footprint.inner_rings_lonlat,
                inner_rings_lonlat=(),
                source=footprint.source,
                tags=footprint.tags,
            )
            for footprint in footprints
        )
    )
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
        for ring in _inner_rings(footprint):
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
                f'<path d="{path_d}" fill="{SEA_FILL}" fill-opacity="0.55" stroke="{SEA_STROKE}" stroke-width="2.0" vector-effect="non-scaling-stroke"/>'
            )

        if show_labels:
            label_ring = footprint.inner_rings_lonlat[0]
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
            f'<rect x="{legend_x}" y="{legend_y - 18}" rx="8" ry="8" width="360" height="76" fill="#ffffff" opacity="0.78" stroke="#cbd5e1"/>'
        )
        body.append(
            f'<rect x="{legend_x + 14}" y="{legend_y - 2}" width="24" height="24" fill="{SEA_FILL}" fill-opacity="0.55" stroke="{SEA_STROKE}" stroke-width="2"/>'
        )
        body.append(
            f'<text x="{legend_x + 52}" y="{legend_y + 16}" fill="{TEXT_FILL}" font-size="16" font-family="sans-serif">coastline inner rings filled</text>'
        )
        body.append(
            f'<text x="{legend_x + 14}" y="{legend_y + 52}" fill="{TEXT_FILL}" font-size="12" font-family="sans-serif" opacity="0.75">Outer bbox ring is suppressed.</text>'
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
