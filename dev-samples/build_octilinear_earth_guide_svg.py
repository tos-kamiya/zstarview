#!/usr/bin/env python3
"""Generate a coarse world-map SVG with octilinear land outlines.

This script is for guide-layer prototyping, not for production rendering.
It accepts Polygon / MultiPolygon GeoJSON, simplifies the geometry with a
pure-Python Douglas-Peucker pass, and then deforms each edge so the final SVG
uses only horizontal, vertical, or 45 degree segments.

Typical source data:

- Natural Earth land polygons converted to GeoJSON
- Any other world landmass GeoJSON in EPSG:4326 lon/lat coordinates

The output stays 2D on purpose so the contour quality can be evaluated before
the runtime below-horizon projection is implemented.
"""

from __future__ import annotations

import argparse
import html
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = REPO_ROOT / "dev-samples" / "octilinear_earth_guide.svg"
DEFAULT_SIMPLIFY_DEG = 2.0
LARGE_RING_AREA_DEG2 = 200.0
LARGE_RING_MAX_SIMPLIFY_DEG = 1.8
SMALL_RING_AREA_DEG2 = 80.0
SMALL_RING_MAX_SIMPLIFY_DEG = 1.4
TINY_RING_AREA_DEG2 = 20.0
TINY_RING_MAX_SIMPLIFY_DEG = 0.8


@dataclass(frozen=True)
class Ring:
    points: tuple[tuple[float, float], ...]
    source_name: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render a simplified SVG preview from GeoJSON world land polygons. "
            "Input coordinates must be lon/lat degrees."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="GeoJSON input path containing Polygon or MultiPolygon features.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"SVG output path. Default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument("--width", type=int, default=1440, help="SVG width in px.")
    parser.add_argument("--height", type=int, default=720, help="SVG height in px.")
    parser.add_argument(
        "--seam-longitude-deg",
        type=float,
        default=-180.0,
        help=(
            "Longitude placed at the left/right map seam in degrees. "
            "Example: around -169 for the Bering Strait."
        ),
    )
    parser.add_argument(
        "--simplify-deg",
        type=float,
        default=DEFAULT_SIMPLIFY_DEG,
        help="Douglas-Peucker tolerance in lon/lat degrees.",
    )
    parser.add_argument(
        "--min-ring-area-deg2",
        type=float,
        default=12.0,
        help="Drop rings whose approximate lon/lat area is below this threshold.",
    )
    parser.add_argument(
        "--min-points",
        type=int,
        default=4,
        help="Minimum vertex count to keep after simplification.",
    )
    parser.add_argument(
        "--transparent-background",
        action="store_true",
        help="Leave the background transparent instead of drawing an ocean rectangle.",
    )
    parser.add_argument(
        "--label",
        action="store_true",
        help="Draw source feature labels at rough ring centroids.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        help="Optional runtime-oriented JSON output path for simplified rings.",
    )
    return parser.parse_args()


def wrap_longitude_relative_to_seam(lon_deg: float, seam_longitude_deg: float) -> float:
    return (lon_deg - seam_longitude_deg) % 360.0


def unwrap_ring_to_seam(
    points: tuple[tuple[float, float], ...],
    seam_longitude_deg: float,
) -> tuple[tuple[float, float], ...]:
    if not points:
        return ()
    wrapped = [wrap_longitude_relative_to_seam(lon_deg, seam_longitude_deg) for lon_deg, _lat_deg in points]
    unwrapped = [wrapped[0]]
    previous = wrapped[0]
    for wrapped_lon in wrapped[1:]:
        candidate = wrapped_lon
        while candidate - previous > 180.0:
            candidate -= 360.0
        while candidate - previous < -180.0:
            candidate += 360.0
        unwrapped.append(candidate)
        previous = candidate
    return tuple((lon_deg, lat_deg) for lon_deg, (_orig_lon, lat_deg) in zip(unwrapped, points))


def unwrap_ring_to_reference_lon(
    points: tuple[tuple[float, float], ...],
    reference_lon_deg: float,
) -> tuple[tuple[float, float], ...]:
    if not points:
        return ()
    wrapped = [wrap_longitude_relative_to_seam(lon_deg, reference_lon_deg) for lon_deg, _lat_deg in points]
    unwrapped = [wrapped[0]]
    previous = wrapped[0]
    for wrapped_lon in wrapped[1:]:
        candidate = wrapped_lon
        while candidate - previous > 180.0:
            candidate -= 360.0
        while candidate - previous < -180.0:
            candidate += 360.0
        unwrapped.append(candidate)
        previous = candidate
    return tuple((lon_deg, lat_deg) for lon_deg, (_orig_lon, lat_deg) in zip(unwrapped, points))


def lonlat_to_xy(
    lon_deg: float,
    lat_deg: float,
    width: int,
    height: int,
    *,
    seam_longitude_deg: float,
) -> tuple[float, float]:
    x = (wrap_longitude_relative_to_seam(lon_deg, seam_longitude_deg) / 360.0) * width
    y = ((90.0 - lat_deg) / 180.0) * height
    return x, y


def lonlat_to_xy_unwrapped(
    lon_deg: float,
    lat_deg: float,
    width: int,
    height: int,
    *,
    seam_longitude_deg: float,
) -> tuple[float, float]:
    x = ((lon_deg - seam_longitude_deg) / 360.0) * width
    y = ((90.0 - lat_deg) / 180.0) * height
    return x, y


def normalize_unwrapped_ring(
    points: tuple[tuple[float, float], ...],
    *,
    seam_longitude_deg: float,
    reference_lon_deg: float | None = None,
) -> tuple[tuple[float, float], ...]:
    if not points:
        return ()
    lons = [lon for lon, _lat in points]
    if reference_lon_deg is None:
        reference_lon_deg = (min(lons) + max(lons)) * 0.5
    center_lon = float(reference_lon_deg)
    shift_turns = math.floor((center_lon - seam_longitude_deg) / 360.0)
    shift_deg = shift_turns * 360.0
    return tuple((lon_deg - shift_deg, lat_deg) for lon_deg, lat_deg in points)


def polygon_area_deg2(points: tuple[tuple[float, float], ...]) -> float:
    if len(points) < 3:
        return 0.0
    total = 0.0
    for index, (x1, y1) in enumerate(points):
        x2, y2 = points[(index + 1) % len(points)]
        total += (x1 * y2) - (x2 * y1)
    return abs(total) * 0.5


def centroid_of_ring(points: tuple[tuple[float, float], ...]) -> tuple[float, float]:
    if not points:
        return 0.0, 0.0
    area2 = 0.0
    cx = 0.0
    cy = 0.0
    for index, (x1, y1) in enumerate(points):
        x2, y2 = points[(index + 1) % len(points)]
        cross = (x1 * y2) - (x2 * y1)
        area2 += cross
        cx += (x1 + x2) * cross
        cy += (y1 + y2) * cross
    if abs(area2) < 1.0e-12:
        mean_x = sum(x for x, _ in points) / len(points)
        mean_y = sum(y for _, y in points) / len(points)
        return mean_x, mean_y
    scale = 1.0 / (3.0 * area2)
    return cx * scale, cy * scale


def point_segment_distance(point: tuple[float, float], start: tuple[float, float], end: tuple[float, float]) -> float:
    px, py = point
    x1, y1 = start
    x2, y2 = end
    dx = x2 - x1
    dy = y2 - y1
    if dx == 0.0 and dy == 0.0:
        return math.hypot(px - x1, py - y1)
    t = ((px - x1) * dx + (py - y1) * dy) / ((dx * dx) + (dy * dy))
    t = max(0.0, min(1.0, t))
    proj_x = x1 + (t * dx)
    proj_y = y1 + (t * dy)
    return math.hypot(px - proj_x, py - proj_y)


def douglas_peucker(points: tuple[tuple[float, float], ...], tolerance: float) -> tuple[tuple[float, float], ...]:
    if len(points) <= 2:
        return points
    start = points[0]
    end = points[-1]
    max_distance = -1.0
    split_index = -1
    for index in range(1, len(points) - 1):
        dist = point_segment_distance(points[index], start, end)
        if dist > max_distance:
            max_distance = dist
            split_index = index
    if max_distance <= tolerance or split_index < 0:
        return (start, end)
    left = douglas_peucker(points[: split_index + 1], tolerance)
    right = douglas_peucker(points[split_index:], tolerance)
    return left[:-1] + right


def close_ring(points: Iterable[tuple[float, float]]) -> tuple[tuple[float, float], ...]:
    pts = tuple((float(lon), float(lat)) for lon, lat in points)
    if len(pts) >= 2 and pts[0] == pts[-1]:
        pts = pts[:-1]
    return pts


def simplify_ring(
    points: tuple[tuple[float, float], ...],
    *,
    tolerance: float,
    min_points: int,
) -> tuple[tuple[float, float], ...]:
    pts = close_ring(points)
    if len(pts) < min_points:
        return ()
    open_pts = pts + (pts[0],)
    simplified = douglas_peucker(open_pts, tolerance)
    simplified = close_ring(simplified)
    if len(simplified) < min_points:
        return ()
    return simplified


def adaptive_simplify_tolerance(area_deg2: float, base_tolerance: float) -> float:
    if area_deg2 < TINY_RING_AREA_DEG2:
        return min(base_tolerance, TINY_RING_MAX_SIMPLIFY_DEG)
    if area_deg2 < SMALL_RING_AREA_DEG2:
        return min(base_tolerance, SMALL_RING_MAX_SIMPLIFY_DEG)
    if area_deg2 < LARGE_RING_AREA_DEG2:
        return min(base_tolerance, LARGE_RING_MAX_SIMPLIFY_DEG)
    return base_tolerance


def adaptive_min_ring_area(area_deg2: float, base_min_area_deg2: float) -> float:
    if area_deg2 < TINY_RING_AREA_DEG2:
        return min(base_min_area_deg2, 2.0)
    if area_deg2 < SMALL_RING_AREA_DEG2:
        return min(base_min_area_deg2, 4.0)
    if area_deg2 < LARGE_RING_AREA_DEG2:
        return min(base_min_area_deg2, 8.0)
    return base_min_area_deg2


def _iter_lonlat_pairs(raw_points: list[Any]) -> Iterable[tuple[float, float]]:
    for item in raw_points:
        if not isinstance(item, list) or len(item) < 2:
            continue
        yield float(item[0]), float(item[1])


def feature_name(feature: dict[str, Any], default_index: int) -> str:
    props = feature.get("properties")
    if isinstance(props, dict):
        for key in ("name", "NAME", "name_en", "continent", "CONTINENT", "featurecla"):
            value = props.get(key)
            if isinstance(value, str) and value.strip():
                return value.strip()
    return f"feature-{default_index}"


def should_render_label(name: str) -> bool:
    lowered = name.strip().casefold()
    return lowered not in {"", "land", "geometry"}


def rings_from_geometry(geometry: dict[str, Any], source_name: str) -> list[Ring]:
    geom_type = geometry.get("type")
    coords = geometry.get("coordinates")
    rings: list[Ring] = []
    if geom_type == "Polygon" and isinstance(coords, list):
        if coords:
            outer = coords[0]
            if isinstance(outer, list):
                rings.append(Ring(points=close_ring(_iter_lonlat_pairs(outer)), source_name=source_name))
    elif geom_type == "MultiPolygon" and isinstance(coords, list):
        for polygon in coords:
            if isinstance(polygon, list) and polygon:
                outer = polygon[0]
                if isinstance(outer, list):
                    rings.append(Ring(points=close_ring(_iter_lonlat_pairs(outer)), source_name=source_name))
    else:
        raise ValueError(f"Unsupported geometry type: {geom_type!r}")
    return rings


def load_geojson_rings(path: Path) -> list[Ring]:
    text = path.read_text(encoding="utf-8")
    stripped = text.lstrip()
    if stripped.startswith("<!DOCTYPE html") or stripped.startswith("<html"):
        raise ValueError(
            "Input file looks like HTML, not GeoJSON. "
            "You likely saved a GitHub page instead of the raw GeoJSON payload."
        )
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        snippet = stripped[:80].replace("\n", "\\n")
        raise ValueError(
            f"Input is not valid JSON near: {snippet!r}"
        ) from exc
    if not isinstance(payload, dict):
        raise ValueError("GeoJSON payload must be a dict.")

    payload_type = payload.get("type")
    rings: list[Ring] = []
    if payload_type == "FeatureCollection":
        features = payload.get("features")
        if not isinstance(features, list):
            raise ValueError("FeatureCollection.features must be a list.")
        for index, feature in enumerate(features):
            if not isinstance(feature, dict):
                continue
            geometry = feature.get("geometry")
            if not isinstance(geometry, dict):
                continue
            rings.extend(rings_from_geometry(geometry, feature_name(feature, index)))
    elif payload_type == "Feature":
        geometry = payload.get("geometry")
        if not isinstance(geometry, dict):
            raise ValueError("Feature.geometry must be a dict.")
        rings.extend(rings_from_geometry(geometry, feature_name(payload, 0)))
    elif payload_type in {"Polygon", "MultiPolygon"}:
        rings.extend(rings_from_geometry(payload, "geometry"))
    else:
        raise ValueError(f"Unsupported GeoJSON top-level type: {payload_type!r}")
    return rings


def octilinearize_segment(
    start: tuple[float, float],
    end: tuple[float, float],
) -> list[tuple[float, float]]:
    x0, y0 = start
    x1, y1 = end
    dx = x1 - x0
    dy = y1 - y0
    adx = abs(dx)
    ady = abs(dy)
    eps = 1.0e-9
    if adx <= eps and ady <= eps:
        return []
    if adx <= eps or ady <= eps or abs(adx - ady) <= eps:
        return [(x1, y1)]

    sx = 1.0 if dx >= 0.0 else -1.0
    sy = 1.0 if dy >= 0.0 else -1.0
    if adx > ady:
        corner = (x0 + (sx * ady), y0 + (sy * ady))
    else:
        corner = (x0 + (sx * adx), y0 + (sy * adx))
    if abs(corner[0] - x1) <= eps and abs(corner[1] - y1) <= eps:
        return [(x1, y1)]
    return [corner, (x1, y1)]


def _clip_polygon_against_edge(
    points: list[tuple[float, float]],
    *,
    inside: Callable[[tuple[float, float]], bool],
    intersect: Callable[[tuple[float, float], tuple[float, float]], tuple[float, float]],
) -> list[tuple[float, float]]:
    if not points:
        return []
    clipped: list[tuple[float, float]] = []
    previous = points[-1]
    previous_inside = inside(previous)
    for current in points:
        current_inside = inside(current)
        if current_inside:
            if not previous_inside:
                clipped.append(intersect(previous, current))
            clipped.append(current)
        elif previous_inside:
            clipped.append(intersect(previous, current))
        previous = current
        previous_inside = current_inside
    return clipped


def clip_polygon_to_viewport(
    points: list[tuple[float, float]],
    width: int,
    height: int,
) -> list[tuple[float, float]]:
    clipped = points
    if not clipped:
        return []

    def clip_left(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
        ax, ay = a
        bx, by = b
        t = (0.0 - ax) / ((bx - ax) or 1.0e-12)
        return 0.0, ay + ((by - ay) * t)

    def clip_right(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
        ax, ay = a
        bx, by = b
        t = (float(width) - ax) / ((bx - ax) or 1.0e-12)
        return float(width), ay + ((by - ay) * t)

    def clip_top(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
        ax, ay = a
        bx, by = b
        t = (0.0 - ay) / ((by - ay) or 1.0e-12)
        return ax + ((bx - ax) * t), 0.0

    def clip_bottom(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
        ax, ay = a
        bx, by = b
        t = (float(height) - ay) / ((by - ay) or 1.0e-12)
        return ax + ((bx - ax) * t), float(height)

    clipped = _clip_polygon_against_edge(
        clipped,
        inside=lambda point: point[0] >= 0.0,
        intersect=clip_left,
    )
    clipped = _clip_polygon_against_edge(
        clipped,
        inside=lambda point: point[0] <= float(width),
        intersect=clip_right,
    )
    clipped = _clip_polygon_against_edge(
        clipped,
        inside=lambda point: point[1] >= 0.0,
        intersect=clip_top,
    )
    clipped = _clip_polygon_against_edge(
        clipped,
        inside=lambda point: point[1] <= float(height),
        intersect=clip_bottom,
    )
    return clipped


def reorder_ring_to_southern_endpoints(
    points: list[tuple[float, float]],
) -> tuple[list[tuple[float, float]], float]:
    if len(points) < 2:
        return points[:], 0.0

    southern = sorted(
        range(len(points)),
        key=lambda index: (points[index][1], points[index][0]),
    )[:2]
    first_index, second_index = southern[0], southern[1]

    def build_arc(step: int) -> list[tuple[float, float]]:
        arc = [points[first_index]]
        index = first_index
        while index != second_index:
            index = (index + step) % len(points)
            arc.append(points[index])
        return arc

    forward = build_arc(1)
    backward = build_arc(-1)

    def arc_score(arc: list[tuple[float, float]]) -> tuple[float, float]:
        lats = [lat for _lon, lat in arc]
        return max(lats), sum(lats) / len(lats)

    chosen = forward if arc_score(forward) >= arc_score(backward) else backward
    reference_lon = (points[first_index][0] + points[second_index][0]) * 0.5
    return chosen, reference_lon


def polygon_paths_d(
    ring: tuple[tuple[float, float], ...],
    width: int,
    height: int,
    *,
    seam_longitude_deg: float,
) -> list[str]:
    unwrapped = unwrap_ring_to_seam(ring, seam_longitude_deg)
    if not unwrapped:
        return []
    lat_min = min(lat for _lon, lat in unwrapped)
    ordered = list(unwrapped)
    south_reference_lon: float | None = None
    if lat_min < -60.0:
        ordered, south_reference_lon = reorder_ring_to_southern_endpoints(list(unwrapped))
        ordered = list(
            unwrap_ring_to_reference_lon(tuple(ordered), south_reference_lon)
        )
    normalized = normalize_unwrapped_ring(
        tuple(ordered),
        seam_longitude_deg=seam_longitude_deg,
        reference_lon_deg=south_reference_lon if south_reference_lon is not None else None,
    )
    subrings: list[list[tuple[float, float]]] = [list(normalized)]

    path_strings: list[str] = []
    for subring in subrings:
        if len(subring) < 3:
            continue
        start_x, start_y = lonlat_to_xy_unwrapped(
            subring[0][0],
            subring[0][1],
            width,
            height,
            seam_longitude_deg=seam_longitude_deg,
        )
        raw_points: list[tuple[float, float]] = [(start_x, start_y)]
        last_point = (start_x, start_y)
        projected = [
            lonlat_to_xy_unwrapped(
                lon_deg,
                lat_deg,
                width,
                height,
                seam_longitude_deg=seam_longitude_deg,
            )
            for lon_deg, lat_deg in subring[1:]
        ]
        for point in projected:
            for next_point in octilinearize_segment(last_point, point):
                raw_points.append(next_point)
                last_point = next_point
        for next_point in octilinearize_segment(last_point, (start_x, start_y)):
            raw_points.append(next_point)
            last_point = next_point
        clipped = clip_polygon_to_viewport(raw_points, width, height)
        if len(clipped) < 3:
            continue
        clipped_commands = [f"M {clipped[0][0]:.2f} {clipped[0][1]:.2f}"]
        for x, y in clipped[1:]:
            clipped_commands.append(f"L {x:.2f} {y:.2f}")
        clipped_commands.append("Z")
        path_strings.append(" ".join(clipped_commands))
    return path_strings


def build_svg(
    *,
    rings: list[Ring],
    width: int,
    height: int,
    seam_longitude_deg: float,
    transparent_background: bool,
    label: bool,
) -> str:
    lines: list[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
            f'style="overflow:hidden">'
        ),
        "<title>Octilinear Earth Guide Preview</title>",
        (
            "<desc>"
            "Simplified land polygons with horizontal, vertical, or 45 degree edges."
            "</desc>"
        ),
    ]
    if not transparent_background:
        lines.append(
            f'<rect x="0" y="0" width="{width}" height="{height}" fill="#10263d" />'
        )

    lines.append('<g id="graticule" stroke="#4f6f89" stroke-width="1" opacity="0.55">')
    for offset_deg in range(30, 361, 30):
        lon_deg = seam_longitude_deg + float(offset_deg)
        x, _ = lonlat_to_xy(
            lon_deg,
            0.0,
            width,
            height,
            seam_longitude_deg=seam_longitude_deg,
        )
        lines.append(f'<line x1="{x:.2f}" y1="0" x2="{x:.2f}" y2="{height}" />')
    for lat_deg in range(-60, 61, 30):
        _, y = lonlat_to_xy(
            seam_longitude_deg,
            float(lat_deg),
            width,
            height,
            seam_longitude_deg=seam_longitude_deg,
        )
        lines.append(f'<line x1="0" y1="{y:.2f}" x2="{width}" y2="{y:.2f}" />')
    lines.append("</g>")

    lines.append(
        '<g id="land" fill="none" stroke="#d7ddd5" '
        'stroke-width="1.2" stroke-linejoin="round" stroke-linecap="round">'
    )
    for ring in rings:
        for path_d in polygon_paths_d(
            ring.points,
            width,
            height,
            seam_longitude_deg=seam_longitude_deg,
        ):
            if not path_d:
                continue
            lines.append(
                f'<path d="{path_d}" data-name="{html.escape(ring.source_name)}" />'
            )
    lines.append("</g>")

    if label:
        lines.append(
            '<g id="labels" fill="#dce9f5" font-family="monospace" font-size="16" opacity="0.85">'
        )
        for ring in rings:
            if not should_render_label(ring.source_name):
                continue
            lon_deg, lat_deg = centroid_of_ring(ring.points)
            x, y = lonlat_to_xy(
                lon_deg,
                lat_deg,
                width,
                height,
                seam_longitude_deg=seam_longitude_deg,
            )
            escaped = html.escape(ring.source_name)
            lines.append(f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="middle">{escaped}</text>')
        lines.append("</g>")

    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def ring_bbox(points: tuple[tuple[float, float], ...]) -> dict[str, float]:
    lons = [lon for lon, _ in points]
    lats = [lat for _, lat in points]
    return {
        "min_lon_deg": min(lons),
        "max_lon_deg": max(lons),
        "min_lat_deg": min(lats),
        "max_lat_deg": max(lats),
    }


def build_runtime_payload(
    *,
    source_path: Path,
    rings: list[Ring],
    simplify_deg: float,
    min_ring_area_deg2: float,
    min_points: int,
    seam_longitude_deg: float,
) -> dict[str, Any]:
    ring_items = []
    for ring in rings:
        area_deg2 = polygon_area_deg2(ring.points)
        ring_items.append(
            {
                "source_name": ring.source_name,
                "label_name": ring.source_name if should_render_label(ring.source_name) else None,
                "point_count": len(ring.points),
                "approx_area_deg2": round(area_deg2, 6),
                "bbox": ring_bbox(ring.points),
                "points_lonlat_deg": [[round(lon, 6), round(lat, 6)] for lon, lat in ring.points],
            }
        )
    ring_items.sort(
        key=lambda item: (-float(item["approx_area_deg2"]), str(item["source_name"]))
    )
    return {
        "format": "zstarview-earth-guide-v1",
        "source": {
            "path": str(source_path),
            "coordinate_order": "lon_lat_deg",
        },
        "simplification": {
            "algorithm": "douglas_peucker",
            "tolerance_deg": simplify_deg,
            "min_ring_area_deg2": min_ring_area_deg2,
            "min_points": min_points,
        },
        "preview_projection": {
            "type": "equirectangular",
            "seam_longitude_deg": seam_longitude_deg,
        },
        "rendering": {
            "line_mode": "octilinear_45deg",
        },
        "ring_count": len(ring_items),
        "rings": ring_items,
    }


def main() -> None:
    args = parse_args()
    if args.width <= 0 or args.height <= 0:
        raise SystemExit("width and height must be positive")
    if args.simplify_deg < 0.0:
        raise SystemExit("simplify-deg must be non-negative")
    if args.min_points < 3:
        raise SystemExit("min-points must be at least 3")

    raw_rings = load_geojson_rings(args.input)
    simplified_rings: list[Ring] = []
    dropped_small = 0
    dropped_short = 0
    for ring in raw_rings:
        raw_area = polygon_area_deg2(ring.points)
        tolerance = adaptive_simplify_tolerance(
            raw_area,
            float(args.simplify_deg),
        )
        simplified = simplify_ring(
            ring.points,
            tolerance=tolerance,
            min_points=int(args.min_points),
        )
        if not simplified:
            dropped_short += 1
            continue
        min_ring_area = adaptive_min_ring_area(raw_area, float(args.min_ring_area_deg2))
        if raw_area < min_ring_area:
            dropped_small += 1
            continue
        simplified_rings.append(Ring(points=simplified, source_name=ring.source_name))

    if not simplified_rings:
        raise SystemExit("No rings left after simplification/filtering.")

    svg = build_svg(
        rings=simplified_rings,
        width=int(args.width),
        height=int(args.height),
        seam_longitude_deg=float(args.seam_longitude_deg),
        transparent_background=bool(args.transparent_background),
        label=bool(args.label),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(svg, encoding="utf-8")
    if args.output_json is not None:
        payload = build_runtime_payload(
            source_path=args.input,
            rings=simplified_rings,
            simplify_deg=float(args.simplify_deg),
            min_ring_area_deg2=float(args.min_ring_area_deg2),
            min_points=int(args.min_points),
            seam_longitude_deg=float(args.seam_longitude_deg),
        )
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
    print(
        "Wrote octilinear Earth guide SVG to "
        f"{args.output} from {args.input} "
        f"({len(raw_rings)} raw rings -> {len(simplified_rings)} simplified; "
        f"dropped_short={dropped_short}, dropped_small={dropped_small})"
    )
    if args.output_json is not None:
        print(f"Wrote runtime Earth guide JSON to {args.output_json}")


if __name__ == "__main__":
    main()
