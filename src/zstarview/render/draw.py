import math
import re
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Tuple

import hashlib
import logging
import numpy as np
from zoneinfo import ZoneInfo

from PySide6.QtCore import QPoint, QPointF, QRectF, Qt
from PySide6.QtGui import (
    QFont,
    QFontMetrics,
    QImage,
    QColor,
    QPainter,
    QPen,
    QPolygonF,
    QRadialGradient,
)

from ..paths import (
    FIELD_OF_VIEW_DEG,
    BACKGROUND_FIELD_OF_VIEW_DEG1,
    BACKGROUND_FIELD_OF_VIEW_DEG2,
    TEXT_STYLES_BY_PRESET,
    CELESTIAL_EQUATOR_COLOR,
    DIRECTIONS,
    ECLIPTIC_COLOR,
    HORIZON_LINE_COLOR,
    TERRAIN_HORIZON_LINE_COLOR,
    URBAN_OUTLINE_LAYER_LINE_COLOR,
)
from ..aircraft.types import AircraftOverlayPoint
from ..aircraft_constants import (
    AIRCRAFT_FADE_START_SECONDS,
    AIRCRAFT_OVERLAY_LINE_COLOR_RGB,
)
from ..satellite_constants import (
    SATELLITE_OVERLAY_MARKER_COLOR_RGB,
    SATELLITE_OVERLAY_MARKER_MAX_ALPHA,
)
from ..satellites.types import SatelliteOverlayPoint
from ..types import ScreenGeometry, CelestialData, ViewerData, CelestialObject, PlanetBody, UrbanOutlinePolyline
from ..astro import (
    altaz_to_normalized_xy,
    lookup_star_name,
    lookup_star_source_id,
    resolve_star_names,
    resolve_star_source_ids,
    is_in_fov,
    is_in_fov_vectorized,
    calculate_moon_render_data,
)
from ..asterisms import ASTERISMS, ASTERISM_REQUIRED_SOURCE_IDS, pick_rotating_asterism
from . import geometry as render_geometry
from .geometry import (
    _altaz_to_normalized_xy_vectorized,
    _normalized_to_screen_xy_vectorized,
    normalized_to_screen_xy,
)
from .photometry import (
    compute_flare_profile as _compute_flare_profile,
    flare_strength_from_vmag as _flare_strength_from_vmag,
    body_label_text,
    bv_to_rgb_vectorized,
    planet_bloom_profile_from_vmag,
    planet_disc_style_from_vmag,
    planet_marker_color,
)
from . import text as render_text
from .text import (
    get_text_style,
    _rect_overlap_count,
    _text_bounds_at_baseline,
    draw_outlined_text,
)
from ..utils.image import generate_moon_phase_rgba
from .qt_image import np_rgba_to_qimage

logger = logging.getLogger(__name__)

# Backward-compatible re-export for callers importing geometry helpers via render.draw.
get_screen_geometry = render_geometry.get_screen_geometry
get_text_style = render_text.get_text_style
draw_label_candidates = render_text._draw_label_candidates
draw_status_line_text = render_text._draw_status_line_text
# Keep these re-exports because focused render tests import photometry helpers via render.draw.
compute_flare_profile = _compute_flare_profile
flare_strength_from_vmag = _flare_strength_from_vmag

_star_render_cache: tuple[tuple, QImage] | None = None
_MAG2_TO_MAG1_SIZE_SCALE = 10.0 ** 0.12
_DIAMOND_OVERLAY_GAIN = 0.85
_SINGLE_STAR_GAUSSIAN_STRENGTH = 0.12
_DSO_SHAPE_SIZE_GAIN = 1.0
_DSO_HOVER_SIZE_GAIN = 3.0
_DSO_SHAPE_MIN_MAJOR_ARCMIN = 15.0
_DSO_CATALOG_LIKE_NAME_RE = re.compile(r"^(M\d+|NGC\d+|IC\d+|MEL\d+|MWSC\d+)$", re.IGNORECASE)


def format_overlay_info_lines(
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    vmag_limit: float,
) -> list[str]:
    """Return the static top-left overlay text lines."""

    def format_height_m(value_m: float) -> str:
        rounded = round(float(value_m))
        if abs(float(value_m) - rounded) < 0.05:
            return f"{int(rounded)} m"
        return f"{float(value_m):.1f} m"

    lines = [viewer_data.city_name]
    if viewer_data.location_height_label and viewer_data.location_height_m is not None:
        lines.append(f"{viewer_data.location_height_label} {format_height_m(viewer_data.location_height_m)}")
    if viewer_data.show_observer_height:
        lines.append(f"Observer height {format_height_m(viewer_data.observer_height_m)}")

    utc_time = celestial_data.time
    tz_name = viewer_data.timezone_name
    try:
        local_tz = ZoneInfo(tz_name)
        local_dt = utc_time.to_datetime(timezone=local_tz)
        lines.append(local_dt.strftime("%Y-%m-%d %H:%M:%S %Z"))
    except Exception:
        lines.append(utc_time.to_datetime().strftime("%Y-%m-%d %H:%M:%S UTC"))

    alt_deg, az_deg = viewer_data.view_center

    def az_to_compass(az: float) -> str:
        """Converts azimuth in degrees to a compass direction string."""
        names = tuple(DIRECTIONS.keys())
        sector = 360.0 / len(names)
        idx = int(((az % 360.0) + sector / 2.0) // sector) % len(names)
        return names[idx]

    compass = az_to_compass(az_deg)
    deg = "\N{DEGREE SIGN}"
    lines.append(f"Alt {alt_deg:.0f}{deg}  Az {az_deg:.0f}{deg} ({compass})")
    lines.append(f"Vmag limit {vmag_limit:.1f}")
    return lines


def _content_fov_deg_from_viewer(viewer_data: ViewerData) -> float:
    return float(viewer_data.content_fov_deg)


def _is_named_dso(name: object, obj_id: object) -> bool:
    name_s = "" if name is None else str(name).strip()
    if not name_s:
        return False
    if _DSO_CATALOG_LIKE_NAME_RE.match(name_s):
        return False
    id_s = "" if obj_id is None else str(obj_id).strip()
    return name_s.upper().replace(" ", "") != id_s.upper().replace(" ", "")


def _dso_ellipse_polygon(
    *,
    alt_deg: float,
    az_deg: float,
    major_arcmin: float,
    minor_arcmin: float,
    pa_deg: float,
    viewer_data: ViewerData,
    geometry: ScreenGeometry,
    gain: float,
    samples: int = 48,
) -> QPolygonF:
    """Sample an ellipse on the local sky frame, then project each point to screen."""
    a_deg = 0.5 * (float(major_arcmin) / 60.0) * float(gain)
    b_deg = 0.5 * (float(minor_arcmin) / 60.0) * float(gain)
    pa_rad = math.radians(float(pa_deg) if math.isfinite(pa_deg) else 0.0)
    c = math.cos(pa_rad)
    s = math.sin(pa_rad)
    cos_alt = max(0.15, math.cos(math.radians(float(alt_deg))))

    pts: List[QPointF] = []
    n = max(16, int(samples))
    for i in range(n):
        t = 2.0 * math.pi * (i / n)
        # Local north/east offsets [deg] on tangent plane around center.
        local_major = a_deg * math.cos(t)
        local_minor = b_deg * math.sin(t)
        d_north = (local_major * c) + (local_minor * -s)
        d_east = (local_major * s) + (local_minor * c)

        alt_i = float(alt_deg) + d_north
        alt_i = max(-89.8, min(89.8, alt_i))
        az_i = (float(az_deg) + (d_east / cos_alt)) % 360.0
        nx, ny = altaz_to_normalized_xy(alt_i, az_i, viewer_data.view_center)
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        pts.append(QPointF(px, py))
    return QPolygonF(pts)


def _array_hash(arr: np.ndarray) -> str:
    if arr.size == 0:
        return "empty"
    return hashlib.md5(arr.tobytes()).hexdigest()


def _apply_weak_gaussian3x3_rgb(arr: np.ndarray, strength: float) -> np.ndarray:
    """Apply a very weak 3x3 Gaussian only to reduce salt-like 1px star noise."""
    s = float(np.clip(strength, 0.0, 1.0))
    if s <= 0.0:
        return arr
    p = np.pad(arr, ((1, 1), (1, 1), (0, 0)), mode="constant")
    g = (
        p[:-2, :-2, :]
        + 2.0 * p[:-2, 1:-1, :]
        + p[:-2, 2:, :]
        + 2.0 * p[1:-1, :-2, :]
        + 4.0 * p[1:-1, 1:-1, :]
        + 2.0 * p[1:-1, 2:, :]
        + p[2:, :-2, :]
        + 2.0 * p[2:, 1:-1, :]
        + p[2:, 2:, :]
    ) / 16.0
    return arr * (1.0 - s) + g * s


def _star_cache_key(
    alt: np.ndarray,
    az: np.ndarray,
    size_factor: np.ndarray,
    color_factor_base: np.ndarray,
    celestial_time_value: float,
    view_center: Tuple[float, float],
    geometry: ScreenGeometry,
    star_base_radius: float,
    visibility_boost: float,
    draw_vmag_limit: float | None,
) -> tuple:
    return (
        _array_hash(alt),
        _array_hash(az),
        _array_hash(size_factor),
        _array_hash(color_factor_base),
        float(celestial_time_value),
        view_center[0],
        view_center[1],
        geometry.center,
        geometry.radius,
        float(star_base_radius),
        float(visibility_boost),
        None if draw_vmag_limit is None else float(draw_vmag_limit),
    )


def find_highlighted_object(
    celestial_data: Optional[CelestialData],
    viewer_data: ViewerData,
    mouse_pos: QPoint,
    geometry: ScreenGeometry,
) -> Optional[Tuple[CelestialObject, QPointF]]:
    """
    Find the nearest celestial object to the mouse cursor.

    This function searches for the closest star or planet to the given mouse
    position. It performs a vectorized search for stars and a scalar search for
    planets to optimize performance.

    Args:
        celestial_data: A CelestialData object containing information about celestial objects.
        viewer_data: A ViewerData object containing the viewer's location and time.
        mouse_pos: The QPoint representing the current mouse cursor position.
        geometry: A ScreenGeometry object for coordinate transformations.

    Returns:
        A tuple containing the highlighted celestial object (as a dictionary or
        a CelestialObject) and its screen position (QPointF), or None if no
        object is close enough to the cursor.
    """
    min_dist_sq = 30**2  # squared pixels
    highlighted_object: Optional[Tuple[CelestialObject, QPointF]] = None

    def _has_named_star(name: object) -> bool:
        if name is None:
            return False
        return bool(str(name).strip())

    def _is_asterism_member(source_id: object) -> bool:
        if source_id is None:
            return False
        return str(source_id).strip() in ASTERISM_REQUIRED_SOURCE_IDS

    if not celestial_data:
        return None

    # Handle stars first (vectorized)
    stars = celestial_data.stars
    if stars["alt"].size > 0:
        names = resolve_star_names(stars, celestial_data.star_catalog_meta)
        source_ids = resolve_star_source_ids(stars, celestial_data.star_catalog_meta)
        hover_mask = np.array(
            [_has_named_star(name) or _is_asterism_member(source_id) for name, source_id in zip(names, source_ids)],
            dtype=bool,
        )
        if np.any(hover_mask):
            valid_indices = np.nonzero(hover_mask)[0]
            alt_named = stars["alt"][hover_mask]
            az_named = stars["az"][hover_mask]
            nx, ny = _altaz_to_normalized_xy_vectorized(alt_named, az_named, viewer_data.view_center)
            x, y = _normalized_to_screen_xy_vectorized(nx, ny, geometry)
            dist_sq = (mouse_pos.x() - x) ** 2 + (mouse_pos.y() - y) ** 2
            closest_idx = np.argmin(dist_sq)
            if dist_sq[closest_idx] < min_dist_sq:
                min_dist_sq = dist_sq[closest_idx]
                original_idx = valid_indices[closest_idx]
                highlighted_star: CelestialObject = {key: val[original_idx] for key, val in stars.items()}
                highlighted_star["name"] = names[original_idx]
                highlighted_star["source_id"] = source_ids[original_idx]
                highlighted_object = (highlighted_star, QPointF(x[closest_idx], y[closest_idx]))

    content_fov_deg = _content_fov_deg_from_viewer(viewer_data)

    # Handle planets (scalar)
    for body in celestial_data.planets:
        # For planets, allow below-horizon display as long as they are in the current FOV.
        if not is_in_fov(body.alt, body.az, viewer_data.view_center, fov_deg=content_fov_deg):
            continue
        nx, ny = altaz_to_normalized_xy(body.alt, body.az, viewer_data.view_center)
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        dist_sq = (mouse_pos.x() - px) ** 2 + (mouse_pos.y() - py) ** 2
        if dist_sq < min_dist_sq:
            min_dist_sq = dist_sq
            highlighted_object = (body, QPointF(px, py))  # body is already an object/dataclass

    return highlighted_object


def find_highlighted_dso(
    celestial_data: Optional[CelestialData],
    viewer_data: ViewerData,
    mouse_pos: QPoint,
    geometry: ScreenGeometry,
) -> Optional[Tuple[CelestialObject, QPointF]]:
    """Find highlighted named DSO independently from star/planet highlight."""
    if not celestial_data:
        return None

    dso = celestial_data.deep_sky_objects
    if dso["alt"].size == 0:
        return None

    names = np.asarray(dso["name"], dtype=object)
    ids = np.asarray(dso["id"], dtype=object)
    major = np.asarray(dso["major_arcmin"], dtype=float)
    minor = np.asarray(dso["minor_arcmin"], dtype=float)
    pa = np.asarray(dso["pa_deg"], dtype=float)
    named_mask = np.array([_is_named_dso(n, i) for n, i in zip(names, ids)], dtype=bool)
    shape_mask = np.isfinite(major) & np.isfinite(minor) & (major >= _DSO_SHAPE_MIN_MAJOR_ARCMIN) & (minor > 0.0)
    valid = np.nonzero(named_mask & shape_mask)[0]
    if valid.size == 0:
        return None

    mouse_pf = QPointF(float(mouse_pos.x()), float(mouse_pos.y()))
    best: Optional[Tuple[CelestialObject, QPointF]] = None
    best_dist = float("inf")
    near_best: Optional[Tuple[CelestialObject, QPointF]] = None
    near_best_dist = 30.0**2
    for idx in valid:
        alt = float(dso["alt"][idx])
        az = float(dso["az"][idx])
        nx, ny = altaz_to_normalized_xy(alt, az, viewer_data.view_center)
        x, y = normalized_to_screen_xy(nx, ny, geometry)
        dist_sq = (mouse_pf.x() - x) ** 2 + (mouse_pf.y() - y) ** 2
        poly = _dso_ellipse_polygon(
            alt_deg=alt,
            az_deg=az,
            major_arcmin=float(major[idx]),
            minor_arcmin=float(minor[idx]),
            pa_deg=float(pa[idx]) if np.isfinite(pa[idx]) else 0.0,
            viewer_data=viewer_data,
            geometry=geometry,
            gain=_DSO_SHAPE_SIZE_GAIN,
            samples=48,
        )
        candidate: CelestialObject = {
            "name": dso["name"][idx],
            "alt": alt,
            "az": az,
            "id": dso["id"][idx],
            "type": dso["type"][idx],
            "major_arcmin": major[idx],
            "minor_arcmin": minor[idx],
            "pa_deg": pa[idx],
            "_kind": "dso",
        }
        if poly.containsPoint(mouse_pf, Qt.FillRule.OddEvenFill) and dist_sq < best_dist:
            best_dist = dist_sq
            best = (candidate, QPointF(x, y))
        if dist_sq < near_best_dist:
            near_best_dist = dist_sq
            near_best = (candidate, QPointF(x, y))

    return best if best is not None else near_best


def draw_radial_background(
    painter: QPainter,
    rect: QRectF,
    geometry: ScreenGeometry,
    *,
    preset: str = "night",
    content_fov_deg: float = BACKGROUND_FIELD_OF_VIEW_DEG2,
) -> None:
    """
    Draws a radial gradient background to represent the sky.

    This function creates a subtle radial gradient to simulate the sky's
    appearance, with a darker tone towards the zenith and a slightly
    lighter tone towards the horizon.

    Args:
        painter: The QPainter object to use for drawing.
        rect: The QRectF defining the area to fill with the background.
        geometry: A ScreenGeometry object for calculating gradient parameters.
    """
    assert geometry.radius >= 10
    fov_outer = max(float(BACKGROUND_FIELD_OF_VIEW_DEG1), float(content_fov_deg))
    r_content = float(geometry.radius * (fov_outer / 90.0))
    cx = float(geometry.center[0])
    cy = float(geometry.center[1])
    corners = (
        (float(rect.left()), float(rect.top())),
        (float(rect.right()), float(rect.top())),
        (float(rect.left()), float(rect.bottom())),
        (float(rect.right()), float(rect.bottom())),
    )
    r_window = max(math.hypot(x - cx, y - cy) for x, y in corners)
    r_max = float(max(r_content + 1.0, r_window))

    def pos(r: float) -> float:
        return max(0.0, min(1.0, r / r_max))

    col_params = {
        "transparent": (10, 7, 12, 8, 16, 11, 54, 8),
        "white": (246, 54, 246, 54, 246, 54, 255, 180),
        "black": (12, 9, 12, 9, 12, 9, 255, 180),
        "day": (230, 28, 242, 34, 255, 34, 200, 60),
        "night": (10, 7, 12, 9, 16, 11, 200, 60),
    }
    param = col_params.get(preset, None) or col_params["black"]
    def col(r: float, s: float) -> QColor:
        t = max(0.0, min(1.0, r / max(1.0, r_max)))
        rr = int(param[0] - param[1] * t)
        gg = int(param[2] - param[3] * t)
        bb = int(param[4] - param[5] * t)
        aa = int(param[6] * (1.0 - s) + param[7] * s)
        return QColor(rr, gg, bb, aa)

    c = geometry.center
    g = QRadialGradient(QPointF(c[0], c[1]), r_max)
    inner_color = (
        QColor(4, 4, 4, 112)
        if preset == "transparent"
        else QColor(4, 4, 4, 255) if preset in ("white", "day", "night", "black") else col(0.0, 0.0)
    )
    boundary_color = col(r_content, 0.3)
    edge_color = col(r_max, 1.0)
    g.setColorAt(0.0, inner_color)
    g.setColorAt(pos(r_content), inner_color)
    g.setColorAt(pos(r_content + 1.0), boundary_color)
    g.setColorAt(1.0, edge_color)

    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(g)
    painter.drawRect(rect)
    painter.restore()


def split_by_gaps(points: List[Tuple[float, float]]) -> List[List[Tuple[float, float]]]:
    """
    Split a polyline by large gaps to avoid drawing long, straight lines
    across the screen when a celestial path wraps around.

    Args:
        points: A list of (x, y) tuples representing the polyline.

    Returns:
        A list of polyline fragments, where each fragment is a list of points.
    """

    def dist(p1: Tuple[float, float], p2: Tuple[float, float]) -> float:
        return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

    fragments: List[List[Tuple[float, float]]] = [[]]
    for p in points:
        if not fragments[-1] or dist(p, fragments[-1][-1]) < 0.2:
            fragments[-1].append(p)
        else:
            fragments.append([p])
    return fragments


def _altaz_to_neu_unit(alt_deg: float, az_deg: float) -> np.ndarray:
    """Convert Alt/Az to a unit vector in local North-East-Up coordinates."""
    alt = math.radians(float(alt_deg))
    az = math.radians(float(az_deg))
    c_alt = math.cos(alt)
    return np.array(
        [
            c_alt * math.cos(az),  # north
            c_alt * math.sin(az),  # east
            math.sin(alt),         # up
        ],
        dtype=float,
    )


def _neu_unit_to_altaz(vec: np.ndarray) -> Tuple[float, float]:
    """Convert local North-East-Up unit vector to Alt/Az (deg)."""
    north = float(vec[0])
    east = float(vec[1])
    up = float(np.clip(float(vec[2]), -1.0, 1.0))
    alt_deg = math.degrees(math.asin(up))
    az_deg = math.degrees(math.atan2(east, north)) % 360.0
    return alt_deg, az_deg


def _great_circle_altaz_points(
    start_alt: float,
    start_az: float,
    end_alt: float,
    end_az: float,
) -> List[Tuple[float, float]]:
    """Sample great-circle points from start to end (both included)."""
    v0 = _altaz_to_neu_unit(start_alt, start_az)
    v1 = _altaz_to_neu_unit(end_alt, end_az)
    dot = float(np.clip(np.dot(v0, v1), -1.0, 1.0))
    omega = math.acos(dot)
    if omega < 1.0e-6:
        return [(float(start_alt), float(start_az)), (float(end_alt), float(end_az))]

    # Segment density grows with angular separation while staying bounded.
    sep_deg = math.degrees(omega)
    samples = max(8, min(64, int(sep_deg / 2.0) + 6))
    sin_omega = math.sin(omega)
    if abs(sin_omega) < 1.0e-8:
        return [(float(start_alt), float(start_az)), (float(end_alt), float(end_az))]

    out: List[Tuple[float, float]] = []
    for i in range(samples + 1):
        t = i / samples
        w0 = math.sin((1.0 - t) * omega) / sin_omega
        w1 = math.sin(t * omega) / sin_omega
        v = (w0 * v0) + (w1 * v1)
        norm = float(np.linalg.norm(v))
        if norm <= 1.0e-12:
            continue
        alt_i, az_i = _neu_unit_to_altaz(v / norm)
        out.append((alt_i, az_i))
    return out


def _clip_polyline_to_radius(
    points: List[Tuple[float, float]],
    max_radius: float,
) -> List[List[Tuple[float, float]]]:
    """Clip a normalized polyline to a circle centered at the origin."""
    if not points:
        return []

    radius_sq = max_radius * max_radius

    def _inside(point: Tuple[float, float]) -> bool:
        return (point[0] * point[0]) + (point[1] * point[1]) <= radius_sq

    def _intersections(
        start: Tuple[float, float],
        end: Tuple[float, float],
    ) -> List[Tuple[float, float, float]]:
        x0, y0 = start
        x1, y1 = end
        dx = x1 - x0
        dy = y1 - y0
        a = (dx * dx) + (dy * dy)
        if a <= 1.0e-12:
            return []
        b = 2.0 * ((x0 * dx) + (y0 * dy))
        c = (x0 * x0) + (y0 * y0) - radius_sq
        disc = (b * b) - (4.0 * a * c)
        if disc < 0.0:
            return []
        sqrt_disc = math.sqrt(max(0.0, disc))
        ts = [(-b - sqrt_disc) / (2.0 * a), (-b + sqrt_disc) / (2.0 * a)]
        hits: List[Tuple[float, float, float]] = []
        for t in ts:
            if 0.0 <= t <= 1.0:
                hits.append((t, x0 + (t * dx), y0 + (t * dy)))
        hits.sort(key=lambda hit: hit[0])
        unique: List[Tuple[float, float, float]] = []
        for hit in hits:
            if unique and abs(hit[0] - unique[-1][0]) < 1.0e-9:
                continue
            unique.append(hit)
        return unique

    fragments: List[List[Tuple[float, float]]] = []
    current: List[Tuple[float, float]] = [points[0]] if _inside(points[0]) else []

    for prev, cur in zip(points, points[1:]):
        prev_inside = _inside(prev)
        cur_inside = _inside(cur)
        hits = _intersections(prev, cur)

        if prev_inside and cur_inside:
            if not current:
                current = [prev]
            current.append(cur)
            continue

        if prev_inside and not cur_inside:
            if hits:
                hit = hits[-1]
                if not current:
                    current = [prev]
                current.append((hit[1], hit[2]))
            if len(current) >= 2:
                fragments.append(current)
            current = []
            continue

        if not prev_inside and cur_inside:
            if hits:
                hit = hits[0]
                current = [(hit[1], hit[2]), cur]
            else:
                current = [cur]
            continue

        if hits:
            first_hit = hits[0]
            second_hit = hits[-1]
            frag = [(first_hit[1], first_hit[2]), (second_hit[1], second_hit[2])]
            if len(frag) >= 2:
                fragments.append(frag)

    if len(current) >= 2:
        fragments.append(current)
    return fragments


def draw_sky_reference_lines(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    *,
    content_fov_deg: float | None = None,
) -> None:
    """
    Draw celestial reference lines like the equator, ecliptic, and horizon.

    Args:
        painter: The QPainter to use for drawing.
        geometry: The screen geometry for coordinate conversion.
        celestial_data: The data containing the points for the reference lines.
    """
    effective_fov_deg = _content_fov_deg_from_viewer(viewer_data) if content_fov_deg is None else float(content_fov_deg)
    point_list_pen_styles: List[Tuple[List[Tuple[float, float]], Tuple[QColor, int, List[int]]]] = [
        (celestial_data.celestial_equator_points, (CELESTIAL_EQUATOR_COLOR, 1, [8, 4])),
        (celestial_data.ecliptic_points, (ECLIPTIC_COLOR, 1, [3, 3])),
        (celestial_data.horizon_points, (HORIZON_LINE_COLOR, 1, [10, 1])),
    ]

    painter.save()
    for altaz_points, (color, width, style) in point_list_pen_styles:
        points: List[Tuple[float, float]] = []
        for alt, az in altaz_points:
            if not is_in_fov(float(alt), float(az), viewer_data.view_center, fov_deg=effective_fov_deg):
                continue
            nx, ny = altaz_to_normalized_xy(float(alt), float(az), viewer_data.view_center)
            points.append((nx, ny))
        for frag in split_by_gaps(points):
            if len(frag) < 2:
                continue
            pts = [QPointF(*normalized_to_screen_xy(nx, ny, geometry)) for nx, ny in frag]
            poly = QPolygonF(pts)

            # Use a tinted outline (same hue) instead of black to avoid dark gaps
            # on bright or highly transparent backgrounds.
            base_color = QColor(*color)
            base_color.setAlpha(70)
            base = QPen(base_color, width + 2, Qt.PenStyle.SolidLine)
            base.setCosmetic(True)
            # Keep the outline only under visible dash segments.
            # This avoids black "gaps" on bright/transparent backgrounds.
            base.setDashPattern(style)
            base.setCapStyle(Qt.PenCapStyle.RoundCap)
            base.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(base)
            painter.drawPolyline(poly)

            fg_color = QColor(*color)
            fg = QPen(fg_color, width)
            fg.setCosmetic(True)
            fg.setDashPattern(style)
            fg.setCapStyle(Qt.PenCapStyle.RoundCap)
            fg.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(fg)
            painter.drawPolyline(poly)
    painter.restore()


def draw_terrain_horizon_line(
    painter: QPainter,
    geometry: ScreenGeometry,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    view_center: tuple[float, float],
    *,
    opacity: float = 1.0,
    line_width_scale: float = 1.0,
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> None:
    """Draw a terrain horizon polyline as an extra overlay over the geometric horizon."""
    if not terrain_profile_altaz or opacity <= 0.0:
        return
    effective_opacity = max(0.0, min(1.0, opacity * 0.7))
    if effective_opacity <= 0.0:
        return

    points: list[tuple[float, float]] = []
    for alt, az in terrain_profile_altaz:
        if not is_in_fov(float(alt), float(az), view_center, fov_deg=content_fov_deg):
            continue
        nx, ny = altaz_to_normalized_xy(float(alt), float(az), view_center)
        points.append((nx, ny))

    if len(points) < 2:
        return

    color = QColor(*TERRAIN_HORIZON_LINE_COLOR)
    color.setAlphaF(max(effective_opacity, min(1.0, 0.42 + (opacity * 0.95))))
    outline = QColor(*TERRAIN_HORIZON_LINE_COLOR)
    outline.setAlpha(max(0, min(255, int(round(135.0 * effective_opacity + 35.0)))))
    width_scale = max(1.0, float(line_width_scale))
    painter.save()
    for frag in split_by_gaps(points):
        if len(frag) < 2:
            continue
        pts = [QPointF(*normalized_to_screen_xy(nx, ny, geometry)) for nx, ny in frag]
        poly = QPolygonF(pts)

        base = QPen(outline, 3.0 * width_scale, Qt.PenStyle.SolidLine)
        base.setCosmetic(True)
        base.setCapStyle(Qt.PenCapStyle.RoundCap)
        base.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        painter.setPen(base)
        painter.drawPolyline(poly)

        fg = QPen(color, 1.0 * width_scale, Qt.PenStyle.SolidLine)
        fg.setCosmetic(True)
        fg.setCapStyle(Qt.PenCapStyle.RoundCap)
        fg.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        painter.setPen(fg)
        painter.drawPolyline(poly)
    painter.restore()


def draw_urban_outlines(
    painter: QPainter,
    geometry: ScreenGeometry,
    urban_outlines: list[UrbanOutlinePolyline] | None,
    view_center: tuple[float, float],
    *,
    opacity: float = 0.2,
    line_width_scale: float = 1.0,
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> None:
    """Draw sampled building-top outlines directly on the sky dome."""
    if not urban_outlines:
        return
    if float(opacity) <= 0.0:
        return

    painter.save()
    width_scale = max(1.0, float(line_width_scale))
    for outline_entry in urban_outlines:
        outline = list(outline_entry.points)
        height_m = float(outline_entry.height_m)
        if len(outline) < 2:
            continue
        color = QColor(*URBAN_OUTLINE_LAYER_LINE_COLOR)
        effective_opacity = float(opacity) * _urban_outline_height_alpha_scale(height_m)
        color.setAlpha(max(0, min(255, int(round(255.0 * effective_opacity)))))
        pen = QPen(color, 1.5 * width_scale, Qt.PenStyle.SolidLine)
        pen.setCosmetic(True)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        simplified_pen = QPen(color, 4.0, Qt.PenStyle.SolidLine)
        simplified_pen.setCosmetic(True)
        simplified_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        simplified_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        azimuth_deg = [float(az) % 360.0 for _alt, az in outline]
        az_start_deg, az_end_deg, az_span_deg = _minimal_azimuth_cover(azimuth_deg)
        if az_span_deg < 0.5:
            representative_alt_deg = float(sum(float(alt) for alt, _az in outline) / len(outline))
            if representative_alt_deg >= -60.0:
                start_nx, start_ny = altaz_to_normalized_xy(representative_alt_deg, az_start_deg, view_center)
                end_nx, end_ny = altaz_to_normalized_xy(representative_alt_deg, az_end_deg, view_center)
                x1, y1 = normalized_to_screen_xy(start_nx, start_ny, geometry)
                x2, y2 = normalized_to_screen_xy(end_nx, end_ny, geometry)
                painter.setPen(simplified_pen)
                y = (float(y1) + float(y2)) * 0.5
                painter.drawPolyline(QPolygonF([QPointF(float(x1), y), QPointF(float(x2), y)]))
            continue

        painter.setPen(pen)
        points: list[tuple[float, float]] = []
        for alt, az in outline:
            if float(alt) < -60.0 or not is_in_fov(float(alt), float(az), view_center, fov_deg=content_fov_deg):
                if len(points) >= 2:
                    for frag in split_by_gaps(points):
                        if len(frag) < 2:
                            continue
                        screen_points = [QPointF(*normalized_to_screen_xy(nx, ny, geometry)) for nx, ny in frag]
                        painter.drawPolyline(QPolygonF(screen_points))
                points = []
                continue
            nx, ny = altaz_to_normalized_xy(float(alt), float(az), view_center)
            points.append((nx, ny))
        if len(points) >= 2:
            for frag in split_by_gaps(points):
                if len(frag) < 2:
                    continue
                screen_points = [QPointF(*normalized_to_screen_xy(nx, ny, geometry)) for nx, ny in frag]
                painter.drawPolyline(QPolygonF(screen_points))
    painter.restore()


def _urban_outline_height_alpha_scale(height_m: float) -> float:
    clamped_height_m = max(0.0, min(50.0, float(height_m)))
    return 0.25 + 0.75 * (clamped_height_m / 50.0)


def draw_satellite_overlay(
    painter: QPainter,
    geometry: ScreenGeometry,
    satellite_points: list[SatelliteOverlayPoint] | None,
    view_center: tuple[float, float],
    *,
    opacity: float = 1.0,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    preset: str = "night",
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> None:
    layer_opacity = max(0.0, min(1.0, float(opacity)))
    if not satellite_points or layer_opacity <= 0.0:
        return

    painter.save()
    marker_color = QColor(
        *SATELLITE_OVERLAY_MARKER_COLOR_RGB,
        max(0, min(255, int(round(SATELLITE_OVERLAY_MARKER_MAX_ALPHA * layer_opacity)))),
    )
    label_text_color, label_outline_color = get_text_style(preset)
    label_outline_width = 3.0
    label_text_color = QColor(label_text_color)
    label_outline_color = QColor(label_outline_color)
    label_text_color.setAlpha(max(0, min(255, int(round(label_text_color.alpha() * layer_opacity)))))
    label_outline_color.setAlpha(max(0, min(255, int(round(label_outline_color.alpha() * layer_opacity)))))
    for point in satellite_points:
        alt = float(point.alt_deg)
        az = float(point.az_deg)
        if not is_in_fov(alt, az, view_center, fov_deg=content_fov_deg):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        pos = QPointF(float(px), float(py))
        draw_gauge_cross(
            painter,
            marker_color,
            pos,
            scale=float(point.marker_scale),
            pen_width=1.0,
        )
        if bool(point.show_label):
            label_text = str(point.satellite_name).strip()
            if label_text:
                if label_candidates is not None:
                    label_candidates.append(
                        {
                            "text": label_text,
                            "pos": QPointF(pos.x() + 10.0, pos.y() - 8.0),
                            "text_color": label_text_color,
                            "outline_color": label_outline_color,
                            "outline_width": label_outline_width,
                            "priority": 42,
                            "hide_on_overlap": True,
                        }
                    )
    painter.restore()


def draw_aircraft_overlay(
    painter: QPainter,
    geometry: ScreenGeometry,
    aircraft_points: list[AircraftOverlayPoint] | None,
    view_center: tuple[float, float],
    *,
    opacity: float = 1.0,
    line_width_scale: float = 1.0,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    preset: str = "night",
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> None:
    """Draw aircraft as orange predicted-motion polylines."""
    layer_opacity = max(0.0, min(1.0, float(opacity)))
    if not aircraft_points or layer_opacity <= 0.0:
        return

    painter.save()
    width_scale = max(1.0, float(line_width_scale))
    line_color = QColor(*AIRCRAFT_OVERLAY_LINE_COLOR_RGB, 255)
    label_text_color, label_outline_color = get_text_style(preset)
    label_outline_width = 3.0
    label_text_color = QColor(label_text_color)
    label_outline_color = QColor(label_outline_color)
    label_text_color.setAlpha(max(0, min(255, int(round(label_text_color.alpha() * layer_opacity)))))
    label_outline_color.setAlpha(max(0, min(255, int(round(label_outline_color.alpha() * layer_opacity)))))
    line_pen = QPen(line_color, 1.0, Qt.PenStyle.SolidLine)
    line_pen.setCosmetic(True)
    line_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    line_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    painter.setPen(line_pen)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    for point in aircraft_points:
        alt = float(point.alt_deg)
        az = float(point.az_deg)
        distance_km = float(point.distance_km)
        if distance_km > _AIRCRAFT_MAX_DRAW_DISTANCE_KM:
            continue
        if alt <= 0.0 or not is_in_fov(alt, az, view_center, fov_deg=content_fov_deg):
            continue
        min_line_alpha = int(round(30.0 * layer_opacity))
        line_alpha = max(
            min_line_alpha,
            min(255, int(round(255.0 * float(point.alpha_scale) * layer_opacity))),
        )
        line_color.setAlpha(line_alpha)
        line_pen.setColor(line_color)
        line_pen.setWidthF(_aircraft_line_width_px(distance_km, width_scale=width_scale))
        painter.setPen(line_pen)
        trail_points = tuple(
            (float(sample_alt_deg), float(sample_az_deg))
            for sample_alt_deg, sample_az_deg in point.trail_alt_az_points
        )
        if any(
            is_in_fov(sample_alt_deg, sample_az_deg, view_center, fov_deg=content_fov_deg)
            for sample_alt_deg, sample_az_deg in trail_points
        ):
            polyline_points: list[QPointF] = []
            for sample_alt_deg, sample_az_deg in trail_points:
                sample_nx, sample_ny = altaz_to_normalized_xy(
                    sample_alt_deg,
                    sample_az_deg,
                    view_center,
                )
                sample_px, sample_py = normalized_to_screen_xy(sample_nx, sample_ny, geometry)
                polyline_points.append(QPointF(float(sample_px), float(sample_py)))
            if len(polyline_points) >= 2:
                painter.drawPolyline(QPolygonF(polyline_points))
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        pos = QPointF(float(px), float(py))
        callsign = (point.callsign or "").strip()
        if (
            callsign
            and float(point.distance_km) <= _AIRCRAFT_CALLSIGN_MAX_DISTANCE_KM
            and float(point.age_seconds) <= AIRCRAFT_FADE_START_SECONDS
        ):
            if label_candidates is not None:
                label_candidates.append(
                    {
                        "text": callsign,
                        "pos": QPointF(pos.x() + 8.0, pos.y() - 6.0),
                        "text_color": label_text_color,
                        "outline_color": label_outline_color,
                        "outline_width": label_outline_width,
                        "priority": 45,
                        "hide_on_overlap": True,
                    }
                )
    painter.restore()


_AIRCRAFT_CALLSIGN_MAX_DISTANCE_KM = 10.0
_AIRCRAFT_MAX_DRAW_DISTANCE_KM = 50.0


def _aircraft_line_width_px(distance_km: float, *, width_scale: float = 1.0) -> float:
    """Return a cosmetic line width where nearer aircraft appear thicker."""
    d = max(0.0, float(distance_km))
    scale = max(1.0, float(width_scale))
    aircraft_scale = 2.4 * scale
    if d <= 1.0:
        return 3.0 * aircraft_scale
    if d <= 3.0:
        return 2.2 * aircraft_scale
    if d <= 5.0:
        return 1.6 * aircraft_scale
    if d <= 10.0:
        return 1.0 * aircraft_scale
    if d <= 20.0:
        return 0.8 * aircraft_scale
    return 0.6 * aircraft_scale


def _minimal_azimuth_cover(azimuth_deg: List[float]) -> Tuple[float, float, float]:
    if not azimuth_deg:
        return 0.0, 0.0, 0.0
    if len(azimuth_deg) == 1:
        value = float(azimuth_deg[0]) % 360.0
        return value, value, 0.0

    values = sorted(float(value) % 360.0 for value in azimuth_deg)
    augmented = values + [values[0] + 360.0]
    largest_gap = -1.0
    gap_index = 0
    for index in range(len(values)):
        gap = augmented[index + 1] - augmented[index]
        if gap > largest_gap:
            largest_gap = gap
            gap_index = index
    start = augmented[gap_index + 1] % 360.0
    end = augmented[gap_index] % 360.0
    span = max(0.0, 360.0 - largest_gap)
    return start, end, span


def draw_asterisms(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    highlighted_object: Optional[Tuple[CelestialObject, QPointF]],
    text_font: QFont,
    label_reservations: Optional[List[QRectF]] = None,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    *,
    preset: str = "night",
    line_width_scale: float = 1.0,
    content_fov_deg: float | None = None,
    draw_base: bool = True,
    draw_highlight: bool = True,
) -> None:
    """Draw dim asterisms always, and brighten the hovered selection with a label."""

    stars = celestial_data.stars
    source_ids = resolve_star_source_ids(stars, celestial_data.star_catalog_meta)
    if source_ids.size == 0:
        return

    star_altaz_by_source: Dict[str, Tuple[float, float]] = {}
    for idx, raw_source in enumerate(source_ids):
        source_id = str(raw_source).strip()
        if not source_id:
            continue
        if source_id in star_altaz_by_source:
            continue
        star_altaz_by_source[source_id] = (float(stars["alt"][idx]), float(stars["az"][idx]))
    if not star_altaz_by_source:
        return

    if preset in ("white", "day"):
        base_line_color = QColor(26, 114, 214, 25)
        base_outline_color = QColor(190, 220, 250, 12)
        highlight_line_color = QColor(26, 114, 214, 124)
        highlight_outline_color = QColor(190, 220, 250, 52)
    else:
        base_line_color = QColor(82, 142, 214, 21)
        base_outline_color = QColor(24, 48, 86, 9)
        highlight_line_color = QColor(120, 190, 255, 92)
        highlight_outline_color = QColor(32, 76, 130, 44)

    painter.save()
    effective_fov_deg = _content_fov_deg_from_viewer(viewer_data) if content_fov_deg is None else float(content_fov_deg)
    clip_radius = effective_fov_deg / 90.0
    width_scale = max(1.0, float(line_width_scale))

    def _make_pen(color: QColor, width: float) -> QPen:
        pen = QPen(color, width)
        pen.setCosmetic(True)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        return pen

    def _draw_segments(
        segments: Iterable[Tuple[str, str]],
        outline_pen: QPen,
        line_pen: QPen,
    ) -> List[QPointF]:
        label_points: List[QPointF] = []
        for source_a, source_b in segments:
            pos_a = star_altaz_by_source.get(source_a)
            pos_b = star_altaz_by_source.get(source_b)
            if pos_a is None or pos_b is None:
                continue
            arc_altaz = _great_circle_altaz_points(pos_a[0], pos_a[1], pos_b[0], pos_b[1])
            arc_points: List[Tuple[float, float]] = []
            for alt_i, az_i in arc_altaz:
                nx_i, ny_i = altaz_to_normalized_xy(alt_i, az_i, viewer_data.view_center)
                arc_points.append((nx_i, ny_i))
            for raw_frag in split_by_gaps(arc_points):
                clipped_frags = _clip_polyline_to_radius(raw_frag, clip_radius)
                for frag in clipped_frags:
                    if len(frag) < 2:
                        continue
                    poly = QPolygonF([QPointF(*normalized_to_screen_xy(nx, ny, geometry)) for nx, ny in frag])
                    painter.setPen(outline_pen)
                    painter.drawPolyline(poly)
                    painter.setPen(line_pen)
                    painter.drawPolyline(poly)
                    label_points.extend(poly)
        return label_points

    def _draw_one_asterism(asterism: Any, outline_pen: QPen, line_pen: QPen) -> List[QPointF]:
        return _draw_segments(asterism.segments(), outline_pen, line_pen)

    if draw_base:
        base_outline_pen = _make_pen(base_outline_color, 4.0 * width_scale)
        base_line_pen = _make_pen(base_line_color, 2.5 * width_scale)
        base_segments: set[Tuple[str, str]] = set()
        for asterism in ASTERISMS:
            for source_a, source_b in asterism.segments():
                if source_a == source_b:
                    continue
                base_segments.add(tuple(sorted((source_a, source_b))))
        _draw_segments(sorted(base_segments), base_outline_pen, base_line_pen)

    highlighted_asterism = None
    if draw_highlight and highlighted_object is not None:
        hovered_obj, _ = highlighted_object
        if isinstance(hovered_obj, dict):
            hovered_source_id = str(hovered_obj.get("source_id", "")).strip()
            if hovered_source_id:
                second_slot = int(datetime.now().timestamp()) // 3
                highlighted_asterism = pick_rotating_asterism(hovered_source_id, second_slot)

    label_points: List[QPointF] = []
    if highlighted_asterism is not None:
        highlight_outline_pen = _make_pen(highlight_outline_color, 3.2 * width_scale)
        highlight_line_pen = _make_pen(highlight_line_color, 2.0 * width_scale)
        label_points = _draw_one_asterism(highlighted_asterism, highlight_outline_pen, highlight_line_pen)

    if label_points:
        cx = sum(pt.x() for pt in label_points) / len(label_points)
        cy = sum(pt.y() for pt in label_points) / len(label_points)
        label_pos = QPointF(cx + 8.0, cy - 8.0)
        if preset in ("white", "day"):
            text_color = QColor(*TEXT_STYLES_BY_PRESET["white"].text, 228)
        else:
            text_color = QColor(110, 195, 255, 230)
        _, outline_text_color = get_text_style(preset)
        outline_width = 3.0
        if label_candidates is not None:
            label_candidates.append(
                {
                    "text": highlighted_asterism.name,
                    "pos": label_pos,
                    "text_color": text_color,
                    "outline_color": outline_text_color,
                    "outline_width": outline_width,
                    "priority": 30,
                }
            )
        else:
            draw_outlined_text(
                painter,
                highlighted_asterism.name,
                label_pos,
                text_font,
                text_color,
                outline_text_color,
                outline_width=outline_width,
            )
            if label_reservations is not None:
                label_reservations.append(_text_bounds_at_baseline(highlighted_asterism.name, text_font, label_pos))

    painter.restore()


def draw_deep_sky_shapes(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    *,
    preset: str = "night",
) -> None:
    """Draw softly filled ellipses for named DSO rows that have shape metadata."""
    dso = celestial_data.deep_sky_objects
    if dso["alt"].size == 0:
        return

    names = np.asarray(dso["name"], dtype=object)
    ids = np.asarray(dso["id"], dtype=object)
    named_mask = np.array([_is_named_dso(n, i) for n, i in zip(names, ids)], dtype=bool)

    major = dso["major_arcmin"]
    minor = dso["minor_arcmin"]
    pa = dso["pa_deg"]
    finite_shape = (
        named_mask
        & np.isfinite(major)
        & np.isfinite(minor)
        & (major >= _DSO_SHAPE_MIN_MAJOR_ARCMIN)
        & (minor > 0.0)
    )
    if not np.any(finite_shape):
        return

    # Blue, soft fill style for named DSO extents.
    if preset in ("white", "day"):
        fill_rgb = (70, 140, 230)
    else:
        fill_rgb = (110, 185, 255)

    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)

    indices = np.nonzero(finite_shape)[0]
    for idx in indices:
        alt = float(dso["alt"][idx])
        az = float(dso["az"][idx])
        nx, ny = altaz_to_normalized_xy(alt, az, viewer_data.view_center)
        x, y = normalized_to_screen_xy(nx, ny, geometry)

        major_deg = float(major[idx]) / 60.0
        # Smaller/medium DSO need more opacity to remain recognizable;
        # very large objects (e.g., Magellanic Clouds) stay softer.
        alpha = int(np.clip(round(95.0 - 14.0 * math.sqrt(max(0.0, major_deg))), 56, 110))
        painter.setBrush(QColor(fill_rgb[0], fill_rgb[1], fill_rgb[2], alpha))
        poly = _dso_ellipse_polygon(
            alt_deg=alt,
            az_deg=az,
            major_arcmin=float(major[idx]),
            minor_arcmin=float(minor[idx]),
            pa_deg=float(pa[idx]) if np.isfinite(pa[idx]) else 0.0,
            viewer_data=viewer_data,
            geometry=geometry,
            gain=_DSO_SHAPE_SIZE_GAIN,
        )
        painter.drawPolygon(poly)

    painter.restore()


def draw_dso_hover_info(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    highlighted_dso: Optional[Tuple[CelestialObject, QPointF]],
    text_font: QFont,
    *,
    preset: str = "night",
) -> None:
    """Draw DSO hover overlay shape (intended to be placed behind stars).

    Hover style keeps the filled body at catalog size and adds a larger
    (3x) blue outline ring for emphasis.
    """
    if not highlighted_dso:
        return
    obj, _ = highlighted_dso
    major_arcmin = float(obj.get("major_arcmin", 0.0))
    minor_arcmin = float(obj.get("minor_arcmin", 0.0))
    pa_deg = float(obj.get("pa_deg", 0.0))
    alt = float(obj.get("alt", 0.0))
    az = float(obj.get("az", 0.0))
    if preset in ("white", "day"):
        hover_pen = QColor(40, 122, 220, 220)
        hover_fill = QColor(70, 140, 230, 70)
    else:
        hover_pen = QColor(110, 195, 255, 230)
        hover_fill = QColor(110, 185, 255, 62)
    base_poly = _dso_ellipse_polygon(
        alt_deg=alt,
        az_deg=az,
        major_arcmin=major_arcmin,
        minor_arcmin=minor_arcmin,
        pa_deg=pa_deg if math.isfinite(pa_deg) else 0.0,
        viewer_data=viewer_data,
        geometry=geometry,
        gain=_DSO_SHAPE_SIZE_GAIN,
        samples=60,
    )
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(hover_fill)
    painter.drawPolygon(base_poly)

    hover_poly = _dso_ellipse_polygon(
        alt_deg=alt,
        az_deg=az,
        major_arcmin=major_arcmin,
        minor_arcmin=minor_arcmin,
        pa_deg=pa_deg if math.isfinite(pa_deg) else 0.0,
        viewer_data=viewer_data,
        geometry=geometry,
        gain=_DSO_HOVER_SIZE_GAIN,
        samples=60,
    )
    painter.setPen(QPen(hover_pen, 2.6))
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawPolygon(hover_poly)


def draw_stars(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    star_base_radius: float,
    *,
    visibility_boost: float = 1.0,
    draw_vmag_limit: Optional[float] = None,
    viewport_size: Tuple[int, int] | None = None,
    content_fov_deg: float | None = None,
) -> None:
    """
    Draw stars using a numpy canvas that paints uniformly sized rectangles.

    The visual area of each star is derived from its magnitude (a 1st-magnitude
    star maps to a 6x6 rectangle), and the color intensity is reduced for stars
    whose calculated area falls below 1 pixel^2. The resulting RGB canvas is
    converted into a `QImage` and blended additively to the painter before
    drawing other objects.

    Args:
        painter: The QPainter object for drawing.
        geometry: The screen geometry for coordinate conversions.
        celestial_data: The data containing star information (positions, magnitudes, etc.).
        viewer_data: The viewer's data, including the view center.
        star_base_radius: A base radius for scaling star sizes.
        visibility_boost: Multiplier for the overall star colors.
        draw_vmag_limit: Optional override to skip fainter stars entirely.
        viewport_size: Optional `(width, height)` of the drawing area, used to create the numpy canvas; if omitted, defaults to the painter's clip rect.
    """
    stars = celestial_data.stars
    effective_fov_deg = _content_fov_deg_from_viewer(viewer_data) if content_fov_deg is None else float(content_fov_deg)
    visibility_boost = float(np.clip(visibility_boost, 0.7, 2.0))

    if draw_vmag_limit is not None:
        draw_mask = stars["vmag"] <= float(draw_vmag_limit)
        if not np.any(draw_mask):
            return
        alt = stars["alt"][draw_mask]
        az = stars["az"][draw_mask]
        vmag = stars["vmag"][draw_mask]
        bv = stars["bv"][draw_mask]
        size_factor = stars["size_factor"][draw_mask]
        color_factor_base = stars["color_factor_base"][draw_mask]
    else:
        alt = stars["alt"]
        az = stars["az"]
        vmag = stars["vmag"]
        bv = stars["bv"]
        size_factor = stars["size_factor"]
        color_factor_base = stars["color_factor_base"]

    width_px = None
    height_px = None
    if viewport_size:
        width_px = max(1, int(viewport_size[0]))
        height_px = max(1, int(viewport_size[1]))
    else:
        width_px = painter.viewport().width()
        height_px = painter.viewport().height()

    if width_px <= 0 or height_px <= 0:
        return

    # Convert directions to pixel centers.
    nx, ny = _altaz_to_normalized_xy_vectorized(alt, az, viewer_data.view_center)
    x, y = _normalized_to_screen_xy_vectorized(nx, ny, geometry)

    bv_clamped = np.nan_to_num(bv, nan=0.45)
    rgb_colors = bv_to_rgb_vectorized(bv_clamped).astype(np.float32)

    # `star_base_radius` is defined as the apparent square size of a 2nd-magnitude star.
    max_size = max(12.0, float(max(1.0, star_base_radius)))
    size_float = float(star_base_radius) * _MAG2_TO_MAG1_SIZE_SCALE * size_factor
    size_px = np.clip(np.round(size_float), 1, int(max_size)).astype(int)

    color_factor = np.clip(0.15 + 0.85 * color_factor_base, 0.0, 1.0) * visibility_boost
    color_factor = np.clip(color_factor, 0.0, 1.0)
    drawn_area = np.maximum(1.0, size_px.astype(float) ** 2)
    expected_area = size_float**2
    area_ratio = np.clip(expected_area / drawn_area, 0.0, 1.0)
    star_colors = rgb_colors * color_factor[:, None] * area_ratio[:, None]

    cache_key = _star_cache_key(
        alt,
        az,
        size_factor,
        color_factor_base,
        celestial_data.time.jd,
        viewer_data.view_center,
        geometry,
        star_base_radius,
        visibility_boost,
        draw_vmag_limit,
    )
    global _star_render_cache
    if _star_render_cache and _star_render_cache[0] == cache_key:
        painter.save()
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
        painter.drawImage(0, 0, _star_render_cache[1])
        painter.restore()
        return

    canvas = np.zeros((height_px, width_px, 3), dtype=np.float32)

    ix = np.round(x).astype(int)
    iy = np.round(y).astype(int)

    half = (size_px // 2).astype(int)
    x0 = ix - half
    y0 = iy - half
    x1 = x0 + size_px
    y1 = y0 + size_px

    x0_clamped = np.clip(x0, 0, width_px)
    y0_clamped = np.clip(y0, 0, height_px)
    x1_clamped = np.clip(x1, 0, width_px)
    y1_clamped = np.clip(y1, 0, height_px)

    outside_content = ~is_in_fov_vectorized(alt, az, viewer_data.view_center, fov_deg=effective_fov_deg)
    valid_base = (
        (x1_clamped > x0_clamped)
        & (y1_clamped > y0_clamped)
        & (size_px > 0)
        & (~outside_content)
    )
    size_one = size_px == 1
    size_two = size_px == 2
    single_mask = valid_base & size_one
    size2_full_fit = (x0 >= 0) & (y0 >= 0) & (x1 <= width_px) & (y1 <= height_px)
    size2_mask = valid_base & size_two & size2_full_fit
    multi_mask = valid_base & ~(size_one | size_two)
    valid = single_mask | size2_mask | multi_mask

    single_indices = np.nonzero(single_mask)[0]
    if single_indices.size > 0:
        single_layer = np.zeros_like(canvas)
        flat_single = single_layer.reshape(-1, 3)
        x_single = x0_clamped[single_indices]
        y_single = y0_clamped[single_indices]
        flat_idx = y_single * width_px + x_single
        np.add.at(flat_single, flat_idx, star_colors[single_indices])
        canvas += _apply_weak_gaussian3x3_rgb(single_layer, _SINGLE_STAR_GAUSSIAN_STRENGTH)

    size2_indices = np.nonzero(size2_mask)[0]
    if size2_indices.size > 0:
        size2_layer = np.zeros_like(canvas)
        flat_size2 = size2_layer.reshape(-1, 3)
        x_size2 = x0[size2_indices]
        y_size2 = y0[size2_indices]
        base_idx = y_size2 * width_px + x_size2
        colors_size2 = star_colors[size2_indices]
        np.add.at(flat_size2, base_idx, colors_size2)
        np.add.at(flat_size2, base_idx + 1, colors_size2)
        np.add.at(flat_size2, base_idx + width_px, colors_size2)
        np.add.at(flat_size2, base_idx + width_px + 1, colors_size2)
        canvas += size2_layer

    multi_indices = np.nonzero(multi_mask)[0]
    for idx in multi_indices:
        canvas[y0_clamped[idx]:y1_clamped[idx], x0_clamped[idx]:x1_clamped[idx], :] += star_colors[idx]

    # Overlay a rotated square (diamond) for bright stars to emphasize stars above 2nd magnitude.
    bright_indices = np.nonzero(valid & (vmag < 2.0))[0]
    if bright_indices.size > 0:
        bright_scale = np.power(10.0, 0.12 * np.clip(2.0 - vmag[bright_indices], 0.0, 4.0))
        for local_i, idx in enumerate(bright_indices):
            cx = float(ix[idx])
            cy = float(iy[idx])
            half_diag = max(0.5, 0.5 * float(size_px[idx]) * float(bright_scale[local_i]))
            xmin = max(0, int(math.floor(cx - half_diag)))
            xmax = min(width_px, int(math.ceil(cx + half_diag + 1.0)))
            ymin = max(0, int(math.floor(cy - half_diag)))
            ymax = min(height_px, int(math.ceil(cy + half_diag + 1.0)))
            if xmin >= xmax or ymin >= ymax:
                continue
            xs = np.arange(xmin, xmax, dtype=np.float32)[None, :]
            ys = np.arange(ymin, ymax, dtype=np.float32)[:, None]
            mask = (np.abs(xs - cx) + np.abs(ys - cy)) <= half_diag
            if not np.any(mask):
                continue
            canvas[ymin:ymax, xmin:xmax, :][mask] += star_colors[idx] * _DIAMOND_OVERLAY_GAIN

    np.clip(canvas, 0.0, 255.0, out=canvas)
    canvas_uint8 = np.ascontiguousarray(canvas.astype(np.uint8))
    if canvas_uint8.size == 0:
        return

    rgba = np.zeros((height_px, width_px, 4), dtype=np.uint8)
    alpha = np.max(canvas_uint8, axis=2)
    nonzero_alpha = alpha > 0
    if np.any(nonzero_alpha):
        rgb_float = canvas_uint8.astype(np.float32) / 255.0
        alpha_float = alpha.astype(np.float32) / 255.0
        rgba_rgb = np.zeros_like(rgb_float, dtype=np.float32)
        rgba_rgb[nonzero_alpha] = rgb_float[nonzero_alpha] / alpha_float[nonzero_alpha, None]
        np.clip(rgba_rgb, 0.0, 1.0, out=rgba_rgb)
        rgba[:, :, :3] = np.round(rgba_rgb * 255.0).astype(np.uint8)
    rgba[:, :, 3] = alpha
    image = QImage(rgba.data, width_px, height_px, width_px * 4, QImage.Format_RGBA8888).copy()
    painter.save()
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
    painter.drawImage(0, 0, image)
    logger.debug("Generated star buffer (%d stars) and cached image", len(size_px))
    _star_render_cache = (cache_key, image)

    painter.restore()
    # Reset composition mode.
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)


def draw_gauge_cross(
    painter: QPainter,
    color: QColor,
    center: QPointF,
    *,
    scale: float = 1.0,
    pen_width: float = 1.0,
) -> None:
    """
    Draws a cross-shaped gauge marker.

    This is used to indicate the position of certain celestial objects, like
    the Sun or the Moon.

    Args:
        painter: The QPainter to use for drawing.
        color: The color of the cross.
        center: The center point (QPointF) of the cross.
    """
    scale = max(0.05, float(scale))
    cross_outer_len = max(1, int(round(15 * scale)))
    cross_inner_len = max(0, min(cross_outer_len - 1, int(round(4 * scale))))
    x, y = center.x(), center.y()
    painter.setPen(QPen(color, float(pen_width)))
    painter.drawLine(QPointF(x - cross_outer_len, y), QPointF(x - cross_inner_len, y))
    painter.drawLine(QPointF(x + cross_inner_len, y), QPointF(x + cross_outer_len, y))
    painter.drawLine(QPointF(x, y - cross_outer_len), QPointF(x, y - cross_inner_len))
    painter.drawLine(QPointF(x, y + cross_inner_len), QPointF(x, y + cross_outer_len))


def draw_zenith_marker(
    painter: QPainter,
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    *,
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> None:
    """
    Draws markers at zenith and nadir.

    Args:
        painter: The QPainter to use for drawing.
        geometry: The screen geometry for coordinate conversion.
        view_center: The current view center (altitude, azimuth).
    """
    az_ref = view_center[1]
    s = 7
    painter.setPen(QPen(QColor(*TEXT_STYLES_BY_PRESET["night"].text), 1))
    for alt in (90.0, -90.0):
        if not is_in_fov(alt, az_ref, view_center, fov_deg=content_fov_deg):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az_ref, view_center)
        x, y = normalized_to_screen_xy(nx, ny, geometry)
        painter.drawLine(QPointF(x - s, y - s), QPointF(x + s, y + s))
        painter.drawLine(QPointF(x - s, y + s), QPointF(x + s, y - s))


def draw_moon(
    painter: QPainter,
    center: QPointF,
    radius_px: float,
    sun_dir_in_moon_frame: np.ndarray,
    screen_rotation_deg: float,
    opacity: float = 1.0,
    base_color: Optional[QColor] = None,
) -> None:
    """
    Draw the moon with its correct phase and orientation.

    The moon's phase is determined by the sun's direction vector relative to
    the moon. The image is also rotated to match the screen orientation.

    Args:
        painter: The QPainter for drawing.
        center: The center position of the moon on the screen.
        radius_px: The radius of the moon in pixels.
        sun_dir_in_moon_frame: The direction vector of the sun in the moon's reference frame.
        screen_rotation_deg: The rotation angle of the screen in degrees.
        opacity: The opacity of the moon image.
        base_color: An optional QColor to tint the moon, used for eclipses.
    """
    img_size = max(5, int(math.ceil(radius_px * 2.0)))
    view_dir = np.array([0, 0, 1], dtype=float)
    if base_color is not None:
        tint_rgba = (base_color.red(), base_color.green(), base_color.blue(), base_color.alpha())
    else:
        tint_rgba = None
    moon_rgba = generate_moon_phase_rgba(img_size, sun_dir_in_moon_frame, view_dir, tint_color=tint_rgba)
    moon_image = np_rgba_to_qimage(moon_rgba)
    source_rect = QRectF(0, 0, img_size, img_size)

    painter.save()
    painter.setOpacity(opacity)
    painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    painter.translate(center)
    if abs(screen_rotation_deg) > 0.1:
        painter.rotate(screen_rotation_deg)
    target_rect = QRectF(-img_size / 2, -img_size / 2, img_size, img_size)
    painter.drawImage(target_rect, moon_image, source_rect)

    if base_color is not None:
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceAtop)
        painter.setBrush(base_color)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawEllipse(QPointF(0.0, 0.0), radius_px, radius_px)

    painter.restore()


def _collect_sun_moon_context(planets: List[PlanetBody]) -> Tuple[Optional[PlanetBody], Optional[Tuple[float, float]], Optional[Tuple[float, float]]]:
    """Collect moon body and sun/moon alt-az pairs from planet list."""
    moon_body: Optional[PlanetBody] = None
    sun_altaz: Optional[Tuple[float, float]] = None
    moon_altaz: Optional[Tuple[float, float]] = None
    for body in planets:
        if body.name == "sun":
            sun_altaz = (body.alt, body.az)
        elif body.name == "moon":
            moon_altaz = (body.alt, body.az)
            moon_body = body
    return moon_body, sun_altaz, moon_altaz


def _moon_eclipse_overlay_color(body: PlanetBody) -> Optional[QColor]:
    """Return overlay color for lunar eclipse rendering, if needed."""
    eclipse = body.lunar_eclipse_info
    if not eclipse or not eclipse.is_eclipse:
        return None
    if eclipse.eclipse_type == "partial":
        return QColor(30, 0, 0, 60)
    if eclipse.eclipse_type == "total":
        return QColor(40, 10, 10, 180)
    return None


def _draw_moon_planet(
    painter: QPainter,
    pos: QPointF,
    geometry: ScreenGeometry,
    body: PlanetBody,
    viewer_data: ViewerData,
    sun_altaz: Tuple[float, float],
    moon_altaz: Tuple[float, float],
    enlarge_moon: bool,
    cross_color: QColor,
) -> None:
    """Draw moon phase disc (with eclipse tint) and its gauge marker."""
    moon_zoom = 5 if enlarge_moon else 1
    sun_dir_in_moon_frame, screen_rotation_deg = calculate_moon_render_data(sun_altaz, moon_altaz, viewer_data.view_center)
    base_moon_radius_px = max((0.25 / 90.0) * geometry.radius, 2.5)
    moon_radius_px = base_moon_radius_px * moon_zoom
    draw_moon(
        painter,
        pos,
        moon_radius_px,
        sun_dir_in_moon_frame=sun_dir_in_moon_frame,
        screen_rotation_deg=screen_rotation_deg,
        opacity=1.0 if not enlarge_moon else 0.7,
        base_color=_moon_eclipse_overlay_color(body),
    )
    draw_gauge_cross(painter, cross_color, pos)


def _marker_intersects_viewport(painter: QPainter, pos: QPointF, radius_px: float) -> bool:
    """Return whether a body marker would be visible inside the current viewport."""
    viewport_getter = getattr(painter, "viewport", None)
    if not callable(viewport_getter):
        return True
    viewport = viewport_getter()
    if viewport.isNull():
        return True
    r = max(1.0, float(radius_px))
    marker_bounds = QRectF(pos.x() - r, pos.y() - r, r * 2.0, r * 2.0)
    return marker_bounds.intersects(QRectF(viewport))


def draw_planet_disc(
    painter: QPainter,
    pos: QPointF,
    color: QColor,
    *,
    radius_px: float,
    alpha: int,
) -> None:
    """Draw a no-flare solid circular marker for planets."""
    r = max(1.5, float(radius_px))
    fill = QColor(color)
    fill.setAlpha(int(np.clip(alpha, 1, 255)))
    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(fill)
    painter.drawEllipse(pos, r, r)
    painter.restore()


def draw_planet_bloom(
    painter: QPainter,
    pos: QPointF,
    color: QColor,
    *,
    core_radius_px: float,
    vmag: Optional[float],
) -> None:
    """Draw a soft additive bloom around a planet marker."""
    bloom_radius, center_alpha, mid_alpha = planet_bloom_profile_from_vmag(vmag, core_radius_px)
    if bloom_radius <= 0.0 or center_alpha <= 0:
        return

    c0 = QColor(color)
    c0.setAlpha(center_alpha)
    c1 = QColor(color)
    c1.setAlpha(mid_alpha)
    c2 = QColor(color)
    c2.setAlpha(0)

    gradient = QRadialGradient(pos, bloom_radius)
    gradient.setColorAt(0.0, c0)
    gradient.setColorAt(0.45, c1)
    gradient.setColorAt(1.0, c2)

    painter.save()
    # Additive blend to emulate bloom-like luminous spill.
    painter.setCompositionMode(QPainter.CompositionMode_Plus)
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(gradient)
    painter.drawEllipse(pos, bloom_radius, bloom_radius)
    painter.restore()


def draw_solar_system_bodies(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    enlarge_moon: bool,
    *,
    text_font: Optional[QFont] = None,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    label_reservations: Optional[List[QRectF]] = None,
    draw_markers: bool = True,
    draw_labels: bool = True,
    preset: str = "night",
    content_fov_deg: float | None = None,
) -> None:
    """
    Draw major solar system bodies (Sun, Moon, and planets).

    This function iterates through bodies in the sky data, calculates their
    screen positions, and draws them. The Sun is drawn as a gauge cross, the
    Moon is drawn with its phase. Other planets are represented by circular markers.

    Args:
        painter: The QPainter for drawing.
        geometry: The screen geometry for coordinate conversion.
        celestial_data: The data containing solar system body information.
        viewer_data: The viewer's data for position calculations.
        enlarge_moon: A boolean indicating whether to draw the moon larger.
    """
    moon_body, sun_altaz, moon_altaz = _collect_sun_moon_context(celestial_data.planets)
    effective_fov_deg = _content_fov_deg_from_viewer(viewer_data) if content_fov_deg is None else float(content_fov_deg)

    # Keep body markers in a stable high-contrast color over the sky disc.
    text_color = QColor(*TEXT_STYLES_BY_PRESET["night"].text)
    label_text_color, label_outline_color = get_text_style(preset)
    label_outline_width = 3.0
    if text_font is not None:
        painter.setFont(text_font)
        label_font = text_font
    else:
        label_font = painter.font() if hasattr(painter, "font") else QFont()

    for body in celestial_data.planets:
        # Draw planets if they are in-view, even below horizon.
        if not is_in_fov(body.alt, body.az, viewer_data.view_center, fov_deg=effective_fov_deg):
            continue

        pos = QPointF(
            *normalized_to_screen_xy(
                *altaz_to_normalized_xy(body.alt, body.az, viewer_data.view_center),
                geometry,
            )
        )
        marker_visible = True
        if body.name == "moon":
            base_moon_radius_px = max((0.25 / 90.0) * geometry.radius, 2.5)
            moon_zoom = 5 if enlarge_moon else 1
            marker_visible = _marker_intersects_viewport(
                painter,
                pos,
                base_moon_radius_px * moon_zoom,
            )
        else:
            radius_px, _alpha = planet_disc_style_from_vmag(body.vmag)
            bloom_radius, _center_alpha, _mid_alpha = planet_bloom_profile_from_vmag(body.vmag, radius_px)
            marker_visible = _marker_intersects_viewport(painter, pos, max(radius_px, bloom_radius))

        if draw_markers:
            if body.name == "sun":
                draw_gauge_cross(painter, text_color, pos)
            elif body.name == "moon" and moon_body and sun_altaz and moon_altaz:
                _draw_moon_planet(
                    painter,
                    pos,
                    geometry,
                    body,
                    viewer_data,
                    sun_altaz,
                    moon_altaz,
                    enlarge_moon,
                    text_color,
                )
            else:
                radius_px, alpha = planet_disc_style_from_vmag(body.vmag)
                marker_color = planet_marker_color(body.name)
                draw_planet_bloom(painter, pos, marker_color, core_radius_px=radius_px, vmag=body.vmag)
                marker_color.setAlpha(alpha)
                draw_planet_disc(painter, pos, marker_color, radius_px=radius_px, alpha=alpha)
                # Keep planet cross markers shorter than the Moon/Sun gauge marker.
                draw_gauge_cross(painter, text_color, pos, scale=0.55, pen_width=1.0)

        if draw_labels and body.name != "sun" and marker_visible:
            label_text = body_label_text(body.name)
            label_pos = QPointF(pos.x() + 12.0, pos.y() - 10.0)
            if label_candidates is not None:
                label_candidates.append(
                    {
                        "text": label_text,
                        "pos": label_pos,
                        "text_color": label_text_color,
                        "outline_color": label_outline_color,
                        "outline_width": label_outline_width,
                        "priority": 40,
                        "hide_on_overlap": True,
                    }
                )
                continue
            if label_reservations is not None:
                rect = _text_bounds_at_baseline(label_text, label_font, label_pos)
                if _rect_overlap_count(rect, label_reservations) > 0:
                    continue
                label_reservations.append(rect)
            draw_outlined_text(
                painter,
                label_text,
                label_pos,
                label_font,
                label_text_color,
                label_outline_color,
                    outline_width=label_outline_width,
                )


def draw_hovered_moon_overlay(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    highlighted_object: Optional[Tuple[CelestialObject, QPointF]],
) -> None:
    """Redraw the hovered moon enlarged on top of the normal scene."""
    if highlighted_object is None:
        return
    obj, pos = highlighted_object
    obj_name = getattr(obj, "name", None) if hasattr(obj, "name") else obj.get("name")
    if str(obj_name).strip().lower() != "moon":
        return
    moon_body, sun_altaz, moon_altaz = _collect_sun_moon_context(celestial_data.planets)
    if moon_body is None or sun_altaz is None or moon_altaz is None:
        return
    text_color = QColor(*TEXT_STYLES_BY_PRESET["night"].text)
    _draw_moon_planet(
        painter,
        pos,
        geometry,
        moon_body,
        viewer_data,
        sun_altaz,
        moon_altaz,
        True,
        text_color,
    )


def draw_direction_labels(
    painter: QPainter,
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    text_font: QFont,
    mouse_pos: Optional[QPoint] = None,
    *,
    preset: str = "night",
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> None:
    """
    Draw compass direction labels and horizon markers on the horizon.

    Args:
        painter: The QPainter for drawing.
        geometry: The screen geometry for coordinate conversion.
        view_center: The current view center to determine which labels are visible.
        text_font: The QFont to use for the labels.
    """
    # Match direction labels to the horizon line color.
    text_color = QColor(*HORIZON_LINE_COLOR)
    outline_color = QColor.fromRgbF(0, 0, 0, 0.3)
    marker_color = QColor(*HORIZON_LINE_COLOR)
    marker_pen = QPen(marker_color, 1.6)
    marker_pen.setCosmetic(True)
    marker_half_len_px = 6.0  # constant visual length regardless of view altitude
    marker_hit_radius_px = 4.0
    tangent_probe_deg = 0.6
    label_outward_offset_px = 4.0
    hide_near_mouse_px = 28.0
    mouse_x = float(mouse_pos.x()) if mouse_pos is not None else None
    mouse_y = float(mouse_pos.y()) if mouse_pos is not None else None
    painter.setFont(text_font)
    fm = QFontMetrics(text_font)
    alt = 0.0
    for label, az in DIRECTIONS.items():
        if not is_in_fov(alt, az, view_center, fov_deg=content_fov_deg):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        pos = QPointF(*normalized_to_screen_xy(nx, ny, geometry))

        # Estimate local horizon tangent from azimuth finite difference.
        az_prev = (az - tangent_probe_deg + 360.0) % 360.0
        az_next = (az + tangent_probe_deg) % 360.0
        p_prev = QPointF(*normalized_to_screen_xy(*altaz_to_normalized_xy(alt, az_prev, view_center), geometry))
        p_next = QPointF(*normalized_to_screen_xy(*altaz_to_normalized_xy(alt, az_next, view_center), geometry))
        tx = p_next.x() - p_prev.x()
        ty = p_next.y() - p_prev.y()
        t_norm = math.hypot(tx, ty)
        if t_norm <= 1e-6:
            rx = pos.x() - geometry.center[0]
            ry = pos.y() - geometry.center[1]
            tx, ty = -ry, rx
            t_norm = math.hypot(tx, ty)
        if t_norm <= 1e-6:
            tx, ty, t_norm = 1.0, 0.0, 1.0
        ux, uy = tx / t_norm, ty / t_norm
        # Draw a marker line perpendicular to the local horizon tangent.
        nxp, nyp = -uy, ux
        painter.save()
        painter.setPen(marker_pen)
        painter.drawLine(
            QPointF(pos.x() - nxp * marker_half_len_px, pos.y() - nyp * marker_half_len_px),
            QPointF(pos.x() + nxp * marker_half_len_px, pos.y() + nyp * marker_half_len_px),
        )
        painter.restore()

        label_pos = QPointF(pos)
        # Baseline-relative bounds for draw_outlined_text(path.addText baseline semantics).
        bounds = fm.tightBoundingRect(label)
        label_rect = QRectF(
            label_pos.x() + bounds.x(),
            label_pos.y() + bounds.y(),
            bounds.width(),
            bounds.height(),
        )
        # If label and marker overlap, push label slightly outward from disc center.
        nearest_x = min(max(pos.x(), label_rect.left()), label_rect.right())
        nearest_y = min(max(pos.y(), label_rect.top()), label_rect.bottom())
        dx0 = pos.x() - nearest_x
        dy0 = pos.y() - nearest_y
        overlap = (dx0 * dx0 + dy0 * dy0) <= (marker_hit_radius_px**2)
        if overlap:
            ox = pos.x() - geometry.center[0]
            oy = pos.y() - geometry.center[1]
            norm = math.hypot(ox, oy)
            if norm > 1e-6:
                label_pos = QPointF(
                    label_pos.x() + (ox / norm) * label_outward_offset_px,
                    label_pos.y() + (oy / norm) * label_outward_offset_px,
                )

        if mouse_x is not None and mouse_y is not None:
            dx = label_pos.x() - mouse_x
            dy = label_pos.y() - mouse_y
            if (dx * dx + dy * dy) <= (hide_near_mouse_px * hide_near_mouse_px):
                continue

        draw_outlined_text(
            painter,
            label,
            label_pos,
            text_font,
            text_color,
            outline_color,
            outline_width=2.5,
        )


def draw_overlay_info(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    vmag_limit: float,
    enlarge_moon: bool,
    highlighted_dso: Optional[Tuple[CelestialObject, QPointF]],
    highlighted_object: Optional[Tuple[CelestialObject, QPointF]],
    text_font: QFont,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    label_reservations: Optional[List[QRectF]] = None,
    *,
    preset: str = "night",
    draw_static_info: bool = True,
    draw_hover_info: bool = True,
) -> None:
    """
    Draws overlay text information on the screen.

    This includes the current time, location, view direction, magnitude limit,
    and information about any highlighted celestial object.

    Args:
        painter: The QPainter for drawing.
        celestial_data: The data containing the current time.
        viewer_data: The data containing viewer location and view direction.
        vmag_limit: The current visual magnitude limit for stars.
        enlarge_moon: A boolean indicating if the moon is enlarged.
        highlighted_object: The currently highlighted object, if any.
        text_font: The QFont to use for the text.
    """
    text_color, outline_color = get_text_style(preset)
    outline_width = 3.0

    line_spacing = QFontMetrics(text_font).lineSpacing()
    line_height = int(line_spacing * 1.2)
    line_x = line_spacing
    line_y = 0

    def print_line(message: str):
        nonlocal line_x, line_y
        line_y += line_height
        draw_outlined_text(
            painter,
            message,
            QPointF(line_x, line_y),
            text_font,
            text_color,
            outline_color,
            outline_width=outline_width,
        )

    if draw_static_info:
        for line in format_overlay_info_lines(celestial_data, viewer_data, vmag_limit):
            print_line(line)

    # ---- DSO label (draw first so star/planet labels stay in front among labels) ----
    if draw_hover_info and highlighted_dso:
        dso_obj, _ = highlighted_dso
        major_arcmin = float(dso_obj.get("major_arcmin", 0.0))
        minor_arcmin = float(dso_obj.get("minor_arcmin", 0.0))
        pa_deg = float(dso_obj.get("pa_deg", 0.0))
        alt = float(dso_obj.get("alt", 0.0))
        az = float(dso_obj.get("az", 0.0))
        if preset in ("white", "day"):
            dso_label_color = QColor(*TEXT_STYLES_BY_PRESET["white"].text, 228)
        else:
            dso_label_color = QColor(110, 195, 255, 230)
        hover_poly = _dso_ellipse_polygon(
            alt_deg=alt,
            az_deg=az,
            major_arcmin=major_arcmin,
            minor_arcmin=minor_arcmin,
            pa_deg=pa_deg if math.isfinite(pa_deg) else 0.0,
            viewer_data=viewer_data,
            geometry=geometry,
            gain=_DSO_HOVER_SIZE_GAIN,
            samples=60,
        )
        dso_name = str(dso_obj.get("name", "")).strip()
        if dso_name:
            bounds = hover_poly.boundingRect()
            label_pos = QPointF(bounds.right() + 10.0, bounds.top() - 6.0)
            if label_candidates is not None:
                label_candidates.append(
                    {
                        "text": dso_name,
                        "pos": label_pos,
                        "text_color": dso_label_color,
                        "outline_color": outline_color,
                        "outline_width": outline_width,
                        "priority": 20,
                    }
                )
            else:
                draw_outlined_text(
                    painter,
                    dso_name,
                    label_pos,
                    text_font,
                    dso_label_color,
                    outline_color,
                    outline_width=outline_width,
                )
                if label_reservations is not None:
                    label_reservations.append(_text_bounds_at_baseline(dso_name, text_font, label_pos))

    # ---- Star/planet highlight ----
    if draw_hover_info and highlighted_object:
        obj, pos = highlighted_object
        painter.setPen(QPen(text_color, 2))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawEllipse(pos, 10, 10)

        # PlanetBody(dataclass) or star(dict)
        if not isinstance(obj, PlanetBody):
            name = getattr(obj, "name", "") if hasattr(obj, "name") else obj.get("name", "")
            label_pos = QPointF(pos.x() + 15, pos.y() - 15)
            if label_candidates is not None:
                label_candidates.append(
                    {
                        "text": str(name),
                        "pos": label_pos,
                        "text_color": text_color,
                        "outline_color": outline_color,
                        "outline_width": outline_width,
                        "priority": 10,
                    }
                )
            else:
                draw_outlined_text(
                    painter,
                    name,
                    label_pos,
                    text_font,
                    text_color,
                    outline_color,
                    outline_width=outline_width,
                )
                if label_reservations is not None:
                    label_reservations.append(_text_bounds_at_baseline(str(name), text_font, label_pos))
