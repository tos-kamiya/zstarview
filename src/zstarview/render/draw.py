import math
import re
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Tuple

import hashlib
import logging
import numpy as np

from PySide6.QtCore import QPoint, QPointF, QRectF, Qt
from PySide6.QtGui import (
    QFont,
    QImage,
    QColor,
    QPainter,
    QPen,
    QPolygonF,
)

from ..paths import (
    FIELD_OF_VIEW_DEG,
    TEXT_STYLES_BY_PRESET,
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
from ..types import ScreenGeometry, CelestialData, ViewerData, CelestialObject
from ..astro import altaz_to_normalized_xy, resolve_star_names, resolve_star_source_ids, is_in_fov, is_in_fov_vectorized
from ..asterisms import ASTERISMS, ASTERISM_REQUIRED_SOURCE_IDS, pick_rotating_asterism
from .geometry import (
    _altaz_to_normalized_xy_vectorized,
    _normalized_to_screen_xy_vectorized,
    normalized_to_screen_xy,
)
from .guides import (
    _clip_polyline_to_radius,
    _great_circle_altaz_points,
    draw_gauge_cross,
    split_by_gaps,
)
from .photometry import compute_flare_profile as _compute_flare_profile, flare_strength_from_vmag as _flare_strength_from_vmag, bv_to_rgb_vectorized
from .text import (
    get_text_style,
    _text_bounds_at_baseline,
    draw_outlined_text,
)

logger = logging.getLogger(__name__)

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
