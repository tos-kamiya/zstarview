import math
from typing import Callable, List, Tuple

import numpy as np

from PySide6.QtCore import QPoint, QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import (
    CELESTIAL_EQUATOR_COLOR,
    DIRECTIONS,
    ECLIPTIC_COLOR,
    FIELD_OF_VIEW_DEG,
    HORIZON_LINE_COLOR,
    TEXT_STYLES_BY_PRESET,
)
from ..types import CelestialData, ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy
from .text import draw_outlined_text


def _content_fov_deg_from_viewer(viewer_data: ViewerData) -> float:
    return float(viewer_data.content_fov_deg)


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
    is_in_fov_func: Callable[..., bool] = is_in_fov,
    altaz_to_normalized_xy_func: Callable[[float, float, Tuple[float, float]], Tuple[float, float]] = altaz_to_normalized_xy,
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
            if not is_in_fov_func(float(alt), float(az), viewer_data.view_center, fov_deg=effective_fov_deg):
                continue
            nx, ny = altaz_to_normalized_xy_func(float(alt), float(az), viewer_data.view_center)
            points.append((nx, ny))
        for frag in split_by_gaps(points):
            if len(frag) < 2:
                continue
            pts = [QPointF(*normalized_to_screen_xy(nx, ny, geometry)) for nx, ny in frag]
            poly = QPolygonF(pts)

            base_color = QColor(*color)
            base_color.setAlpha(70)
            base = QPen(base_color, width + 2, Qt.PenStyle.SolidLine)
            base.setCosmetic(True)
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


def draw_direction_labels(
    painter: QPainter,
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    text_font: QFont,
    mouse_pos: QPoint | None = None,
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
    text_color = QColor(*HORIZON_LINE_COLOR)
    outline_color = QColor.fromRgbF(0, 0, 0, 0.3)
    marker_color = QColor(*HORIZON_LINE_COLOR)
    marker_pen = QPen(marker_color, 1.6)
    marker_pen.setCosmetic(True)
    marker_half_len_px = 6.0
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
        nxp, nyp = -uy, ux
        painter.save()
        painter.setPen(marker_pen)
        painter.drawLine(
            QPointF(pos.x() - nxp * marker_half_len_px, pos.y() - nyp * marker_half_len_px),
            QPointF(pos.x() + nxp * marker_half_len_px, pos.y() + nyp * marker_half_len_px),
        )
        painter.restore()

        label_pos = QPointF(pos)
        bounds = fm.tightBoundingRect(label)
        label_rect = QRectF(
            label_pos.x() + bounds.x(),
            label_pos.y() + bounds.y(),
            bounds.width(),
            bounds.height(),
        )
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
