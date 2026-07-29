import math
from collections.abc import Callable
from dataclasses import dataclass

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
    ThemeStyle,
)
from ..types import CelestialData, ScreenGeometry, ViewerData, ViewProjection
from .geometry import normalized_to_screen_xy
from .text import (
    ResolvedTextStyle,
    _clamp_baseline_pos_to_viewport,
    draw_outlined_text,
    resolve_label_text_style,
)

REFERENCE_LINE_OUTER_WIDTH = 1.1
REFERENCE_LINE_MID_WIDTH = 0.75
REFERENCE_LINE_FG_WIDTH = 0.55
REFERENCE_LINE_OUTER_ALPHA = 18
REFERENCE_LINE_MID_ALPHA = 30
_DIRECTION_MARKER_HIT_RADIUS_PX = 9.0
GRID_LINE_WIDTH = 0.51
GRID_LINE_ALPHA = 190
GRID_MAJOR_LINE_WIDTH_SCALE = 1.0
GRID_MINOR_CROSS_HALF_LEN_SCALE = 0.006
GRID_MINOR_CROSS_WIDTH = 0.51
GRID_MINOR_CROSS_SAMPLE_DEG = 1.5
GRID_ALTITUDE_SAMPLES = 73
GRID_ALTITUDE_BANDS = (0.0,)
GRID_FINE_ALTITUDE_BANDS = tuple(
    float(value)
    for value in range(-80, 81, 10)
    if value != 0
)
GRID_PARALLEL_ALTITUDES = GRID_ALTITUDE_BANDS + GRID_FINE_ALTITUDE_BANDS
GRID_FINE_AZIMUTHS = tuple(
    float(value)
    for value in range(0, 360, 10)
)
GRID_PARALLEL_AZ_SAMPLES = 145
AXIS_MARKER_HALF_SIZE_PX = 7


def _is_major_grid_step(value: float) -> bool:
    return math.isclose(abs(float(value)) % 30.0, 0.0, abs_tol=1.0e-6)


def _direction_grid_minor_cross_half_len(surface_size: tuple[int, int]) -> float:
    width, height = surface_size
    return max(1.5, float(min(max(1, int(width)), max(1, int(height)))) * GRID_MINOR_CROSS_HALF_LEN_SCALE)


@dataclass(frozen=True, slots=True)
class DirectionMarkerHover:
    label: str
    az_deg: float
    screen_pos: QPointF


def _content_fov_deg_from_viewer(viewer_data: ViewerData) -> float:
    return float(viewer_data.content_fov_deg)


def _draw_cross_marker_at_altaz(
    painter: QPainter,
    geometry: ScreenGeometry,
    projection: ViewProjection,
    *,
    alt_deg: float,
    az_deg: float,
    color: tuple[int, int, int],
) -> None:
    if not is_in_fov(alt_deg, az_deg, tuple(float(value) for value in projection.view_center), fov_deg=float(projection.content_fov_deg)):
        return
    nx, ny = altaz_to_normalized_xy(
        alt_deg,
        az_deg,
        tuple(float(value) for value in projection.view_center),
        edge_fov_deg=float(projection.edge_fov_deg),
    )
    x, y = normalized_to_screen_xy(nx, ny, geometry)
    s = AXIS_MARKER_HALF_SIZE_PX
    painter.setPen(QPen(QColor(*color), 1))
    painter.drawLine(QPointF(x - s, y - s), QPointF(x + s, y + s))
    painter.drawLine(QPointF(x - s, y + s), QPointF(x + s, y - s))


def resolve_direction_marker_hover(
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    mouse_pos: QPoint | None,
) -> DirectionMarkerHover | None:
    if mouse_pos is None:
        return None

    view_center = tuple(float(value) for value in viewer_data.view_center)
    edge_fov_deg = float(viewer_data.edge_fov_deg)
    content_fov_deg = float(viewer_data.content_fov_deg)
    mouse_x = float(mouse_pos.x())
    mouse_y = float(mouse_pos.y())
    best: DirectionMarkerHover | None = None
    best_dist_sq = _DIRECTION_MARKER_HIT_RADIUS_PX * _DIRECTION_MARKER_HIT_RADIUS_PX
    marker_alt = 0.0

    for label, az in DIRECTIONS.items():
        if not is_in_fov(marker_alt, az, view_center, fov_deg=content_fov_deg):
            continue
        marker_nx, marker_ny = altaz_to_normalized_xy(
            marker_alt,
            az,
            view_center,
            edge_fov_deg=edge_fov_deg,
        )
        screen_x, screen_y = normalized_to_screen_xy(marker_nx, marker_ny, geometry)
        dx = float(screen_x) - mouse_x
        dy = float(screen_y) - mouse_y
        dist_sq = (dx * dx) + (dy * dy)
        if dist_sq <= best_dist_sq:
            best = DirectionMarkerHover(
                label=label,
                az_deg=float(az),
                screen_pos=QPointF(float(screen_x), float(screen_y)),
            )
            best_dist_sq = dist_sq

    return best


def _project_reference_altaz_point(
    alt_deg: float,
    az_deg: float,
    *,
    projection: ViewProjection,
    altaz_to_normalized_xy_func: Callable[..., tuple[float, float]] | None,
) -> tuple[float, float]:
    project_xy = (
        altaz_to_normalized_xy
        if altaz_to_normalized_xy_func is None
        else altaz_to_normalized_xy_func
    )
    try:
        nx, ny = project_xy(
            float(alt_deg),
            float(az_deg),
            tuple(float(value) for value in projection.view_center),
            edge_fov_deg=float(projection.edge_fov_deg),
        )
    except TypeError:
        nx, ny = project_xy(float(alt_deg), float(az_deg), tuple(float(value) for value in projection.view_center))
    return float(nx), float(ny)


def _altaz_to_neu_unit(alt_deg: float, az_deg: float) -> np.ndarray:
    """Convert Alt/Az to a unit vector in local North-East-Up coordinates."""
    alt = math.radians(float(alt_deg))
    az = math.radians(float(az_deg))
    c_alt = math.cos(alt)
    return np.array(
        [
            c_alt * math.cos(az),  # north
            c_alt * math.sin(az),  # east
            math.sin(alt),  # up
        ],
        dtype=float,
    )


def _neu_unit_to_altaz(vec: np.ndarray) -> tuple[float, float]:
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
) -> list[tuple[float, float]]:
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

    out: list[tuple[float, float]] = []
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


def split_by_gaps(points: list[tuple[float, float]]) -> list[list[tuple[float, float]]]:
    """
    Split a polyline by large gaps to avoid drawing long, straight lines
    across the screen when a celestial path wraps around.

    Args:
        points: A list of (x, y) tuples representing the polyline.

    Returns:
        A list of polyline fragments, where each fragment is a list of points.
    """

    def dist(p1: tuple[float, float], p2: tuple[float, float]) -> float:
        return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

    fragments: list[list[tuple[float, float]]] = [[]]
    for p in points:
        if not fragments[-1] or dist(p, fragments[-1][-1]) < 0.2:
            fragments[-1].append(p)
        else:
            fragments.append([p])
    return fragments


def _clip_polyline_to_radius(
    points: list[tuple[float, float]],
    max_radius: float,
) -> list[list[tuple[float, float]]]:
    """Clip a normalized polyline to a circle centered at the origin."""
    if not points:
        return []

    radius_sq = max_radius * max_radius

    def _inside(point: tuple[float, float]) -> bool:
        return (point[0] * point[0]) + (point[1] * point[1]) <= radius_sq

    def _intersections(
        start: tuple[float, float],
        end: tuple[float, float],
    ) -> list[tuple[float, float, float]]:
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
        hits: list[tuple[float, float, float]] = []
        for t in ts:
            if 0.0 <= t <= 1.0:
                hits.append((t, x0 + (t * dx), y0 + (t * dy)))
        hits.sort(key=lambda hit: hit[0])
        unique: list[tuple[float, float, float]] = []
        for hit in hits:
            if unique and abs(hit[0] - unique[-1][0]) < 1.0e-9:
                continue
            unique.append(hit)
        return unique

    fragments: list[list[tuple[float, float]]] = []
    current: list[tuple[float, float]] = [points[0]] if _inside(points[0]) else []

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
    viewer_data: ViewerData,
    celestial_data: CelestialData,
    *,
    altaz_to_normalized_xy_func: Callable[..., tuple[float, float]] | None = None,
    theme: ThemeStyle | None = None,
) -> None:
    """
    Draw celestial reference lines like the equator, ecliptic, and horizon.

    Args:
        painter: The QPainter to use for drawing.
        geometry: The screen geometry for coordinate conversion.
        celestial_data: The data containing the points for the reference lines.
    """
    effective_fov_deg = _content_fov_deg_from_viewer(viewer_data)
    projection = ViewProjection(
        view_center=tuple(float(value) for value in viewer_data.view_center),
        edge_fov_deg=float(viewer_data.edge_fov_deg),
        content_fov_deg=effective_fov_deg,
    )
    painter.save()
    guide_style = theme.guide_style if theme is not None else None

    def _make_reference_pen(color: tuple[int, int, int], width: float, alpha: int, style: Qt.PenStyle | None = None) -> QPen:
        pen_color = QColor(*color)
        pen_color.setAlpha(alpha)
        pen = QPen(pen_color, width, style) if style is not None else QPen(pen_color, width)
        pen.setCosmetic(True)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        return pen

    def _draw_reference_line(
        altaz_points: list[tuple[float, float]],
        color: tuple[int, int, int],
        dash_pattern: list[int] | None,
        *,
        width_scale: float = 1.0,
        fg_width: float | None = None,
        simple_dash_pattern: tuple[float, ...] | None = None,
    ) -> None:
        width_scale = max(1.0, float(width_scale))
        if len(altaz_points) == 1:
            _project_reference_altaz_point(
                float(altaz_points[0][0]),
                float(altaz_points[0][1]),
                projection=projection,
                altaz_to_normalized_xy_func=altaz_to_normalized_xy_func,
            )
        projected_points: list[tuple[float, float]] = []
        for alt_deg, az_deg in altaz_points:
            nx, ny = _project_reference_altaz_point(
                float(alt_deg),
                float(az_deg),
                projection=projection,
                altaz_to_normalized_xy_func=altaz_to_normalized_xy_func,
            )
            projected_points.append((nx, ny))
        fragments = split_by_gaps(projected_points)
        for frag in fragments:
            if len(frag) < 2:
                continue
            clipped_frags = _clip_polyline_to_radius(frag, effective_fov_deg / max(1.0e-6, float(viewer_data.edge_fov_deg)))
            for clipped_frag in clipped_frags:
                if len(clipped_frag) < 2:
                    continue
                pts = [
                    QPointF(*normalized_to_screen_xy(x, y, geometry))
                    for x, y in clipped_frag
                ]
                poly = QPolygonF(pts)

                if guide_style is not None and guide_style.simple_reference_lines:
                    simple = _make_reference_pen(
                        guide_style.reference_rgb,
                        guide_style.reference_width,
                        guide_style.reference_alpha,
                        Qt.PenStyle.SolidLine,
                    )
                    if simple_dash_pattern:
                        simple.setDashPattern(list(simple_dash_pattern))
                    painter.setPen(simple)
                    painter.drawPolyline(poly)
                else:
                    outer = _make_reference_pen(
                        color,
                        REFERENCE_LINE_OUTER_WIDTH * width_scale,
                        REFERENCE_LINE_OUTER_ALPHA,
                        Qt.PenStyle.SolidLine,
                    )
                    painter.setPen(outer)
                    painter.drawPolyline(poly)

                    mid = _make_reference_pen(
                        color,
                        REFERENCE_LINE_MID_WIDTH * width_scale,
                        REFERENCE_LINE_MID_ALPHA,
                        Qt.PenStyle.SolidLine,
                    )
                    painter.setPen(mid)
                    painter.drawPolyline(poly)

                    fg_pen_width = REFERENCE_LINE_FG_WIDTH if fg_width is None else float(fg_width)
                    fg = _make_reference_pen(color, fg_pen_width * width_scale, 255)
                    if dash_pattern:
                        fg.setDashPattern(dash_pattern)
                    painter.setPen(fg)
                    painter.drawPolyline(poly)

    # Keep the ecliptic dash cadence visible while restoring the equator as a longer dash.
    _draw_reference_line(
        celestial_data.celestial_equator_points,
        CELESTIAL_EQUATOR_COLOR if guide_style is None else guide_style.equator_rgb,
        [16, 6],
        width_scale=1.0,
        fg_width=GRID_LINE_WIDTH,
    )
    _draw_reference_line(
        celestial_data.ecliptic_points,
        ECLIPTIC_COLOR if guide_style is None else guide_style.ecliptic_rgb,
        [4, 6],
        width_scale=1.14,
        simple_dash_pattern=(
            None
            if guide_style is None
            else guide_style.ecliptic_dash_pattern
        ),
    )
    _draw_reference_line(
        celestial_data.horizon_points,
        HORIZON_LINE_COLOR if guide_style is None else guide_style.horizon_rgb,
        [10, 1],
    )
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
    viewer_data: ViewerData,
    *,
    theme: ThemeStyle | None = None,
) -> None:
    """
    Draws markers at zenith and nadir.

    Args:
        painter: The QPainter to use for drawing.
        geometry: The screen geometry for coordinate conversion.
        view_center: The current view center (altitude, azimuth).
    """
    guide_style = theme.guide_style if theme is not None else None
    # Match the direction-grid color so the zenith/nadir markers stay visually
    # aligned with the compass overlay in every theme.
    projection = ViewProjection(
        view_center=tuple(float(value) for value in viewer_data.view_center),
        edge_fov_deg=float(viewer_data.edge_fov_deg),
        content_fov_deg=float(viewer_data.content_fov_deg),
    )
    view_center = tuple(float(value) for value in projection.view_center)
    az_ref = view_center[1]
    for alt in (90.0, -90.0):
        _draw_cross_marker_at_altaz(
            painter,
            geometry,
            projection,
            alt_deg=alt,
            az_deg=az_ref,
            color=(
                HORIZON_LINE_COLOR
                if guide_style is None
                else (guide_style.grid_rgb or guide_style.horizon_rgb)
            ),
        )


def draw_celestial_pole_markers(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    *,
    theme: ThemeStyle | None = None,
) -> None:
    """Draw X markers at the north and south celestial poles."""
    guide_style = theme.guide_style if theme is not None else None
    lat_deg = float(viewer_data.lat_deg)
    projection = ViewProjection(
        view_center=tuple(float(value) for value in viewer_data.view_center),
        edge_fov_deg=float(viewer_data.edge_fov_deg),
        content_fov_deg=float(viewer_data.content_fov_deg),
    )
    pole_specs = (
        (lat_deg, 0.0),
        (-lat_deg, 180.0),
    )
    for alt_deg, az_deg in pole_specs:
        _draw_cross_marker_at_altaz(
            painter,
            geometry,
            projection,
            alt_deg=alt_deg,
            az_deg=az_deg,
            color=CELESTIAL_EQUATOR_COLOR if guide_style is None else guide_style.equator_rgb,
        )


def draw_direction_labels(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    text_font: QFont,
    mouse_pos: QPoint | None = None,
    *,
    theme: ThemeStyle,
) -> None:
    """
    Draw compass direction labels and horizon markers on the horizon.

    Args:
        painter: The QPainter for drawing.
        geometry: The screen geometry for coordinate conversion.
        view_center: The current view center to determine which labels are visible.
        text_font: The QFont to use for the labels.
    """
    view_center = tuple(float(value) for value in viewer_data.view_center)
    edge_fov_deg = float(viewer_data.edge_fov_deg)
    content_fov_deg = float(viewer_data.content_fov_deg)
    direction_label_font = QFont(text_font)
    direction_label_font.setWeight(QFont.Weight.DemiBold)
    label_style = resolve_label_text_style(theme, direction_label_font)
    label_style = ResolvedTextStyle(
        font=label_style.font,
        text_color=QColor(
            *(theme.guide_style.label_rgb or theme.guide_style.horizon_rgb)
        ),
        outline_color=QColor(0, 0, 0, 0),
        outline_width=0.0,
    )
    marker_color = QColor(*theme.guide_style.horizon_rgb)
    marker_pen = QPen(marker_color, theme.guide_style.marker_width)
    marker_pen.setCosmetic(True)
    marker_half_len_px = 6.0
    marker_hit_radius_px = 4.0
    tangent_probe_deg = 0.6
    label_outward_offset_px = 4.0
    hide_near_mouse_px = 28.0
    mouse_x = float(mouse_pos.x()) if mouse_pos is not None else None
    mouse_y = float(mouse_pos.y()) if mouse_pos is not None else None
    painter.setFont(direction_label_font)
    fm = QFontMetrics(direction_label_font)
    marker_alt = 0.0
    label_alt = -2.0
    for label, az in DIRECTIONS.items():
        if not is_in_fov(marker_alt, az, view_center, fov_deg=content_fov_deg):
            continue
        marker_nx, marker_ny = altaz_to_normalized_xy(
            marker_alt,
            az,
            view_center,
            edge_fov_deg=edge_fov_deg,
        )
        pos = QPointF(*normalized_to_screen_xy(marker_nx, marker_ny, geometry))
        label_nx, label_ny = altaz_to_normalized_xy(
            label_alt,
            az,
            view_center,
            edge_fov_deg=edge_fov_deg,
        )
        label_pos = QPointF(*normalized_to_screen_xy(label_nx, label_ny, geometry))

        az_prev = (az - tangent_probe_deg + 360.0) % 360.0
        az_next = (az + tangent_probe_deg) % 360.0
        p_prev = QPointF(
            *normalized_to_screen_xy(
                *altaz_to_normalized_xy(
                    marker_alt,
                    az_prev,
                    view_center,
                    edge_fov_deg=edge_fov_deg,
                ),
                geometry,
            )
        )
        p_next = QPointF(
            *normalized_to_screen_xy(
                *altaz_to_normalized_xy(
                    marker_alt,
                    az_next,
                    view_center,
                    edge_fov_deg=edge_fov_deg,
                ),
                geometry,
            )
        )
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

        label_pos = _clamp_baseline_pos_to_viewport(
            label,
            direction_label_font,
            label_pos,
            QRectF(painter.viewport()),
        )

        draw_outlined_text(
            painter,
            label,
            label_pos,
            direction_label_font,
            style=label_style,
        )


def _draw_direction_polyline(
    painter: QPainter,
    points: list[tuple[float, float]],
    geometry: ScreenGeometry,
    *,
    width: float,
    alpha: int,
    color_rgb: tuple[int, int, int] = HORIZON_LINE_COLOR,
) -> None:
    fragments = split_by_gaps(points)
    if not fragments:
        return

    def _make_pen(width: float, alpha: int) -> QPen:
        color = QColor(*color_rgb)
        color.setAlpha(alpha)
        pen = QPen(color, width)
        pen.setCosmetic(True)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        return pen

    for frag in fragments:
        if len(frag) < 2:
            continue
        poly = QPolygonF([QPointF(*normalized_to_screen_xy(x, y, geometry)) for x, y in frag])
        painter.setPen(_make_pen(width, alpha))
        painter.drawPolyline(poly)


def _draw_direction_cross_marker(
    painter: QPainter,
    alt_deg: float,
    az_deg: float,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    *,
    edge_fov_deg: float,
    width: float,
    alpha: int,
    half_len: float,
    color_rgb: tuple[int, int, int] = HORIZON_LINE_COLOR,
) -> None:
    if not is_in_fov(float(alt_deg), float(az_deg), view_center, fov_deg=FIELD_OF_VIEW_DEG):
        return

    delta = float(GRID_MINOR_CROSS_SAMPLE_DEG)
    alt_lo = max(-90.0, float(alt_deg) - delta)
    alt_hi = min(90.0, float(alt_deg) + delta)
    az_lo = (float(az_deg) - delta) % 360.0
    az_hi = (float(az_deg) + delta) % 360.0

    alt_lo_nx, alt_lo_ny = altaz_to_normalized_xy(
        alt_lo,
        float(az_deg),
        view_center,
        edge_fov_deg=edge_fov_deg,
    )
    alt_hi_nx, alt_hi_ny = altaz_to_normalized_xy(
        alt_hi,
        float(az_deg),
        view_center,
        edge_fov_deg=edge_fov_deg,
    )
    az_lo_nx, az_lo_ny = altaz_to_normalized_xy(
        float(alt_deg),
        az_lo,
        view_center,
        edge_fov_deg=edge_fov_deg,
    )
    az_hi_nx, az_hi_ny = altaz_to_normalized_xy(
        float(alt_deg),
        az_hi,
        view_center,
        edge_fov_deg=edge_fov_deg,
    )

    alt_vec_x = alt_hi_nx - alt_lo_nx
    alt_vec_y = alt_hi_ny - alt_lo_ny
    az_vec_x = az_hi_nx - az_lo_nx
    az_vec_y = az_hi_ny - az_lo_ny
    alt_norm = math.hypot(alt_vec_x, alt_vec_y)
    az_norm = math.hypot(az_vec_x, az_vec_y)
    if alt_norm <= 1.0e-6 or az_norm <= 1.0e-6:
        return

    alt_unit_x = alt_vec_x / alt_norm
    alt_unit_y = alt_vec_y / alt_norm
    az_unit_x = az_vec_x / az_norm
    az_unit_y = az_vec_y / az_norm

    color = QColor(*color_rgb)
    color.setAlpha(alpha)
    pen = QPen(color, width)
    pen.setCosmetic(True)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    painter.setPen(pen)
    cx, cy = normalized_to_screen_xy(
        *altaz_to_normalized_xy(
            float(alt_deg),
            float(az_deg),
            view_center,
            edge_fov_deg=edge_fov_deg,
        ),
        geometry,
    )
    painter.drawLine(
        QPointF(cx - (alt_unit_x * half_len), cy - (alt_unit_y * half_len)),
        QPointF(cx + (alt_unit_x * half_len), cy + (alt_unit_y * half_len)),
    )
    painter.drawLine(
        QPointF(cx - (az_unit_x * half_len), cy - (az_unit_y * half_len)),
        QPointF(cx + (az_unit_x * half_len), cy + (az_unit_y * half_len)),
    )


def draw_direction_grid_overlay(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    surface_size: tuple[int, int],
    *,
    theme: ThemeStyle | None = None,
) -> None:
    view_center = tuple(float(value) for value in viewer_data.view_center)
    edge_fov_deg = float(viewer_data.edge_fov_deg)
    content_fov_deg = float(viewer_data.content_fov_deg)
    grid_color = (
        HORIZON_LINE_COLOR
        if theme is None
        else (theme.guide_style.grid_rgb or theme.guide_style.horizon_rgb)
    )
    grid_width = GRID_LINE_WIDTH if theme is None else theme.guide_style.grid_width
    grid_minor_width = (
        GRID_MINOR_CROSS_WIDTH
        if theme is None
        else theme.guide_style.grid_minor_width
    )
    grid_alpha = GRID_LINE_ALPHA if theme is None else theme.guide_style.grid_alpha
    painter.save()
    try:
        meridian_alt_samples = np.linspace(-90.0, 90.0, GRID_ALTITUDE_SAMPLES)
        parallel_az_samples = np.linspace(
            0.0, 360.0, GRID_PARALLEL_AZ_SAMPLES, endpoint=False
        )

        major_parallel_alts = tuple(
            alt for alt in GRID_PARALLEL_ALTITUDES if _is_major_grid_step(alt)
        )
        minor_parallel_alts = tuple(
            alt for alt in GRID_PARALLEL_ALTITUDES if not _is_major_grid_step(alt)
        )
        major_azimuths = tuple(
            az for az in GRID_FINE_AZIMUTHS if _is_major_grid_step(az)
        )
        minor_azimuths = tuple(
            az for az in GRID_FINE_AZIMUTHS if not _is_major_grid_step(az)
        )

        for alt in major_parallel_alts:
            parallel_points: list[tuple[float, float]] = []
            for az in parallel_az_samples:
                az_norm = float(az) % 360.0
                if not is_in_fov(float(alt), az_norm, view_center, fov_deg=content_fov_deg):
                    continue
                nx, ny = altaz_to_normalized_xy(
                    float(alt),
                    az_norm,
                    view_center,
                    edge_fov_deg=edge_fov_deg,
                )
                parallel_points.append((nx, ny))

            _draw_direction_polyline(
                painter,
                parallel_points,
                geometry,
                width=grid_width
                * (GRID_MAJOR_LINE_WIDTH_SCALE if _is_major_grid_step(alt) else 1.0),
                alpha=grid_alpha,
                color_rgb=grid_color,
            )

        for az in major_azimuths:
            meridian_points = []
            for alt in meridian_alt_samples:
                if not is_in_fov(
                    float(alt), az, view_center, fov_deg=content_fov_deg
                ):
                    continue
                nx, ny = altaz_to_normalized_xy(
                    float(alt),
                    az,
                    view_center,
                    edge_fov_deg=edge_fov_deg,
                )
                meridian_points.append((nx, ny))

            _draw_direction_polyline(
                painter,
                meridian_points,
                geometry,
                width=grid_width
                * (GRID_MAJOR_LINE_WIDTH_SCALE if _is_major_grid_step(az) else 1.0),
                alpha=grid_alpha,
                color_rgb=grid_color,
            )

        for alt in minor_parallel_alts:
            for az in minor_azimuths:
                if not is_in_fov(float(alt), float(az), view_center, fov_deg=content_fov_deg):
                    continue
                nx, ny = altaz_to_normalized_xy(
                    float(alt),
                    float(az),
                    view_center,
                    edge_fov_deg=edge_fov_deg,
                )
                _draw_direction_cross_marker(
                    painter,
                    float(alt),
                    float(az),
                    geometry,
                    view_center,
                    edge_fov_deg=edge_fov_deg,
                    width=grid_minor_width,
                    alpha=grid_alpha,
                    half_len=_direction_grid_minor_cross_half_len(surface_size),
                    color_rgb=grid_color,
                )
    finally:
        painter.restore()
