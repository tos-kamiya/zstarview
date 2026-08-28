from __future__ import annotations

import math

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..precipitation import (
    OBSERVER_PRECIPITATION_MARKER_SCALE,
    PRECIPITATION_COLUMN_COLOR_RGB,
    PRECIPITATION_NEAR_STREAK_HEIGHT_DEG,
    PRECIPITATION_TILE_ALPHA_SCALE,
    PRECIPITATION_TILE_LINE_WIDTH,
    PRECIPITATION_TILE_LINE_WIDTH_RATIO,
    ObserverPrecipitationMarker,
    PrecipitationRenderItem,
    precipitation_distance_opacity_factor,
    precipitation_streak_count,
    precipitation_streak_height_deg,
)
from ..types import ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy


def draw_precipitation_columns(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    columns: list[PrecipitationRenderItem] | None,
    *,
    opacity: float,
    line_width_scale: float = 1.0,
) -> None:
    if not columns or opacity <= 0.0:
        return
    view_center = viewer.view_center
    painter.save()
    for column in columns:
        if isinstance(column, ObserverPrecipitationMarker):
            _draw_observer_precipitation_marker(
                painter,
                geometry,
                column,
                opacity=opacity,
                edge_fov_deg=float(viewer.edge_fov_deg),
                line_width_scale=line_width_scale,
            )
            continue
        if not (
            is_in_fov(
                column.base_alt_deg,
                column.base_az_deg,
                view_center,
                fov_deg=float(viewer.content_fov_deg),
            )
            or is_in_fov(
                column.top_alt_deg,
                column.top_az_deg,
                view_center,
                fov_deg=float(viewer.content_fov_deg),
            )
        ):
            continue
        base_nx, base_ny = altaz_to_normalized_xy(
            column.base_alt_deg,
            column.base_az_deg,
            view_center,
            edge_fov_deg=float(viewer.edge_fov_deg),
        )
        base_x, base_y = normalized_to_screen_xy(base_nx, base_ny, geometry)
        distance_opacity = precipitation_distance_opacity_factor(column.distance_km)
        alpha = int(
            round(255.0 * min(1.0, max(0.0, opacity * distance_opacity)))
        )
        solid_pen = QPen(
            QColor(
                *PRECIPITATION_COLUMN_COLOR_RGB,
                _precipitation_line_alpha(alpha),
            )
        )
        tile_side_px = _precipitation_tile_side_px(
            geometry,
            precipitation_streak_height_deg(column.distance_km),
            float(viewer.edge_fov_deg),
        )
        screen_up_x, screen_up_y = _precipitation_screen_up_vector(
            column.base_alt_deg,
            column.base_az_deg,
            view_center,
            float(viewer.edge_fov_deg),
            geometry,
            float(base_x),
            float(base_y),
        )
        tile_center_x = float(base_x) + screen_up_x * tile_side_px * 0.5
        tile_center_y = float(base_y) + screen_up_y * tile_side_px * 0.5
        screen_rotation_deg = math.degrees(
            math.atan2(screen_up_y, screen_up_x)
        ) + 90.0
        line_width = _precipitation_tile_line_width(tile_side_px, line_width_scale)
        solid_pen.setWidthF(line_width)
        solid_pen.setCosmetic(True)
        solid_pen.setCapStyle(Qt.PenCapStyle.FlatCap)
        streak_count = precipitation_streak_count(column.rate_mm_h)
        painter.save()
        painter.setPen(solid_pen)
        for start, end in _precipitation_tile_lines(
            tile_center_x,
            tile_center_y,
            tile_side_px,
            streak_count,
            line_width_scale,
            rotation_deg=screen_rotation_deg,
        ):
            painter.drawLine(start, end)
        painter.restore()
    painter.restore()


def _draw_observer_precipitation_marker(
    painter: QPainter,
    geometry: ScreenGeometry,
    marker: ObserverPrecipitationMarker,
    *,
    opacity: float,
    edge_fov_deg: float,
    line_width_scale: float,
) -> None:
    marker_scale = OBSERVER_PRECIPITATION_MARKER_SCALE
    alpha = int(round(255.0 * min(1.0, max(0.0, opacity))))
    solid_pen = QPen(
        QColor(
            *PRECIPITATION_COLUMN_COLOR_RGB,
            _precipitation_line_alpha(alpha),
        )
    )
    streak_count = precipitation_streak_count(marker.rate_mm_h)
    tile_side_px = (
        float(geometry.radius)
        * PRECIPITATION_NEAR_STREAK_HEIGHT_DEG
        / max(1.0e-6, float(edge_fov_deg))
        * marker_scale
    )
    line_width = _precipitation_tile_line_width(tile_side_px, line_width_scale)
    solid_pen.setWidthF(line_width)
    solid_pen.setCosmetic(True)
    solid_pen.setCapStyle(Qt.PenCapStyle.FlatCap)
    center_x, center_y = geometry.center
    painter.save()
    painter.setPen(solid_pen)
    for start, end in _precipitation_tile_lines(
        float(center_x), float(center_y), tile_side_px, streak_count, line_width_scale
    ):
        painter.drawLine(start, end)
    painter.restore()


def _precipitation_tile_side_px(
    geometry: ScreenGeometry,
    height_deg: float,
    edge_fov_deg: float,
) -> float:
    return max(
        1.0,
        float(geometry.radius)
        * float(height_deg)
        / max(1.0e-6, float(edge_fov_deg)),
    )


def _precipitation_tile_line_width(side_px: float, line_width_scale: float) -> float:
    scale = max(0.0, float(line_width_scale))
    return max(
        0.1,
        PRECIPITATION_TILE_LINE_WIDTH * scale,
        float(side_px) * PRECIPITATION_TILE_LINE_WIDTH_RATIO * scale,
    )


def _precipitation_line_alpha(base_alpha: int) -> int:
    return min(255, int(round(float(base_alpha) * 0.7 * PRECIPITATION_TILE_ALPHA_SCALE)))


def _precipitation_line_offsets(line_count: int) -> tuple[float, ...]:
    if line_count <= 0:
        return ()
    if line_count % 2:
        center = line_count // 2
        return tuple(float(index - center) for index in range(line_count))
    center = line_count / 2.0
    return tuple(float(index - center + 0.5) for index in range(line_count))


def _clip_precipitation_line_to_vertical_band(
    center_x: float,
    center_y: float,
    side_px: float,
    normal_offset_px: float,
) -> tuple[QPointF, QPointF] | None:
    half_side = max(0.5, float(side_px)) * 0.5
    direction_x = 1.0 / math.sqrt(2.0)
    direction_y = -direction_x
    normal_x = direction_x
    normal_y = -direction_y
    point_x = float(center_x) + normal_x * float(normal_offset_px)
    point_y = float(center_y) + normal_y * float(normal_offset_px)
    if abs(direction_y) < 1.0e-12:
        return None
    first = (center_y - half_side - point_y) / direction_y
    second = (center_y + half_side - point_y) / direction_y
    lower = min(first, second)
    upper = max(first, second)
    return (
        QPointF(point_x + direction_x * lower, point_y + direction_y * lower),
        QPointF(point_x + direction_x * upper, point_y + direction_y * upper),
    )


def _precipitation_tile_lines(
    center_x: float,
    center_y: float,
    side_px: float,
    line_count: int,
    line_width_scale: float,
    *,
    rotation_deg: float = 0.0,
) -> tuple[tuple[QPointF, QPointF], ...]:
    offsets = _precipitation_line_offsets(line_count)
    if not offsets:
        return ()
    max_offset = max(abs(offset) for offset in offsets)
    nominal_spacing = _precipitation_tile_line_spacing(
        side_px, line_width_scale
    )
    max_normal_offset = max(0.0, float(side_px)) / math.sqrt(2.0) * 0.8
    spacing = (
        nominal_spacing
        if max_offset <= 0.0
        else min(nominal_spacing, max_normal_offset / max_offset)
    )
    lines = []
    for offset in offsets:
        line = _clip_precipitation_line_to_vertical_band(
            center_x,
            center_y,
            side_px,
            offset * spacing,
        )
        if line is not None:
            lines.append(
                _rotate_precipitation_line(line, center_x, center_y, rotation_deg)
            )
    return tuple(lines)


def _rotate_precipitation_line(
    line: tuple[QPointF, QPointF],
    center_x: float,
    center_y: float,
    rotation_deg: float,
) -> tuple[QPointF, QPointF]:
    if abs(float(rotation_deg)) <= 1.0e-12:
        return line
    angle_rad = math.radians(float(rotation_deg))
    cos_angle = math.cos(angle_rad)
    sin_angle = math.sin(angle_rad)

    def rotate(point: QPointF) -> QPointF:
        dx = point.x() - center_x
        dy = point.y() - center_y
        return QPointF(
            center_x + cos_angle * dx - sin_angle * dy,
            center_y + sin_angle * dx + cos_angle * dy,
        )

    return rotate(line[0]), rotate(line[1])


def _precipitation_screen_up_vector(
    base_alt_deg: float,
    base_az_deg: float,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    geometry: ScreenGeometry,
    base_x: float,
    base_y: float,
) -> tuple[float, float]:
    step_deg = 0.1
    if float(base_alt_deg) < 90.0 - step_deg:
        probe_alt_deg = float(base_alt_deg) + step_deg
        probe_sign = 1.0
    else:
        probe_alt_deg = float(base_alt_deg) - step_deg
        probe_sign = -1.0
    probe_nx, probe_ny = altaz_to_normalized_xy(
        probe_alt_deg,
        base_az_deg,
        view_center,
        edge_fov_deg=edge_fov_deg,
    )
    probe_x, probe_y = normalized_to_screen_xy(probe_nx, probe_ny, geometry)
    up_x = (probe_x - base_x) * probe_sign
    up_y = (probe_y - base_y) * probe_sign
    length = math.hypot(up_x, up_y)
    if length <= 1.0e-12:
        return 0.0, -1.0
    return up_x / length, up_y / length


def _precipitation_screen_up_rotation_deg(
    base_alt_deg: float,
    base_az_deg: float,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    geometry: ScreenGeometry,
    base_x: float,
    base_y: float,
) -> float:
    up_x, up_y = _precipitation_screen_up_vector(
        base_alt_deg,
        base_az_deg,
        view_center,
        edge_fov_deg,
        geometry,
        base_x,
        base_y,
    )
    return math.degrees(math.atan2(up_y, up_x)) + 90.0


def _precipitation_tile_line_spacing(
    side_px: float,
    line_width_scale: float,
) -> float:
    return max(0.0, float(side_px)) * 0.25 * max(
        0.0, float(line_width_scale)
    )
