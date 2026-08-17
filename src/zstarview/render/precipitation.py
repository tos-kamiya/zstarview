from __future__ import annotations

import math

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QPainter, QPen

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..precipitation import (
    OBSERVER_PRECIPITATION_MARKER_SCALE,
    PRECIPITATION_COLUMN_DARK_COLOR_RGB,
    PRECIPITATION_NEAR_STREAK_HEIGHT_DEG,
    PRECIPITATION_TILE_ALPHA_SCALE,
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
        solid_pen = QPen(QColor(*PRECIPITATION_COLUMN_DARK_COLOR_RGB, _precipitation_line_alpha(alpha)))
        tile_side_px = _precipitation_tile_side_px(
            geometry,
            precipitation_streak_height_deg(column.distance_km),
            float(viewer.edge_fov_deg),
        )
        line_width = _precipitation_tile_line_width(tile_side_px, line_width_scale)
        solid_pen.setWidthF(line_width)
        solid_pen.setCosmetic(True)
        solid_pen.setCapStyle(Qt.PenCapStyle.FlatCap)
        streak_count = precipitation_streak_count(column.rate_mm_h)
        painter.save()
        painter.setClipRect(
            QRectF(
                float(base_x) - tile_side_px * 0.5,
                float(base_y) - tile_side_px * 0.5,
                tile_side_px,
                tile_side_px,
            )
        )
        painter.setPen(solid_pen)
        for start, end in _precipitation_tile_lines(
            float(base_x), float(base_y), tile_side_px, streak_count, line_width_scale
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
    solid_pen = QPen(QColor(*PRECIPITATION_COLUMN_DARK_COLOR_RGB, _precipitation_line_alpha(alpha)))
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
    painter.setClipRect(
        QRectF(
            float(center_x) - tile_side_px * 0.5,
            float(center_y) - tile_side_px * 0.5,
            tile_side_px,
            tile_side_px,
        )
    )
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
    return max(
        0.1,
        float(side_px)
        * PRECIPITATION_TILE_LINE_WIDTH_RATIO
        * max(0.0, float(line_width_scale)),
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


def _clip_precipitation_line_to_tile(
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
    lower = -math.inf
    upper = math.inf
    for point, direction, axis_center in (
        (point_x, direction_x, center_x),
        (point_y, direction_y, center_y),
    ):
        if abs(direction) < 1.0e-12:
            if point < axis_center - half_side or point > axis_center + half_side:
                return None
            continue
        first = (axis_center - half_side - point) / direction
        second = (axis_center + half_side - point) / direction
        lower = max(lower, min(first, second))
        upper = min(upper, max(first, second))
    if lower > upper:
        return None
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
        line = _clip_precipitation_line_to_tile(
            center_x,
            center_y,
            side_px,
            offset * spacing,
        )
        if line is not None:
            lines.append(line)
    return tuple(lines)


def _precipitation_tile_line_spacing(
    side_px: float,
    line_width_scale: float,
) -> float:
    return max(0.0, float(side_px)) * 0.25 * max(
        0.0, float(line_width_scale)
    )
