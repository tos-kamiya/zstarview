from __future__ import annotations

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..precipitation import (
    OBSERVER_PRECIPITATION_MARKER_SCALE,
    PRECIPITATION_COLUMN_DARK_COLOR_RGB,
    PRECIPITATION_COLUMN_COLOR_RGB,
    PRECIPITATION_NEAR_STREAK_HEIGHT_DEG,
    ObserverPrecipitationMarker,
    PrecipitationRenderItem,
    precipitation_distance_opacity_factor,
    precipitation_streak_count,
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
        top_nx, top_ny = altaz_to_normalized_xy(
            column.top_alt_deg,
            column.top_az_deg,
            view_center,
            edge_fov_deg=float(viewer.edge_fov_deg),
        )
        base_x, base_y = normalized_to_screen_xy(base_nx, base_ny, geometry)
        top_x, top_y = normalized_to_screen_xy(top_nx, top_ny, geometry)
        distance_opacity = precipitation_distance_opacity_factor(column.distance_km)
        alpha = int(
            round(255.0 * min(1.0, max(0.0, opacity * distance_opacity)))
        )
        line_width = max(0.5, 2.2 * float(line_width_scale))
        solid_pen = QPen(
            QColor(*PRECIPITATION_COLUMN_DARK_COLOR_RGB, int(round(alpha * 0.7)))
        )
        solid_pen.setWidthF(line_width)
        solid_pen.setCosmetic(True)
        solid_pen.setCapStyle(Qt.PenCapStyle.FlatCap)
        dashed_pen = QPen(
            QColor(*PRECIPITATION_COLUMN_COLOR_RGB, int(round(alpha * 0.7)))
        )
        dashed_pen.setWidthF(line_width * 0.4)
        dashed_pen.setCosmetic(True)
        dashed_pen.setCapStyle(Qt.PenCapStyle.FlatCap)
        dashed_pen.setStyle(Qt.PenStyle.CustomDashLine)
        dashed_pen.setDashPattern([0.2, 2.8])
        streak_count = precipitation_streak_count(column.rate_mm_h)
        spacing_px = 4.0 * float(line_width_scale)
        center_offset = 0.5 * float(streak_count - 1)
        slant_px = max(2.0, abs(float(top_y) - float(base_y)) * 0.3)
        for index in range(streak_count):
            offset_x = (float(index) - center_offset) * spacing_px
            painter.setPen(solid_pen)
            painter.drawLine(
                QPointF(float(base_x) + offset_x, float(base_y)),
                QPointF(float(top_x) + offset_x + slant_px, float(top_y)),
            )
            painter.setPen(dashed_pen)
            painter.drawLine(
                QPointF(float(base_x) + offset_x, float(base_y)),
                QPointF(float(top_x) + offset_x + slant_px, float(top_y)),
            )
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
    line_width = max(0.5, 2.2 * float(line_width_scale) * marker_scale)
    solid_pen = QPen(
        QColor(*PRECIPITATION_COLUMN_DARK_COLOR_RGB, int(round(alpha * 0.7)))
    )
    solid_pen.setWidthF(line_width)
    solid_pen.setCosmetic(True)
    solid_pen.setCapStyle(Qt.PenCapStyle.FlatCap)
    dashed_pen = QPen(
        QColor(*PRECIPITATION_COLUMN_COLOR_RGB, int(round(alpha * 0.7)))
    )
    dashed_pen.setWidthF(line_width * 0.4)
    dashed_pen.setCosmetic(True)
    dashed_pen.setCapStyle(Qt.PenCapStyle.FlatCap)
    dashed_pen.setStyle(Qt.PenStyle.CustomDashLine)
    dashed_pen.setDashPattern([0.2, 2.8])

    streak_count = precipitation_streak_count(marker.rate_mm_h)
    spacing_px = 4.0 * float(line_width_scale) * marker_scale
    center_offset = 0.5 * float(streak_count - 1)
    height_px = (
        float(geometry.radius)
        * PRECIPITATION_NEAR_STREAK_HEIGHT_DEG
        / max(1.0e-6, float(edge_fov_deg))
        * marker_scale
    )
    slant_px = max(2.0, height_px * 0.3)
    center_x, center_y = geometry.center
    for index in range(streak_count):
        offset_x = (float(index) - center_offset) * spacing_px
        painter.setPen(solid_pen)
        painter.drawLine(
            QPointF(
                float(center_x) + offset_x - (slant_px * 0.5),
                float(center_y) + (height_px * 0.5),
            ),
            QPointF(
                float(center_x) + offset_x + (slant_px * 0.5),
                float(center_y) - (height_px * 0.5),
            ),
        )
        painter.setPen(dashed_pen)
        painter.drawLine(
            QPointF(
                float(center_x) + offset_x - (slant_px * 0.5),
                float(center_y) + (height_px * 0.5),
            ),
            QPointF(
                float(center_x) + offset_x + (slant_px * 0.5),
                float(center_y) - (height_px * 0.5),
            ),
        )
