from __future__ import annotations

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..precipitation import (
    PRECIPITATION_COLUMN_COLOR_RGB,
    ProjectedPrecipitationColumn,
    precipitation_distance_opacity_factor,
    precipitation_streak_count,
)
from ..types import ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy


def draw_precipitation_columns(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    columns: list[ProjectedPrecipitationColumn] | None,
    *,
    opacity: float,
    line_width_scale: float = 1.0,
) -> None:
    if not columns or opacity <= 0.0:
        return
    view_center = tuple(viewer.view_center)
    painter.save()
    for column in columns:
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
        pen = QPen(
            QColor(
                *PRECIPITATION_COLUMN_COLOR_RGB,
                int(
                    round(
                        255.0
                        * min(1.0, max(0.0, opacity * distance_opacity))
                    )
                ),
            )
        )
        pen.setWidthF(max(0.5, 1.8 * float(line_width_scale)))
        pen.setCosmetic(True)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        streak_count = precipitation_streak_count(column.rate_mm_h)
        spacing_px = 4.0 * float(line_width_scale)
        center_offset = 0.5 * float(streak_count - 1)
        slant_px = max(2.0, abs(float(top_y) - float(base_y)) * 0.3)
        for index in range(streak_count):
            offset_x = (float(index) - center_offset) * spacing_px
            painter.drawLine(
                QPointF(float(base_x) + offset_x, float(base_y)),
                QPointF(float(top_x) + offset_x + slant_px, float(top_y)),
            )
    painter.restore()
