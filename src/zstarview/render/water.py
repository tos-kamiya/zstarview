from __future__ import annotations

import math
from typing import Callable, Sequence

from PySide6.QtCore import QPointF
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy
from ..paths import FIELD_OF_VIEW_DEG
from ..types import ScreenGeometry
from ..water_surface_mesh import WaterSurfaceMesh
from .geometry import normalized_to_screen_xy

WATER_SURFACE_FILL_RGB = (126, 194, 226)
WATER_SURFACE_STROKE_RGB = (74, 144, 194)


def draw_water_surface_meshes(
    painter: QPainter,
    geometry: ScreenGeometry,
    water_surfaces: Sequence[WaterSurfaceMesh] | None,
    view_center: tuple[float, float],
    *,
    observer_elevation_m: float = 0.0,
    opacity: float = 0.35,
    line_width_scale: float = 1.0,
    edge_fov_deg: float = FIELD_OF_VIEW_DEG,
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
    project_triangle_func: Callable[..., tuple[QPointF, QPointF, QPointF]] | None = None,
) -> None:
    if not water_surfaces:
        return
    layer_opacity = max(0.0, min(1.0, float(opacity)))
    if layer_opacity <= 0.0:
        return
    if project_triangle_func is None:
        project_triangle_func = _default_project_triangle

    painter.save()
    fill_color = QColor(*WATER_SURFACE_FILL_RGB)
    fill_color.setAlpha(int(round(255.0 * layer_opacity)))
    stroke_color = QColor(*WATER_SURFACE_STROKE_RGB)
    stroke_color.setAlpha(int(round(255.0 * min(1.0, layer_opacity * 0.9))))
    pen = QPen(stroke_color, 1.0)
    pen.setCosmetic(True)
    painter.setPen(pen)
    painter.setBrush(fill_color)
    for surface in water_surfaces:
        for triangle in surface.triangles_xy_m:
            screen_triangle = project_triangle_func(
                triangle,
                surface,
                view_center=view_center,
                observer_elevation_m=observer_elevation_m,
                edge_fov_deg=edge_fov_deg,
                content_fov_deg=content_fov_deg,
                geometry=geometry,
            )
            if len(screen_triangle) != 3:
                continue
            painter.drawPolygon(QPolygonF(screen_triangle))
    painter.restore()


def _default_project_triangle(
    triangle: tuple[tuple[float, float], tuple[float, float], tuple[float, float]],
    surface: WaterSurfaceMesh,
    *,
    view_center: tuple[float, float],
    observer_elevation_m: float,
    edge_fov_deg: float,
    content_fov_deg: float,
    geometry: ScreenGeometry,
) -> tuple[QPointF, QPointF, QPointF]:
    projected: list[QPointF] = []
    surface_elevation_m = float(getattr(surface, "surface_elevation_m", 0.0))
    for x_m, y_m in triangle:
        distance_m = math.hypot(float(x_m), float(y_m))
        azimuth_deg = math.degrees(math.atan2(float(x_m), float(y_m))) % 360.0
        altitude_deg = math.degrees(
            math.atan2(surface_elevation_m - float(observer_elevation_m), max(1.0e-6, distance_m))
        )
        nx, ny = altaz_to_normalized_xy(
            altitude_deg,
            azimuth_deg,
            view_center,
            edge_fov_deg=float(content_fov_deg if content_fov_deg > 0.0 else edge_fov_deg),
        )
        screen_x, screen_y = normalized_to_screen_xy(nx, ny, geometry)
        projected.append(QPointF(screen_x, screen_y))
    return tuple(projected)  # type: ignore[return-value]
