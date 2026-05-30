from __future__ import annotations

from dataclasses import dataclass

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..location_resolver.place_projection import project_place_targets_to_altaz
from ..paths import ThemeStyle
from ..types import ScreenGeometry, ViewerData
from ..tropical_cyclones.models import TropicalCycloneSnapshot
from .geometry import normalized_to_screen_xy

TROPICAL_CYCLONE_TARGET_HEIGHT_M = 0.0


def _debug_projection_line(
    *,
    kind: str,
    label: str,
    lat_deg: float,
    lon_deg: float,
    alt_deg: float,
    az_deg: float,
) -> None:
    print(
        (
            f"tropical_cyclone {kind} {label}: "
            f"lat={lat_deg:.3f} lon={lon_deg:.3f} alt={alt_deg:.3f} az={az_deg:.3f}"
        ),
        flush=True,
    )


@dataclass(frozen=True, slots=True)
class _RenderPoint:
    nx: float
    ny: float
    alt_deg: float
    az_deg: float
    distance_km: float


def _project_point(
    lat_deg: float,
    lon_deg: float,
    *,
    viewer: ViewerData,
    height_m: float,
) -> _RenderPoint | None:
    projections = project_place_targets_to_altaz(
        observer_latitude_deg=float(viewer.lat_deg),
        observer_longitude_deg=float(viewer.lon_deg),
        observer_height_m=float(viewer.ground_elevation_m),
        target_latitude_deg=[float(lat_deg)],
        target_longitude_deg=[float(lon_deg)],
        target_height_m=[float(height_m)],
    )
    if not projections:
        return None
    projection = projections[0]
    view_center = tuple(float(value) for value in viewer.view_center)
    if not is_in_fov(
        float(projection.alt_deg),
        float(projection.az_deg),
        view_center,
        fov_deg=float(viewer.content_fov_deg),
    ):
        return None
    nx, ny = altaz_to_normalized_xy(
        float(projection.alt_deg),
        float(projection.az_deg),
        view_center,
        edge_fov_deg=float(viewer.edge_fov_deg),
    )
    _debug_projection_line(
        kind="point",
        label="storm",
        lat_deg=float(lat_deg),
        lon_deg=float(lon_deg),
        alt_deg=float(projection.alt_deg),
        az_deg=float(projection.az_deg),
    )
    return _RenderPoint(
        nx=float(nx),
        ny=float(ny),
        alt_deg=float(projection.alt_deg),
        az_deg=float(projection.az_deg),
        distance_km=float(projection.distance_km),
    )


def _draw_line(
    painter: QPainter,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    geometry: ScreenGeometry,
    color_rgba: tuple[int, int, int, int],
    width_px: float,
) -> None:
    pen = QPen(QColor(*color_rgba), float(width_px), Qt.PenStyle.SolidLine)
    pen.setCosmetic(True)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    painter.setPen(pen)
    painter.drawLine(
        QPointF(*normalized_to_screen_xy(start[0], start[1], geometry)),
        QPointF(*normalized_to_screen_xy(end[0], end[1], geometry)),
    )


def _draw_inverted_triangle(
    painter: QPainter,
    point: _RenderPoint,
    *,
    geometry: ScreenGeometry,
    color_rgba: tuple[int, int, int, int],
) -> None:
    # Scale the symbol down as the storm is farther away from the observer.
    scale = max(0.45, min(1.15, 20.0 / max(1.0, point.distance_km)))
    body_height = 20.0 * scale
    body_width = 14.0 * scale
    tip = QPointF(*normalized_to_screen_xy(point.nx, point.ny, geometry))
    top_left = QPointF(tip.x() - body_width, tip.y() - body_height)
    top_right = QPointF(tip.x() + body_width, tip.y() - body_height)
    polygon = QPolygonF([tip, top_left, top_right])
    pen = QPen(QColor(*color_rgba), 1.6)
    pen.setCosmetic(True)
    painter.setPen(pen)
    painter.setBrush(QColor(*color_rgba))
    painter.drawPolygon(polygon)
    painter.drawLine(top_left, top_right)


def _draw_label(
    painter: QPainter,
    point: _RenderPoint,
    *,
    geometry: ScreenGeometry,
    color_rgba: tuple[int, int, int, int],
    text: str,
) -> None:
    x, y = normalized_to_screen_xy(point.nx, point.ny, geometry)
    painter.setPen(QColor(*color_rgba))
    painter.drawText(
        QPointF(x - 28.0, y + 20.0),
        text,
    )


def draw_tropical_cyclone_overlay(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    snapshot: TropicalCycloneSnapshot | None,
    theme: ThemeStyle,
    enabled: bool = True,
) -> None:
    if not enabled or snapshot is None:
        return
    del theme

    observed = snapshot.observed_position
    point = _project_point(
        observed.lat_deg,
        observed.lon_deg,
        viewer=viewer,
        height_m=TROPICAL_CYCLONE_TARGET_HEIGHT_M,
    )
    if point is None:
        return

    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    triangle_color = (255, 74, 74, 230)
    _draw_inverted_triangle(
        painter,
        point,
        geometry=geometry,
        color_rgba=triangle_color,
    )
    _draw_label(
        painter,
        point,
        geometry=geometry,
        color_rgba=triangle_color,
        text=snapshot.storm_name,
    )
    painter.restore()
