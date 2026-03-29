from typing import Callable, List, Tuple

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import FIELD_OF_VIEW_DEG, TERRAIN_HORIZON_LINE_COLOR, URBAN_OUTLINE_LAYER_LINE_COLOR
from ..types import ScreenGeometry, UrbanOutlinePolyline
from .geometry import normalized_to_screen_xy
from .guides import split_by_gaps


def _urban_outline_height_alpha_scale(height_m: float) -> float:
    clamped_height_m = max(0.0, min(50.0, float(height_m)))
    return 0.25 + 0.75 * (clamped_height_m / 50.0)


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


def draw_terrain_horizon_line(
    painter: QPainter,
    geometry: ScreenGeometry,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    view_center: tuple[float, float],
    *,
    opacity: float = 1.0,
    line_width_scale: float = 1.0,
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
    is_in_fov_func: Callable[..., bool] = is_in_fov,
    altaz_to_normalized_xy_func: Callable[[float, float, Tuple[float, float]], Tuple[float, float]] = altaz_to_normalized_xy,
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], Tuple[float, float]] = normalized_to_screen_xy,
    split_by_gaps_func: Callable[[List[Tuple[float, float]]], List[List[Tuple[float, float]]]] = split_by_gaps,
) -> None:
    """Draw a terrain horizon polyline as an extra overlay over the geometric horizon."""
    if not terrain_profile_altaz or opacity <= 0.0:
        return
    effective_opacity = max(0.0, min(1.0, opacity * 0.7))
    if effective_opacity <= 0.0:
        return

    points: list[tuple[float, float]] = []
    for alt, az in terrain_profile_altaz:
        if not is_in_fov_func(float(alt), float(az), view_center, fov_deg=content_fov_deg):
            continue
        nx, ny = altaz_to_normalized_xy_func(float(alt), float(az), view_center)
        points.append((nx, ny))

    if len(points) < 2:
        return

    color = QColor(*TERRAIN_HORIZON_LINE_COLOR)
    color.setAlphaF(max(effective_opacity, min(1.0, 0.42 + (opacity * 0.95))))
    outline = QColor(*TERRAIN_HORIZON_LINE_COLOR)
    outline.setAlpha(max(0, min(255, int(round(135.0 * effective_opacity + 35.0)))))
    width_scale = max(1.0, float(line_width_scale))
    painter.save()
    for frag in split_by_gaps_func(points):
        if len(frag) < 2:
            continue
        pts = [QPointF(*normalized_to_screen_xy_func(nx, ny, geometry)) for nx, ny in frag]
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
    is_in_fov_func: Callable[..., bool] = is_in_fov,
    altaz_to_normalized_xy_func: Callable[[float, float, Tuple[float, float]], Tuple[float, float]] = altaz_to_normalized_xy,
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], Tuple[float, float]] = normalized_to_screen_xy,
    split_by_gaps_func: Callable[[List[Tuple[float, float]]], List[List[Tuple[float, float]]]] = split_by_gaps,
    minimal_azimuth_cover_func: Callable[[List[float]], Tuple[float, float, float]] = _minimal_azimuth_cover,
    urban_outline_height_alpha_scale_func: Callable[[float], float] = _urban_outline_height_alpha_scale,
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
        effective_opacity = float(opacity) * urban_outline_height_alpha_scale_func(height_m)
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
        az_start_deg, az_end_deg, az_span_deg = minimal_azimuth_cover_func(azimuth_deg)
        if az_span_deg < 0.5:
            representative_alt_deg = float(sum(float(alt) for alt, _az in outline) / len(outline))
            if representative_alt_deg >= -60.0:
                start_nx, start_ny = altaz_to_normalized_xy_func(representative_alt_deg, az_start_deg, view_center)
                end_nx, end_ny = altaz_to_normalized_xy_func(representative_alt_deg, az_end_deg, view_center)
                x1, y1 = normalized_to_screen_xy_func(start_nx, start_ny, geometry)
                x2, y2 = normalized_to_screen_xy_func(end_nx, end_ny, geometry)
                painter.setPen(simplified_pen)
                y = (float(y1) + float(y2)) * 0.5
                painter.drawPolyline(QPolygonF([QPointF(float(x1), y), QPointF(float(x2), y)]))
            continue

        painter.setPen(pen)
        points: list[tuple[float, float]] = []
        for alt, az in outline:
            if float(alt) < -60.0 or not is_in_fov_func(float(alt), float(az), view_center, fov_deg=content_fov_deg):
                if len(points) >= 2:
                    for frag in split_by_gaps_func(points):
                        if len(frag) < 2:
                            continue
                        screen_points = [QPointF(*normalized_to_screen_xy_func(nx, ny, geometry)) for nx, ny in frag]
                        painter.drawPolyline(QPolygonF(screen_points))
                points = []
                continue
            nx, ny = altaz_to_normalized_xy_func(float(alt), float(az), view_center)
            points.append((nx, ny))
        if len(points) >= 2:
            for frag in split_by_gaps_func(points):
                if len(frag) < 2:
                    continue
                screen_points = [QPointF(*normalized_to_screen_xy_func(nx, ny, geometry)) for nx, ny in frag]
                painter.drawPolyline(QPolygonF(screen_points))
    painter.restore()
