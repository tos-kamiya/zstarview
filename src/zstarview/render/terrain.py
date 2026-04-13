import math
from typing import Callable, List, Tuple

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import (
    FIELD_OF_VIEW_DEG,
    TERRAIN_HORIZON_LINE_COLOR,
    URBAN_OUTLINE_LAYER_LINE_COLOR,
)
from ..types import ScreenGeometry, UrbanOutlinePolyline
from .geometry import normalized_to_screen_xy
from .guides import split_by_gaps

URBAN_OUTLINE_FOREGROUND_MIN_WIDTH = 0.66
URBAN_OUTLINE_FOREGROUND_MAX_WIDTH = 1.14
URBAN_OUTLINE_FOREGROUND_CORE_WIDTH = 0.66
URBAN_OUTLINE_UNDERLAY_MIN_WIDTH = 1.2
URBAN_OUTLINE_UNDERLAY_WIDTH = 4.8
URBAN_OUTLINE_UNDERLAY_MID_WIDTH = 2.4
URBAN_OUTLINE_UNDERLAY_OUTER_WIDTH = 5.8
URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM = 0.01
URBAN_OUTLINE_NEAR_DISTANCE_KM = 0.5


def _urban_outline_foreground_alpha(opacity: float) -> float:
    return max(0.0, min(1.0, float(opacity)))


def _urban_outline_underlay_alpha(opacity: float) -> float:
    opacity = max(0.0, min(1.0, float(opacity)))
    return max(0.0, min(1.0, 0.01 + (opacity * 0.05)))


def _urban_outline_uses_underlay(distance_km: float) -> bool:
    return math.isfinite(float(distance_km)) and float(distance_km) <= URBAN_OUTLINE_NEAR_DISTANCE_KM


def _urban_outline_underlay_width(distance_km: float, *, width_scale: float = 1.0) -> float:
    if not math.isfinite(float(distance_km)):
        return 0.0
    d = max(
        URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM,
        min(URBAN_OUTLINE_NEAR_DISTANCE_KM, float(distance_km)),
    )
    span = URBAN_OUTLINE_NEAR_DISTANCE_KM - URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM
    if span <= 0.0:
        return URBAN_OUTLINE_UNDERLAY_WIDTH * max(1.0, float(width_scale))
    t = (URBAN_OUTLINE_NEAR_DISTANCE_KM - d) / span
    base_width = URBAN_OUTLINE_UNDERLAY_MIN_WIDTH + (
        (URBAN_OUTLINE_UNDERLAY_WIDTH - URBAN_OUTLINE_UNDERLAY_MIN_WIDTH) * t
    )
    return base_width * max(1.0, float(width_scale))


def _urban_outline_mid_width(distance_km: float, *, width_scale: float = 1.0) -> float:
    if not math.isfinite(float(distance_km)):
        return 0.0
    d = max(
        URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM,
        min(URBAN_OUTLINE_NEAR_DISTANCE_KM, float(distance_km)),
    )
    span = URBAN_OUTLINE_NEAR_DISTANCE_KM - URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM
    if span <= 0.0:
        return URBAN_OUTLINE_UNDERLAY_MID_WIDTH * max(1.0, float(width_scale))
    t = (URBAN_OUTLINE_NEAR_DISTANCE_KM - d) / span
    base_width = 1.0 + ((URBAN_OUTLINE_UNDERLAY_MID_WIDTH - 1.0) * t)
    return base_width * max(1.0, float(width_scale))


def _urban_outline_outer_width(distance_km: float, *, width_scale: float = 1.0) -> float:
    if not math.isfinite(float(distance_km)):
        return 0.0
    d = max(
        URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM,
        min(URBAN_OUTLINE_NEAR_DISTANCE_KM, float(distance_km)),
    )
    span = URBAN_OUTLINE_NEAR_DISTANCE_KM - URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM
    if span <= 0.0:
        return URBAN_OUTLINE_UNDERLAY_OUTER_WIDTH * max(1.0, float(width_scale))
    t = (URBAN_OUTLINE_NEAR_DISTANCE_KM - d) / span
    base_width = 1.2 + ((URBAN_OUTLINE_UNDERLAY_OUTER_WIDTH - 1.2) * t)
    return base_width * max(1.0, float(width_scale))


def _urban_outline_foreground_width(distance_km: float, *, width_scale: float = 1.0) -> float:
    if not math.isfinite(float(distance_km)):
        return URBAN_OUTLINE_FOREGROUND_MAX_WIDTH * max(1.0, float(width_scale))
    d = max(
        URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM,
        min(URBAN_OUTLINE_NEAR_DISTANCE_KM, float(distance_km)),
    )
    span = URBAN_OUTLINE_NEAR_DISTANCE_KM - URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM
    if span <= 0.0:
        base_width = URBAN_OUTLINE_FOREGROUND_MAX_WIDTH
    else:
        t = (URBAN_OUTLINE_NEAR_DISTANCE_KM - d) / span
        base_width = URBAN_OUTLINE_FOREGROUND_MIN_WIDTH + (
            (URBAN_OUTLINE_FOREGROUND_MAX_WIDTH - URBAN_OUTLINE_FOREGROUND_MIN_WIDTH) * t
        )
    return base_width * max(1.0, float(width_scale))


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


def terrain_horizon_line_alpha(opacity: float) -> float:
    """Return the alpha curve used by the terrain horizon foreground line."""
    opacity = max(0.0, min(1.0, float(opacity)))
    return max(opacity, min(1.0, 0.42 + (opacity * 0.95)))


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
    effective_opacity = max(0.0, min(1.0, float(opacity)))
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
    color.setAlphaF(terrain_horizon_line_alpha(effective_opacity))
    outline = QColor(*TERRAIN_HORIZON_LINE_COLOR)
    outline.setAlpha(max(0, min(255, int(round(135.0 * effective_opacity + 35.0)))))
    width_scale = max(1.0, float(line_width_scale))
    painter.save()
    for frag in split_by_gaps_func(points):
        if len(frag) < 2:
            continue
        pts = [QPointF(*normalized_to_screen_xy_func(nx, ny, geometry)) for nx, ny in frag]
        poly = QPolygonF(pts)

        base = QPen(outline, 3.6 * width_scale, Qt.PenStyle.SolidLine)
        base.setCosmetic(True)
        base.setCapStyle(Qt.PenCapStyle.RoundCap)
        base.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        painter.setPen(base)
        painter.drawPolyline(poly)

        fg = QPen(color, 1.2 * width_scale, Qt.PenStyle.SolidLine)
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
        distance_km = float(getattr(outline_entry, "distance_km", float("inf")))
        if len(outline) < 2:
            continue
        foreground_color = QColor(*URBAN_OUTLINE_LAYER_LINE_COLOR)
        foreground_color.setAlpha(
            max(0, min(255, int(round(255.0 * _urban_outline_foreground_alpha(opacity)))))
        )
        foreground_pen = QPen(
            foreground_color,
            _urban_outline_foreground_width(distance_km, width_scale=width_scale),
            Qt.PenStyle.SolidLine,
        )
        foreground_pen.setCosmetic(True)
        foreground_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        foreground_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)

        underlay_pen = None
        mid_underlay_pen = None
        outer_underlay_pen = None
        if _urban_outline_uses_underlay(distance_km):
            underlay_color = QColor(*URBAN_OUTLINE_LAYER_LINE_COLOR)
            underlay_color.setAlpha(
                max(0, min(255, int(round(255.0 * _urban_outline_underlay_alpha(opacity)))))
            )
            underlay_pen = QPen(
                underlay_color,
                _urban_outline_underlay_width(distance_km, width_scale=width_scale),
                Qt.PenStyle.SolidLine,
            )
            underlay_pen.setCosmetic(True)
            underlay_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            underlay_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            mid_underlay_color = QColor(*URBAN_OUTLINE_LAYER_LINE_COLOR)
            mid_underlay_color.setAlpha(
                max(0, min(255, int(round(255.0 * (0.25 * _urban_outline_underlay_alpha(opacity))))))
            )
            mid_underlay_pen = QPen(
                mid_underlay_color,
                _urban_outline_mid_width(distance_km, width_scale=width_scale),
                Qt.PenStyle.SolidLine,
            )
            mid_underlay_pen.setCosmetic(True)
            mid_underlay_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            mid_underlay_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            outer_underlay_color = QColor(*URBAN_OUTLINE_LAYER_LINE_COLOR)
            outer_underlay_color.setAlpha(
                max(0, min(255, int(round(255.0 * (0.12 * _urban_outline_underlay_alpha(opacity))))))
            )
            outer_underlay_pen = QPen(
                outer_underlay_color,
                _urban_outline_outer_width(distance_km, width_scale=width_scale),
                Qt.PenStyle.SolidLine,
            )
            outer_underlay_pen.setCosmetic(True)
            outer_underlay_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            outer_underlay_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)

        def _draw_fragments(fragments: list[list[tuple[float, float]]], pen: QPen) -> None:
            painter.setPen(pen)
            for frag in fragments:
                if len(frag) < 2:
                    continue
                screen_points = [
                    QPointF(*normalized_to_screen_xy_func(nx, ny, geometry))
                    for nx, ny in frag
                ]
                painter.drawPolyline(QPolygonF(screen_points))

        def _draw_points(points: list[tuple[float, float]]) -> None:
            if len(points) < 2:
                return
            fragments = [points] if len(points) == 2 else split_by_gaps_func(points)
            if outer_underlay_pen is not None:
                _draw_fragments(fragments, outer_underlay_pen)
            if mid_underlay_pen is not None:
                _draw_fragments(fragments, mid_underlay_pen)
            if underlay_pen is not None:
                _draw_fragments(fragments, underlay_pen)
            _draw_fragments(fragments, foreground_pen)

        azimuth_deg = [float(az) % 360.0 for _alt, az in outline]
        az_start_deg, az_end_deg, az_span_deg = minimal_azimuth_cover_func(azimuth_deg)
        if az_span_deg < 0.5:
            representative_alt_deg = float(sum(float(alt) for alt, _az in outline) / len(outline))
            if representative_alt_deg >= -60.0:
                start_nx, start_ny = altaz_to_normalized_xy_func(representative_alt_deg, az_start_deg, view_center)
                end_nx, end_ny = altaz_to_normalized_xy_func(representative_alt_deg, az_end_deg, view_center)
                x1, y1 = normalized_to_screen_xy_func(start_nx, start_ny, geometry)
                x2, y2 = normalized_to_screen_xy_func(end_nx, end_ny, geometry)
                y = (float(y1) + float(y2)) * 0.5
                _draw_points([(float(x1), y), (float(x2), y)])
            continue

        points: list[tuple[float, float]] = []
        for alt, az in outline:
            if float(alt) < -60.0 or not is_in_fov_func(float(alt), float(az), view_center, fov_deg=content_fov_deg):
                if len(points) >= 2:
                    _draw_points(points)
                points = []
                continue
            nx, ny = altaz_to_normalized_xy_func(float(alt), float(az), view_center)
            points.append((nx, ny))
        if len(points) >= 2:
            _draw_points(points)
    painter.restore()
