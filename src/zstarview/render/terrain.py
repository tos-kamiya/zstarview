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

URBAN_OUTLINE_FOREGROUND_MIN_WIDTH = 1.32
URBAN_OUTLINE_FOREGROUND_MAX_WIDTH = 2.28
URBAN_OUTLINE_FOREGROUND_CORE_WIDTH = 0.96
URBAN_OUTLINE_UNDERLAY_MIN_WIDTH = 2.0
URBAN_OUTLINE_UNDERLAY_WIDTH = 4.4
URBAN_OUTLINE_UNDERLAY_MID_WIDTH = 7.2
URBAN_OUTLINE_UNDERLAY_OUTER_WIDTH = 9.2
URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM = 0.01
URBAN_OUTLINE_NEAR_DISTANCE_KM = 0.5
URBAN_OUTLINE_HEIGHT_THICKEN_START_M = 100.0
URBAN_OUTLINE_HEIGHT_THICKEN_FULL_M = 600.0
TERRAIN_HORIZON_FAST_WIDTH = 3.6
TERRAIN_HORIZON_FAR_BASE_WIDTH = 1.9
TERRAIN_DISTANCE_BAND_OUTLINE_WIDTH = 2.0


def _urban_outline_foreground_alpha(opacity: float) -> float:
    return max(0.0, min(1.0, float(opacity)))


def _urban_outline_height_width_scale(height_m: float) -> float:
    height_m = float(height_m)
    if height_m <= URBAN_OUTLINE_HEIGHT_THICKEN_START_M:
        return 1.0
    if height_m >= URBAN_OUTLINE_HEIGHT_THICKEN_FULL_M:
        return 2.0
    span = URBAN_OUTLINE_HEIGHT_THICKEN_FULL_M - URBAN_OUTLINE_HEIGHT_THICKEN_START_M
    if span <= 0.0:
        return 2.0
    t = (height_m - URBAN_OUTLINE_HEIGHT_THICKEN_START_M) / span
    return 1.0 + t


def _urban_outline_underlay_alpha(opacity: float) -> float:
    opacity = max(0.0, min(1.0, float(opacity)))
    return max(0.0, min(1.0, 0.006 + (opacity * 0.03)))


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
        return URBAN_OUTLINE_UNDERLAY_WIDTH * float(width_scale)
    t = (URBAN_OUTLINE_NEAR_DISTANCE_KM - d) / span
    base_width = URBAN_OUTLINE_UNDERLAY_MIN_WIDTH + (
        (URBAN_OUTLINE_UNDERLAY_WIDTH - URBAN_OUTLINE_UNDERLAY_MIN_WIDTH) * t
    )
    return base_width * float(width_scale)


def _urban_outline_mid_width(distance_km: float, *, width_scale: float = 1.0) -> float:
    if not math.isfinite(float(distance_km)):
        return 0.0
    d = max(
        URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM,
        min(URBAN_OUTLINE_NEAR_DISTANCE_KM, float(distance_km)),
    )
    span = URBAN_OUTLINE_NEAR_DISTANCE_KM - URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM
    if span <= 0.0:
        return URBAN_OUTLINE_UNDERLAY_MID_WIDTH * float(width_scale)
    t = (URBAN_OUTLINE_NEAR_DISTANCE_KM - d) / span
    base_width = 2.0 + ((URBAN_OUTLINE_UNDERLAY_MID_WIDTH - 2.0) * t)
    return base_width * float(width_scale)


def _urban_outline_outer_width(distance_km: float, *, width_scale: float = 1.0) -> float:
    if not math.isfinite(float(distance_km)):
        return 0.0
    d = max(
        URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM,
        min(URBAN_OUTLINE_NEAR_DISTANCE_KM, float(distance_km)),
    )
    span = URBAN_OUTLINE_NEAR_DISTANCE_KM - URBAN_OUTLINE_UNDERLAY_MIN_DISTANCE_KM
    if span <= 0.0:
        return URBAN_OUTLINE_UNDERLAY_OUTER_WIDTH * float(width_scale)
    t = (URBAN_OUTLINE_NEAR_DISTANCE_KM - d) / span
    base_width = 2.4 + ((URBAN_OUTLINE_UNDERLAY_OUTER_WIDTH - 2.4) * t)
    return base_width * float(width_scale)


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
    return base_width * float(width_scale)


def _project_altaz_to_normalized_xy(
    alt_deg: float,
    az_deg: float,
    view_center: tuple[float, float],
    *,
    edge_fov_deg: float,
) -> tuple[float, float]:
    return altaz_to_normalized_xy(
        alt_deg,
        az_deg,
        view_center,
        edge_fov_deg=edge_fov_deg,
    )


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
    """Return the alpha curve used by the terrain horizon line."""
    opacity = max(0.0, min(1.0, float(opacity)))
    return max(opacity, min(1.0, 0.42 + (opacity * 0.95)))


def terrain_secondary_ridge_line_alpha(opacity: float) -> float:
    """Return the alpha curve used by the secondary ridge overlay."""
    opacity = max(0.0, min(1.0, float(opacity)))
    return max(0.0, min(1.0, 0.06 + (opacity * 0.22)))


def _distance_band_widths(
    band_index: int,
    band_count: int,
) -> float:
    return float(TERRAIN_DISTANCE_BAND_OUTLINE_WIDTH)


def _distance_band_alpha(
    band_index: int,
    band_count: int,
    opacity: float,
) -> float:
    near_alpha = terrain_horizon_line_alpha(opacity)
    far_alpha = terrain_secondary_ridge_line_alpha(opacity * 0.35)
    if band_count <= 1:
        return near_alpha
    t = max(0.0, min(1.0, float(band_index) / float(band_count - 1)))
    return near_alpha - (t * (near_alpha - far_alpha))


def _draw_terrain_profile_layer(
    painter: QPainter,
    geometry: ScreenGeometry,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    terrain_profile_distances_m: list[float] | None,
    view_center: tuple[float, float],
    *,
    opacity: float,
    base_width: float,
    far_base_width: float,
    fg_alpha: float,
    line_width_scale: float,
    color_rgb: tuple[int, int, int],
    fast_mode: bool,
    distance_widths: bool,
    edge_fov_deg: float,
    content_fov_deg: float,
    is_in_fov_func: Callable[..., bool],
    altaz_to_normalized_xy_func: Callable[[float, float, Tuple[float, float]], Tuple[float, float]],
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], Tuple[float, float]],
    split_by_gaps_func: Callable[[List[Tuple[float, float]]], List[List[Tuple[float, float]]]],
) -> None:
    if not terrain_profile_altaz or opacity <= 0.0:
        return
    effective_opacity = max(0.0, min(1.0, float(opacity)))
    if effective_opacity <= 0.0:
        return
    if terrain_profile_distances_m is not None and len(terrain_profile_distances_m) != len(terrain_profile_altaz):
        terrain_profile_distances_m = None

    points: list[tuple[float, float]] = []
    distances_m: list[float] = []
    for index, (alt, az) in enumerate(terrain_profile_altaz):
        if not is_in_fov_func(float(alt), float(az), view_center, fov_deg=content_fov_deg):
            continue
        try:
            nx, ny = altaz_to_normalized_xy_func(
                float(alt),
                float(az),
                view_center,
                edge_fov_deg=float(edge_fov_deg),
            )
        except TypeError:
            nx, ny = _project_altaz_to_normalized_xy(
                float(alt),
                float(az),
                view_center,
                edge_fov_deg=edge_fov_deg,
            )
        points.append((nx, ny))
        if terrain_profile_distances_m is not None:
            distances_m.append(float(terrain_profile_distances_m[index]))
        else:
            distances_m.append(float("nan"))

    if len(points) < 2:
        return

    color = QColor(*color_rgb)
    color.setAlphaF(max(0.0, min(1.0, float(fg_alpha))))
    width_scale = float(line_width_scale)
    has_distance_widths = distance_widths and (not fast_mode) and terrain_profile_distances_m is not None
    if has_distance_widths:
        valid_distances = [distance for distance in distances_m if math.isfinite(float(distance))]
        max_distance_m = max(valid_distances) if valid_distances else float("nan")
    else:
        max_distance_m = float("nan")
    painter.save()
    point_index = 0
    for frag in split_by_gaps_func(points):
        if len(frag) < 2:
            point_index += len(frag)
            continue
        frag_distances = distances_m[point_index:point_index + len(frag)]
        frag_points = [QPointF(*normalized_to_screen_xy_func(nx, ny, geometry)) for nx, ny in frag]
        if not has_distance_widths:
            poly = QPolygonF(frag_points)
            pen = QPen(color, float(base_width) * width_scale, Qt.PenStyle.SolidLine)
            pen.setCosmetic(True)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(pen)
            painter.drawPolyline(poly)
        else:
            if len(frag_points) >= 2:
                for (start, end, start_dist, end_dist) in zip(
                    frag_points,
                    frag_points[1:],
                    frag_distances,
                    frag_distances[1:],
                ):
                    if not (math.isfinite(float(start_dist)) and math.isfinite(float(end_dist))):
                        continue
                    segment_dist_m = 0.5 * (float(start_dist) + float(end_dist))
                    if not math.isfinite(max_distance_m) or max_distance_m <= 0.0:
                        t = 0.0
                    else:
                        t = max(0.0, min(1.0, segment_dist_m / max_distance_m))
                    base_width_m = float(base_width) - (t * (float(base_width) - float(far_base_width)))
                    pen = QPen(color, base_width_m * width_scale, Qt.PenStyle.SolidLine)
                    pen.setCosmetic(True)
                    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
                    painter.setPen(pen)
                    painter.drawLine(start, end)
        point_index += len(frag)
    painter.restore()


def draw_terrain_horizon_line(
    painter: QPainter,
    geometry: ScreenGeometry,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    terrain_profile_distances_m: list[float] | None,
    view_center: tuple[float, float],
    *,
    opacity: float = 1.0,
    line_width_scale: float = 1.0,
    fast_mode: bool = False,
    edge_fov_deg: float = FIELD_OF_VIEW_DEG,
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
    is_in_fov_func: Callable[..., bool] = is_in_fov,
    altaz_to_normalized_xy_func: Callable[[float, float, Tuple[float, float]], Tuple[float, float]] = altaz_to_normalized_xy,
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], Tuple[float, float]] = normalized_to_screen_xy,
    split_by_gaps_func: Callable[[List[Tuple[float, float]]], List[List[Tuple[float, float]]]] = split_by_gaps,
) -> None:
    """Draw a terrain horizon polyline as an extra overlay over the geometric horizon."""
    _draw_terrain_profile_layer(
        painter,
        geometry,
        terrain_profile_altaz,
        terrain_profile_distances_m,
        view_center,
        opacity=opacity,
        base_width=TERRAIN_HORIZON_FAST_WIDTH,
        far_base_width=TERRAIN_HORIZON_FAR_BASE_WIDTH,
        fg_alpha=terrain_horizon_line_alpha(opacity),
        line_width_scale=line_width_scale,
        color_rgb=TERRAIN_HORIZON_LINE_COLOR,
        fast_mode=fast_mode,
        distance_widths=True,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
        is_in_fov_func=is_in_fov_func,
        altaz_to_normalized_xy_func=altaz_to_normalized_xy_func,
        normalized_to_screen_xy_func=normalized_to_screen_xy_func,
        split_by_gaps_func=split_by_gaps_func,
    )


def draw_terrain_secondary_ridges(
    painter: QPainter,
    geometry: ScreenGeometry,
    terrain_secondary_profile_layers: list[list[tuple[float, float]]] | None,
    terrain_secondary_profile_distances_m_layers: list[list[float]] | None,
    view_center: tuple[float, float],
    *,
    opacity: float = 0.25,
    line_width_scale: float = 1.0,
    fast_mode: bool = False,
    edge_fov_deg: float = FIELD_OF_VIEW_DEG,
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
    is_in_fov_func: Callable[..., bool] = is_in_fov,
    altaz_to_normalized_xy_func: Callable[[float, float, Tuple[float, float]], Tuple[float, float]] = altaz_to_normalized_xy,
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], Tuple[float, float]] = normalized_to_screen_xy,
    split_by_gaps_func: Callable[[List[Tuple[float, float]]], List[List[Tuple[float, float]]]] = split_by_gaps,
) -> None:
    """Draw fixed-width ridge bands grouped by distance interval."""
    if fast_mode or not terrain_secondary_profile_layers or opacity <= 0.0:
        return

    ridge_opacity = terrain_secondary_ridge_line_alpha(opacity)
    if ridge_opacity <= 0.0:
        return

    ridge_color_rgb = TERRAIN_HORIZON_LINE_COLOR
    if terrain_secondary_profile_distances_m_layers is not None and len(terrain_secondary_profile_distances_m_layers) != len(terrain_secondary_profile_layers):
        terrain_secondary_profile_distances_m_layers = None
    layer_count = len(terrain_secondary_profile_layers)
    for layer_index, layer in enumerate(terrain_secondary_profile_layers):
        layer_distances_m = None
        if terrain_secondary_profile_distances_m_layers is not None:
            layer_distances_m = terrain_secondary_profile_distances_m_layers[layer_index]
        base_width = _distance_band_widths(
            layer_index,
            layer_count,
        )
        band_alpha = _distance_band_alpha(
            layer_index,
            layer_count,
            opacity,
        )
        _draw_terrain_profile_layer(
            painter,
            geometry,
            layer,
            layer_distances_m,
            view_center,
            opacity=opacity,
            base_width=base_width,
            far_base_width=base_width,
            fg_alpha=band_alpha,
            line_width_scale=line_width_scale,
            color_rgb=ridge_color_rgb,
            fast_mode=fast_mode,
            distance_widths=False,
            edge_fov_deg=edge_fov_deg,
            content_fov_deg=content_fov_deg,
            is_in_fov_func=is_in_fov_func,
            altaz_to_normalized_xy_func=altaz_to_normalized_xy_func,
            normalized_to_screen_xy_func=normalized_to_screen_xy_func,
            split_by_gaps_func=split_by_gaps_func,
        )


def draw_urban_outlines(
    painter: QPainter,
    geometry: ScreenGeometry,
    urban_outlines: list[UrbanOutlinePolyline] | None,
    view_center: tuple[float, float],
    *,
    opacity: float = 0.2,
    line_width_scale: float = 1.0,
    edge_fov_deg: float = FIELD_OF_VIEW_DEG,
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
    width_scale = float(line_width_scale)
    for outline_entry in urban_outlines:
        outline = list(outline_entry.points)
        distance_km = float(getattr(outline_entry, "distance_km", float("inf")))
        if len(outline) < 2:
            continue
        height_scale = _urban_outline_height_width_scale(float(getattr(outline_entry, "height_m", 0.0)))
        thickened_width_scale = width_scale * height_scale
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
                max(0, min(255, int(round(255.0 * (0.25 * _urban_outline_underlay_alpha(opacity))))))
            )
            underlay_pen = QPen(
                underlay_color,
                _urban_outline_underlay_width(distance_km, width_scale=thickened_width_scale),
                Qt.PenStyle.SolidLine,
            )
            underlay_pen.setCosmetic(True)
            underlay_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            underlay_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            mid_underlay_color = QColor(*URBAN_OUTLINE_LAYER_LINE_COLOR)
            mid_underlay_color.setAlpha(
                max(0, min(255, int(round(255.0 * (0.50 * _urban_outline_underlay_alpha(opacity))))))
            )
            mid_underlay_pen = QPen(
                mid_underlay_color,
                _urban_outline_mid_width(distance_km, width_scale=thickened_width_scale),
                Qt.PenStyle.SolidLine,
            )
            mid_underlay_pen.setCosmetic(True)
            mid_underlay_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            mid_underlay_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            outer_underlay_color = QColor(*URBAN_OUTLINE_LAYER_LINE_COLOR)
            outer_underlay_color.setAlpha(
                max(0, min(255, int(round(255.0 * (0.90 * _urban_outline_underlay_alpha(opacity))))))
            )
            outer_underlay_pen = QPen(
                outer_underlay_color,
                _urban_outline_outer_width(distance_km, width_scale=thickened_width_scale),
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
                try:
                    start_nx, start_ny = altaz_to_normalized_xy_func(
                        representative_alt_deg,
                        az_start_deg,
                        view_center,
                        edge_fov_deg=float(edge_fov_deg),
                    )
                    end_nx, end_ny = altaz_to_normalized_xy_func(
                        representative_alt_deg,
                        az_end_deg,
                        view_center,
                        edge_fov_deg=float(edge_fov_deg),
                    )
                except TypeError:
                    start_nx, start_ny = _project_altaz_to_normalized_xy(
                        representative_alt_deg,
                        az_start_deg,
                        view_center,
                        edge_fov_deg=edge_fov_deg,
                    )
                    end_nx, end_ny = _project_altaz_to_normalized_xy(
                        representative_alt_deg,
                        az_end_deg,
                        view_center,
                        edge_fov_deg=edge_fov_deg,
                    )
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
            try:
                nx, ny = altaz_to_normalized_xy_func(
                    float(alt),
                    float(az),
                    view_center,
                    edge_fov_deg=float(edge_fov_deg),
                )
            except TypeError:
                nx, ny = _project_altaz_to_normalized_xy(
                    float(alt),
                    float(az),
                    view_center,
                    edge_fov_deg=edge_fov_deg,
                )
            points.append((nx, ny))
        if len(points) >= 2:
            _draw_points(points)
    painter.restore()
