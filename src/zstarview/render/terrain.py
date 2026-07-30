import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import (
    TERRAIN_HORIZON_LINE_COLOR,
    URBAN_OUTLINE_LAYER_LINE_COLOR,
    OverlayLayerStyle,
)
from ..types import ScreenGeometry, UrbanOutlinePolyline, ViewerData
from ..water_overlay import WaterOverlayPoint, WaterOverlayPolyline
from .geometry import normalized_to_screen_xy
from .ground_mask import build_ground_mask
from .guides import _clip_polyline_to_radius, split_by_gaps
from .qt_image import np_rgba_to_qimage

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
URBAN_OUTLINE_FILL_ALPHA_FLOOR = 0.04
URBAN_OUTLINE_FILL_ALPHA_SCALE = 0.18
URBAN_OUTLINE_FILL_MAX_ENDPOINT_GAP_PX = 8.0
URBAN_OUTLINE_FILL_MIN_SCREEN_SPAN_PX = 8.0
TERRAIN_HORIZON_FAST_WIDTH = 3.6 / 3.0
TERRAIN_HORIZON_FAR_BASE_WIDTH = 1.9 / 3.0
TERRAIN_DISTANCE_BAND_NEAR_OUTLINE_WIDTH = 2.0
TERRAIN_DISTANCE_BAND_FAR_OUTLINE_WIDTH = 1.1
TERRAIN_DISTANCE_BAND_NEAR_DISTANCE_KM = 0.5
TERRAIN_DISTANCE_BAND_FAR_DISTANCE_KM = 128.0
TERRAIN_DISTANCE_BAND_REFERENCE_ALPHA_DISTANCE_KM = 4.0
TERRAIN_DISTANCE_BAND_WIDTH_DECAY_EXPONENT = 1.35
TERRAIN_DISTANCE_BAND_UNDERLAY_NEAR_SCALE = 1.15
TERRAIN_DISTANCE_BAND_UNDERLAY_FAR_SCALE = 1.55
TERRAIN_DISTANCE_BAND_UNDERLAY_NEAR_ALPHA_SCALE = 0.10
TERRAIN_DISTANCE_BAND_UNDERLAY_FAR_ALPHA_SCALE = 0.06
TERRAIN_DISTANCE_BAND_ALPHA_DECAY_EXPONENT = 1.85
WATER_OVERLAY_POINT_COLOR_RGB = (122, 218, 240)
# Keep sea tones aligned with dev-samples/basic-color-palette.html:
# HSV(191°, 49.2%, 94.1%) -> RGB(122, 218, 240).
WATER_OVERLAY_SEA_COLOR_RGB = (122, 218, 240)
WATER_OVERLAY_SEA_125_COLOR_RGB = (122, 218, 240)
WATER_OVERLAY_SEA_250_COLOR_RGB = (54, 200, 184)
WATER_OVERLAY_SEA_500_COLOR_RGB = (255, 170, 64)
# Inland water shares the river color in the current palette.
WATER_OVERLAY_LAKE_COLOR_RGB = (94, 214, 255)
WATER_OVERLAY_RIVER_COLOR_RGB = (94, 214, 255)
WATER_OVERLAY_POINT_RADIUS_PX = 3.0
WATER_OVERLAY_MARKER_MAJOR_RADIUS_SCALE = 0.59
WATER_OVERLAY_MARKER_MINOR_RADIUS_SCALE = 0.46
WATER_OVERLAY_MARKER_PEN_WIDTH_SCALE = 0.42
WATER_OVERLAY_DISTANCE_ALPHA_REFERENCE_KM = 128.0
WATER_OVERLAY_DISTANCE_ALPHA_REFERENCE_SCALE = 16.0


@dataclass(frozen=True)
class TerrainHorizonRenderSpec:
    opacity: float
    base_width: float
    far_base_width: float
    fg_alpha: float
    line_width_scale: float
    color_rgb: tuple[int, int, int]
    fast_mode: bool
    distance_widths: bool


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


def _urban_outline_fill_alpha(opacity: float) -> int:
    fill_opacity = max(
        0.0,
        min(1.0, URBAN_OUTLINE_FILL_ALPHA_FLOOR + (URBAN_OUTLINE_FILL_ALPHA_SCALE * float(opacity))),
    )
    return int(round(255.0 * fill_opacity))


def _urban_outline_fragment_is_closed_for_fill(screen_points: list[QPointF]) -> bool:
    if len(screen_points) < 4:
        return False
    first_point = screen_points[0]
    last_point = screen_points[-1]
    dx = float(first_point.x()) - float(last_point.x())
    dy = float(first_point.y()) - float(last_point.y())
    return ((dx * dx) + (dy * dy)) <= (URBAN_OUTLINE_FILL_MAX_ENDPOINT_GAP_PX * URBAN_OUTLINE_FILL_MAX_ENDPOINT_GAP_PX)


def _urban_outline_fragment_is_large_enough_for_fill(screen_points: list[QPointF]) -> bool:
    xs = [float(point.x()) for point in screen_points]
    ys = [float(point.y()) for point in screen_points]
    return (
        (max(xs) - min(xs)) >= URBAN_OUTLINE_FILL_MIN_SCREEN_SPAN_PX
        and (max(ys) - min(ys)) >= URBAN_OUTLINE_FILL_MIN_SCREEN_SPAN_PX
    )


def _urban_outline_underlay_alpha(opacity: float) -> float:
    opacity = max(0.0, min(1.0, float(opacity)))
    return max(0.0, min(1.0, 0.006 + (opacity * 0.03)))


def _dampen_alpha_for_narrow_width(alpha: float, width: float) -> float:
    width_alpha_scale = min(1.0, max(0.0, float(width)))
    return max(0.0, min(1.0, float(alpha) * width_alpha_scale))


def _viewer_projection_params(viewer: ViewerData) -> tuple[tuple[float, float], float, float]:
    view_center = tuple(float(value) for value in viewer.view_center)
    edge_fov_deg = float(viewer.edge_fov_deg)
    content_fov_deg = float(viewer.content_fov_deg)
    return view_center, edge_fov_deg, content_fov_deg


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


def _rotate_profile_to_seam_azimuth(
    samples: list[tuple[float, float, float]],
    *,
    seam_az_deg: float,
) -> list[tuple[float, float, float]]:
    if len(samples) < 3:
        return samples
    seam = float(seam_az_deg) % 360.0
    return sorted(
        samples,
        key=lambda item: (float(item[1]) - seam) % 360.0,
    )


def terrain_horizon_line_alpha(opacity: float) -> float:
    """Return the alpha curve used by the terrain horizon line."""
    opacity = max(0.0, min(1.0, float(opacity)))
    return max(opacity, min(1.0, 0.42 + (opacity * 0.95)))


def terrain_secondary_ridge_line_alpha(opacity: float) -> float:
    """Return the alpha curve used by the distance-band secondary-ridge overlay."""
    opacity = max(0.0, min(1.0, float(opacity)))
    return max(0.0, min(1.0, 0.03 + (opacity * 0.12)))


def _azimuth_bin_key(az_deg: float, *, bin_size_deg: float = 1.0) -> int:
    if bin_size_deg <= 0.0:
        return 0
    az = float(az_deg) % 360.0
    return int(math.floor(az / float(bin_size_deg)))


def _circular_midpoint_azimuth_deg(start_az_deg: float, end_az_deg: float) -> float:
    start = math.radians(float(start_az_deg) % 360.0)
    end = math.radians(float(end_az_deg) % 360.0)
    sin_sum = math.sin(start) + math.sin(end)
    cos_sum = math.cos(start) + math.cos(end)
    if abs(sin_sum) < 1.0e-12 and abs(cos_sum) < 1.0e-12:
        return float(start_az_deg) % 360.0
    return math.degrees(math.atan2(sin_sum, cos_sum)) % 360.0


def _solid_pen(color_rgb: tuple[int, int, int], alpha: float, width: float) -> QPen:
    color = QColor(*color_rgb)
    color.setAlphaF(max(0.0, min(1.0, float(alpha))))
    pen = QPen(color, float(width), Qt.PenStyle.SolidLine)
    pen.setCosmetic(True)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    return pen


def _water_overlay_point_color_rgb(water_point: WaterOverlayPoint) -> tuple[int, int, int]:
    category = str(getattr(water_point, "water_category", "")).strip().lower()
    if category in {"sea", "sea-125", "sea-250", "sea-500"}:
        return WATER_OVERLAY_SEA_COLOR_RGB
    if category == "river":
        return WATER_OVERLAY_RIVER_COLOR_RGB
    if category == "lake":
        return WATER_OVERLAY_LAKE_COLOR_RGB
    return WATER_OVERLAY_POINT_COLOR_RGB


def _water_overlay_marker_geometry(
    line_width_scale: float,
    *,
    distance_m: float = 0.0,
) -> tuple[float, float, float]:
    scale = max(1.0, float(line_width_scale))
    base_radius = WATER_OVERLAY_POINT_RADIUS_PX * scale
    distance_scale = max(0.35, _water_overlay_distance_alpha_scale(distance_m))
    major_radius = base_radius * WATER_OVERLAY_MARKER_MAJOR_RADIUS_SCALE * distance_scale
    minor_radius = max(0.6, major_radius * max(0.2, 0.48 * distance_scale))
    pen_width = max(1.0, base_radius * WATER_OVERLAY_MARKER_PEN_WIDTH_SCALE)
    return major_radius, minor_radius, pen_width


def _water_overlay_marker_rotation_deg(
    px: float,
    py: float,
    geometry: ScreenGeometry,
) -> float:
    # Keep the marker horizontal; the projection foreshortens the vertical axis by distance.
    return 0.0


def _water_overlay_distance_alpha_scale(distance_m: float) -> float:
    distance_m = max(0.0, float(distance_m))
    if distance_m <= 0.0:
        return 1.0
    reference = max(1.0, float(WATER_OVERLAY_DISTANCE_ALPHA_REFERENCE_KM) * 1000.0)
    factor = max(1.0, float(WATER_OVERLAY_DISTANCE_ALPHA_REFERENCE_SCALE))
    return float(math.exp(-math.log(factor) * (distance_m / reference)))


def _distance_band_widths(
    *,
    distance_km: float,
    band_count: int,
) -> float:
    if band_count <= 1:
        return float(TERRAIN_DISTANCE_BAND_NEAR_OUTLINE_WIDTH)
    near_km = float(TERRAIN_DISTANCE_BAND_NEAR_DISTANCE_KM)
    far_km = float(TERRAIN_DISTANCE_BAND_FAR_DISTANCE_KM)
    d = max(near_km, min(far_km, float(distance_km)))
    span = math.log(far_km / near_km) if far_km > near_km else 0.0
    if span > 0.0:
        t = max(0.0, min(1.0, math.log(d / near_km) / span))
    else:
        t = 0.0
    eased_t = t ** TERRAIN_DISTANCE_BAND_WIDTH_DECAY_EXPONENT
    outline_width = TERRAIN_DISTANCE_BAND_NEAR_OUTLINE_WIDTH - (
        eased_t * (TERRAIN_DISTANCE_BAND_NEAR_OUTLINE_WIDTH - TERRAIN_DISTANCE_BAND_FAR_OUTLINE_WIDTH)
    )
    return float(outline_width)


def _distance_band_underlay_width(
    *,
    distance_km: float,
    band_count: int,
) -> float:
    if band_count <= 1:
        return float(TERRAIN_DISTANCE_BAND_NEAR_OUTLINE_WIDTH) * TERRAIN_DISTANCE_BAND_UNDERLAY_NEAR_SCALE
    near_km = float(TERRAIN_DISTANCE_BAND_NEAR_DISTANCE_KM)
    far_km = float(TERRAIN_DISTANCE_BAND_FAR_DISTANCE_KM)
    d = max(near_km, min(far_km, float(distance_km)))
    span = math.log(far_km / near_km) if far_km > near_km else 0.0
    if span > 0.0:
        t = max(0.0, min(1.0, math.log(d / near_km) / span))
    else:
        t = 0.0
    scale = TERRAIN_DISTANCE_BAND_UNDERLAY_NEAR_SCALE + (
        t * (TERRAIN_DISTANCE_BAND_UNDERLAY_FAR_SCALE - TERRAIN_DISTANCE_BAND_UNDERLAY_NEAR_SCALE)
    )
    return float(_distance_band_widths(distance_km=distance_km, band_count=band_count)) * scale


def _distance_band_underlay_alpha(
    *,
    distance_km: float,
    band_count: int,
    opacity: float,
) -> float:
    band_alpha = _distance_band_alpha(distance_km=distance_km, band_count=band_count, opacity=opacity)
    if band_count <= 1:
        return band_alpha * TERRAIN_DISTANCE_BAND_UNDERLAY_NEAR_ALPHA_SCALE
    near_km = float(TERRAIN_DISTANCE_BAND_NEAR_DISTANCE_KM)
    far_km = float(TERRAIN_DISTANCE_BAND_FAR_DISTANCE_KM)
    d = max(near_km, min(far_km, float(distance_km)))
    span = math.log(far_km / near_km) if far_km > near_km else 0.0
    if span > 0.0:
        t = max(0.0, min(1.0, math.log(d / near_km) / span))
    else:
        t = 0.0
    eased_t = t ** TERRAIN_DISTANCE_BAND_ALPHA_DECAY_EXPONENT
    scale = TERRAIN_DISTANCE_BAND_UNDERLAY_NEAR_ALPHA_SCALE + (
        eased_t * (TERRAIN_DISTANCE_BAND_UNDERLAY_FAR_ALPHA_SCALE - TERRAIN_DISTANCE_BAND_UNDERLAY_NEAR_ALPHA_SCALE)
    )
    return band_alpha * scale


def _distance_band_alpha(
    *,
    distance_km: float,
    band_count: int,
    opacity: float,
) -> float:
    near_alpha = terrain_horizon_line_alpha(opacity)
    far_alpha = terrain_secondary_ridge_line_alpha(opacity * 0.35)
    if band_count <= 1:
        return near_alpha
    near_km = float(TERRAIN_DISTANCE_BAND_NEAR_DISTANCE_KM)
    far_km = float(TERRAIN_DISTANCE_BAND_FAR_DISTANCE_KM)
    d = max(near_km, min(far_km, float(distance_km)))
    span = math.log(far_km / near_km) if far_km > near_km else 0.0
    if span > 0.0:
        t = max(0.0, min(1.0, math.log(d / near_km) / span))
    else:
        t = 0.0
    eased_t = t ** TERRAIN_DISTANCE_BAND_ALPHA_DECAY_EXPONENT
    return near_alpha - (eased_t * (near_alpha - far_alpha))


def _draw_terrain_profile_layer(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    terrain_profile_distances_m: list[float] | None,
    *,
    spec: TerrainHorizonRenderSpec,
    is_in_fov_func: Callable[..., bool],
    altaz_to_normalized_xy_func: Callable[[float, float, tuple[float, float]], tuple[float, float]],
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], tuple[float, float]],
    split_by_gaps_func: Callable[[list[tuple[float, float]]], list[list[tuple[float, float]]]],
) -> None:
    if not terrain_profile_altaz or spec.opacity <= 0.0:
        return
    effective_opacity = max(0.0, min(1.0, float(spec.opacity)))
    if effective_opacity <= 0.0:
        return
    if terrain_profile_distances_m is not None and len(terrain_profile_distances_m) != len(terrain_profile_altaz):
        terrain_profile_distances_m = None

    view_center, edge_fov_deg, content_fov_deg = _viewer_projection_params(viewer)

    samples: list[tuple[float, float, float]] = []
    for index, (alt, az) in enumerate(terrain_profile_altaz):
        if not is_in_fov_func(float(alt), float(az), view_center, fov_deg=content_fov_deg):
            continue
        if terrain_profile_distances_m is not None:
            distance_m = float(terrain_profile_distances_m[index])
        else:
            distance_m = float("nan")
        samples.append((float(alt), float(az), distance_m))

    if len(samples) < 2:
        return

    samples = _rotate_profile_to_seam_azimuth(
        samples,
        seam_az_deg=(float(view_center[1]) + 180.0) % 360.0,
    )
    projected_points: list[tuple[float, float]] = []
    distances_m: list[float] = []
    for alt, az, distance_m in samples:
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
        projected_points.append((nx, ny))
        distances_m.append(distance_m)
    points = projected_points

    color = QColor(*spec.color_rgb)
    color.setAlphaF(max(0.0, min(1.0, float(spec.fg_alpha))))
    width_scale = float(spec.line_width_scale)
    has_distance_widths = spec.distance_widths and (not spec.fast_mode) and terrain_profile_distances_m is not None
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
            pen = QPen(color, float(spec.base_width) * width_scale, Qt.PenStyle.SolidLine)
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
                    base_width_m = float(spec.base_width) - (t * (float(spec.base_width) - float(spec.far_base_width)))
                    pen = QPen(color, base_width_m * width_scale, Qt.PenStyle.SolidLine)
                    pen.setCosmetic(True)
                    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
                    painter.setPen(pen)
                    painter.drawLine(start, end)
        point_index += len(frag)
    painter.restore()


def draw_ground_tint(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    *,
    opacity: float = 0.15,
    tint_rgb: tuple[int, int, int] = (128, 128, 128),
) -> None:
    """Fill the region below the terrain horizon with a pixel-based tint."""
    if not terrain_profile_altaz or opacity <= 0.0:
        return
    viewport = painter.viewport()
    width = int(viewport.width())
    height = int(viewport.height())
    if width <= 0 or height <= 0:
        return
    view_center, edge_fov_deg, content_fov_deg = _viewer_projection_params(viewer)
    mask = build_ground_mask(
        width,
        height,
        geometry,
        view_center,
        terrain_profile_altaz,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
        origin=(float(viewport.x()), float(viewport.y())),
    )
    if not np.any(mask):
        return
    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    rgba[mask, :3] = np.asarray(tint_rgb, dtype=np.uint8)
    rgba[mask, 3] = int(round(255.0 * max(0.0, min(1.0, float(opacity)))))
    painter.drawImage(viewport.topLeft(), np_rgba_to_qimage(rgba))


def draw_terrain_secondary_ridges(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    terrain_secondary_ridges_layers: list[list[tuple[float, float]]] | None,
    terrain_secondary_ridges_distances_m_layers: list[list[float]] | None,
    *,
    opacity: float = 0.25,
    line_width_scale: float = 1.0,
    layer_style: OverlayLayerStyle | None = None,
    is_in_fov_func: Callable[..., bool] = is_in_fov,
    altaz_to_normalized_xy_func: Callable[[float, float, tuple[float, float]], tuple[float, float]] = altaz_to_normalized_xy,
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], tuple[float, float]] = normalized_to_screen_xy,
    split_by_gaps_func: Callable[[list[tuple[float, float]]], list[list[tuple[float, float]]]] = split_by_gaps,
) -> None:
    """Draw fixed-width ridge bands grouped by distance interval."""
    alpha_scale = 1.0 if layer_style is None else float(layer_style.alpha_scale)
    layer_opacity = float(opacity) * alpha_scale
    if layer_opacity <= 0.0:
        return

    if not terrain_secondary_ridges_layers:
        return

    if terrain_secondary_ridges_distances_m_layers is not None and len(terrain_secondary_ridges_distances_m_layers) != len(terrain_secondary_ridges_layers):
        terrain_secondary_ridges_distances_m_layers = None
    layer_count = len(terrain_secondary_ridges_layers)
    max_visible_alt_by_bin: dict[int, float] = {}
    view_center = tuple(float(value) for value in viewer.view_center)
    edge_fov_deg = float(viewer.edge_fov_deg)
    content_fov_deg = float(viewer.content_fov_deg)
    seam_az_deg = (float(view_center[1]) + 180.0) % 360.0
    alpha_distance_km = TERRAIN_DISTANCE_BAND_REFERENCE_ALPHA_DISTANCE_KM
    occlusion_bin_deg = 1.0
    hidden_altitude_delta_deg = 0.1
    hidden_alpha_scale = 0.48
    hidden_width_scale = 0.82
    layer_rgb = TERRAIN_HORIZON_LINE_COLOR if layer_style is None else layer_style.rgb
    style_width_scale = 1.0 if layer_style is None else float(layer_style.width_scale)

    def _draw_secondary_ridge_run(
        *,
        points: list[QPointF],
        is_hidden: bool,
        base_width: float,
        underlay_width: float,
        band_alpha: float,
        underlay_alpha: float,
    ) -> None:
        if len(points) < 2:
            return
        poly = QPolygonF(points)
        alpha_scale = hidden_alpha_scale if is_hidden else 1.0
        width_scale = hidden_width_scale if is_hidden else 1.0
        if underlay_alpha > 0.0 and underlay_width > base_width:
            painter.setPen(
                _solid_pen(
                    layer_rgb,
                    underlay_alpha * alpha_scale,
                    float(underlay_width)
                    * float(line_width_scale)
                    * style_width_scale
                    * width_scale,
                )
            )
            painter.drawPolyline(poly)
        painter.setPen(
            _solid_pen(
                layer_rgb,
                band_alpha * alpha_scale,
                float(base_width)
                * float(line_width_scale)
                * style_width_scale
                * width_scale,
            )
        )
        painter.drawPolyline(poly)

    for layer_index, layer in enumerate(terrain_secondary_ridges_layers):
        if terrain_secondary_ridges_distances_m_layers is not None:
            layer_distances = terrain_secondary_ridges_distances_m_layers[layer_index]
        else:
            layer_distances = [float("nan")] * len(layer)
        finite_layer_distances_km = [
            float(distance_m) / 1000.0
            for distance_m in layer_distances
            if math.isfinite(float(distance_m))
        ]
        representative_distance_km = max(finite_layer_distances_km) if finite_layer_distances_km else None
        layer_samples = [
            (float(alt), float(az), float(distance_m))
            for (alt, az), distance_m in zip(layer, layer_distances)
        ]
        if len(layer_samples) >= 3:
            layer_samples = _rotate_profile_to_seam_azimuth(
                layer_samples,
                seam_az_deg=seam_az_deg,
            )

        base_width = _distance_band_widths(
            distance_km=representative_distance_km,
            band_count=layer_count,
        )
        underlay_width = _distance_band_underlay_width(
            distance_km=representative_distance_km,
            band_count=layer_count,
        )
        band_alpha = _distance_band_alpha(
            distance_km=alpha_distance_km,
            band_count=layer_count,
            opacity=layer_opacity,
        )
        underlay_alpha = _distance_band_underlay_alpha(
            distance_km=alpha_distance_km,
            band_count=layer_count,
            opacity=layer_opacity,
        )
        visible_points: list[tuple[float, float]] = []
        visible_altaz: list[tuple[float, float]] = []
        for alt, az, _distance_m in layer_samples:
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
            visible_points.append((nx, ny))
            visible_altaz.append((float(alt), float(az)))

        if len(visible_points) < 2:
            continue

        point_fragments = split_by_gaps_func(visible_points) if len(visible_points) > 2 else [visible_points]
        point_offset = 0
        for frag in point_fragments:
            if len(frag) < 2:
                point_offset += len(frag)
                continue
            frag_points = [QPointF(*normalized_to_screen_xy_func(nx, ny, geometry)) for nx, ny in frag]
            frag_altaz = visible_altaz[point_offset:point_offset + len(frag)]
            point_offset += len(frag)
            run_points: list[QPointF] = [frag_points[0]]
            run_is_hidden = False
            for start_idx, (start_point, end_point) in enumerate(zip(frag_points, frag_points[1:])):
                start_alt, start_az = frag_altaz[start_idx]
                end_alt, end_az = frag_altaz[start_idx + 1]
                segment_mid_alt = 0.5 * (float(start_alt) + float(end_alt))
                segment_mid_az = _circular_midpoint_azimuth_deg(float(start_az), float(end_az))
                bin_key = _azimuth_bin_key(segment_mid_az, bin_size_deg=occlusion_bin_deg)
                previous_max_alt = max_visible_alt_by_bin.get(bin_key, float("-inf"))
                segment_is_hidden = math.isfinite(previous_max_alt) and (
                    float(previous_max_alt) - float(segment_mid_alt)
                ) >= hidden_altitude_delta_deg
                if start_idx == 0:
                    run_is_hidden = segment_is_hidden
                if segment_is_hidden != run_is_hidden:
                    _draw_secondary_ridge_run(
                        points=run_points,
                        is_hidden=run_is_hidden,
                        base_width=base_width,
                        underlay_width=underlay_width,
                        band_alpha=band_alpha,
                        underlay_alpha=underlay_alpha,
                    )
                    run_points = [start_point, end_point]
                    run_is_hidden = segment_is_hidden
                else:
                    run_points.append(end_point)
                max_visible_alt_by_bin[bin_key] = max(previous_max_alt, float(segment_mid_alt))
            _draw_secondary_ridge_run(
                points=run_points,
                is_hidden=run_is_hidden,
                base_width=base_width,
                underlay_width=underlay_width,
                band_alpha=band_alpha,
                underlay_alpha=underlay_alpha,
            )


def draw_urban_outlines(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    urban_outlines: list[UrbanOutlinePolyline] | None,
    *,
    opacity: float = 0.2,
    fill_opacity_factor: float = 1.0,
    line_width_scale: float = 1.0,
    layer_style: OverlayLayerStyle | None = None,
    is_in_fov_func: Callable[..., bool] = is_in_fov,
    altaz_to_normalized_xy_func: Callable[[float, float, tuple[float, float]], tuple[float, float]] = altaz_to_normalized_xy,
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], tuple[float, float]] = normalized_to_screen_xy,
    split_by_gaps_func: Callable[[list[tuple[float, float]]], list[list[tuple[float, float]]]] = split_by_gaps,
) -> None:
    """Draw sampled building-top outlines directly on the sky dome."""
    if not urban_outlines:
        return
    alpha_scale = 1.0 if layer_style is None else float(layer_style.alpha_scale)
    layer_opacity = float(opacity) * alpha_scale
    fill_factor = max(0.0, min(1.0, float(fill_opacity_factor)))
    if layer_opacity <= 0.0:
        return

    view_center, edge_fov_deg, content_fov_deg = _viewer_projection_params(viewer)
    painter.save()
    width_scale = float(line_width_scale) * (
        1.0 if layer_style is None else float(layer_style.width_scale)
    )
    layer_rgb = URBAN_OUTLINE_LAYER_LINE_COLOR if layer_style is None else layer_style.rgb
    underlay_rgb = (
        layer_rgb
        if layer_style is None or layer_style.outline_rgba is None
        else layer_style.outline_rgba[:3]
    )
    for outline_entry in urban_outlines:
        outline = list(outline_entry.points)
        distance_km = float(getattr(outline_entry, "distance_km", float("inf")))
        if len(outline) < 2:
            continue
        height_scale = _urban_outline_height_width_scale(float(getattr(outline_entry, "height_m", 0.0)))
        thickened_width_scale = width_scale * height_scale
        fill_color = QColor(*layer_rgb)
        fill_color.setAlpha(
            int(round(_urban_outline_fill_alpha(layer_opacity) * fill_factor))
        )
        foreground_width = _urban_outline_foreground_width(distance_km, width_scale=width_scale)
        foreground_color = QColor(*layer_rgb)
        foreground_color.setAlpha(
            int(
                round(
                    255.0
                    * _dampen_alpha_for_narrow_width(
                        _urban_outline_foreground_alpha(layer_opacity),
                        foreground_width,
                    )
                )
            )
        )
        foreground_pen = QPen(
            foreground_color,
            foreground_width,
            Qt.PenStyle.SolidLine,
        )
        foreground_pen.setCosmetic(True)
        foreground_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        foreground_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)

        underlay_pen = None
        mid_underlay_pen = None
        outer_underlay_pen = None
        if _urban_outline_uses_underlay(distance_km):
            underlay_width = _urban_outline_underlay_width(distance_km, width_scale=thickened_width_scale)
            underlay_color = QColor(*underlay_rgb)
            underlay_color.setAlpha(
                int(
                    round(
                        255.0
                        * _dampen_alpha_for_narrow_width(
                            0.25 * _urban_outline_underlay_alpha(layer_opacity),
                            underlay_width,
                        )
                    )
                )
            )
            underlay_pen = QPen(
                underlay_color,
                underlay_width,
                Qt.PenStyle.SolidLine,
            )
            underlay_pen.setCosmetic(True)
            underlay_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            underlay_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            mid_underlay_width = _urban_outline_mid_width(distance_km, width_scale=thickened_width_scale)
            mid_underlay_color = QColor(*underlay_rgb)
            mid_underlay_color.setAlpha(
                int(
                    round(
                        255.0
                        * _dampen_alpha_for_narrow_width(
                            0.50 * _urban_outline_underlay_alpha(layer_opacity),
                            mid_underlay_width,
                        )
                    )
                )
            )
            mid_underlay_pen = QPen(
                mid_underlay_color,
                mid_underlay_width,
                Qt.PenStyle.SolidLine,
            )
            mid_underlay_pen.setCosmetic(True)
            mid_underlay_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            mid_underlay_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            outer_underlay_width = _urban_outline_outer_width(distance_km, width_scale=thickened_width_scale)
            outer_underlay_color = QColor(*underlay_rgb)
            outer_underlay_color.setAlpha(
                int(
                    round(
                        255.0
                        * _dampen_alpha_for_narrow_width(
                            0.90 * _urban_outline_underlay_alpha(layer_opacity),
                            outer_underlay_width,
                        )
                    )
                )
            )
            outer_underlay_pen = QPen(
                outer_underlay_color,
                outer_underlay_width,
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

        def _fill_fragments(
            fragments: list[list[tuple[float, float]]],
            fill_color: QColor = fill_color,
        ) -> None:
            if fill_color.alpha() <= 0:
                return
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(fill_color)
            for frag in fragments:
                if len(frag) < 3:
                    continue
                screen_points = [
                    QPointF(*normalized_to_screen_xy_func(nx, ny, geometry))
                    for nx, ny in frag
                ]
                if not _urban_outline_fragment_is_closed_for_fill(screen_points):
                    continue
                if not _urban_outline_fragment_is_large_enough_for_fill(screen_points):
                    continue
                polygon = QPolygonF(screen_points)
                if polygon.isEmpty():
                    continue
                painter.drawPolygon(polygon)

        def _draw_points(
            points: list[tuple[float, float]],
            outer_underlay_pen: QPen | None = outer_underlay_pen,
            mid_underlay_pen: QPen | None = mid_underlay_pen,
            underlay_pen: QPen | None = underlay_pen,
            foreground_pen: QPen = foreground_pen,
        ) -> None:
            if len(points) < 2:
                return
            if len(points) == 2:
                clip_radius = float(content_fov_deg) / max(1.0e-6, float(edge_fov_deg))
                fragments = _clip_polyline_to_radius(points, clip_radius)
            else:
                fragments = split_by_gaps_func(points)
            _fill_fragments(fragments)
            if outer_underlay_pen is not None:
                _draw_fragments(fragments, outer_underlay_pen)
            if mid_underlay_pen is not None:
                _draw_fragments(fragments, mid_underlay_pen)
            if underlay_pen is not None:
                _draw_fragments(fragments, underlay_pen)
            _draw_fragments(fragments, foreground_pen)

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


def draw_water_overlay_dots(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    water_dots: list[WaterOverlayPoint] | None,
    *,
    opacity: float = 0.85,
    line_width_scale: float = 1.0,
    layer_style: OverlayLayerStyle | None = None,
    fast_mode: bool = False,
    pairwise_thinning: bool = True,
    is_in_fov_func: Callable[..., bool] = is_in_fov,
    altaz_to_normalized_xy_func: Callable[[float, float, tuple[float, float]], tuple[float, float]] = altaz_to_normalized_xy,
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], tuple[float, float]] = normalized_to_screen_xy,
) -> None:
    """Draw sampled water surface points as small filled circles."""
    alpha_scale = 1.0 if layer_style is None else float(layer_style.alpha_scale)
    layer_opacity = max(0.0, min(1.0, float(opacity) * alpha_scale))
    if not water_dots or layer_opacity <= 0.0:
        return

    view_center, edge_fov_deg, content_fov_deg = _viewer_projection_params(viewer)
    dots_to_draw = (
        _thin_water_overlay_dots_pairwise(water_dots) if pairwise_thinning else list(water_dots)
    )
    visible_points = _visible_water_overlay_dots(
        dots_to_draw,
        view_center=view_center,
        content_fov_deg=float(content_fov_deg),
        is_in_fov_func=is_in_fov_func,
    )
    if not visible_points:
        return
    dot_alpha = max(0, min(255, int(round(255.0 * layer_opacity))))

    painter.save()
    for point in visible_points:
        alt = float(point.alt_deg)
        az = float(point.az_deg)
        try:
            nx, ny = altaz_to_normalized_xy_func(
                alt,
                az,
                view_center,
                edge_fov_deg=float(edge_fov_deg),
            )
        except TypeError:
            nx, ny = _project_altaz_to_normalized_xy(
                alt,
                az,
                view_center,
                edge_fov_deg=float(edge_fov_deg),
            )
        px, py = normalized_to_screen_xy_func(nx, ny, geometry)
        scan_distance_m = float(getattr(point, "scan_distance_m", 0.0) or 0.0)
        distance_alpha = _water_overlay_distance_alpha_scale(scan_distance_m)
        major_radius, _minor_radius, _pen_width = _water_overlay_marker_geometry(
            line_width_scale,
            distance_m=scan_distance_m,
        )
        point_alpha = max(0, min(255, int(round(dot_alpha * distance_alpha))))
        point_rgb = _water_overlay_point_color_rgb(point) if layer_style is None else layer_style.rgb
        outline_color = QColor(*point_rgb, point_alpha)
        painter.save()
        painter.translate(float(px), float(py))
        if layer_style is not None and layer_style.outline_rgba is not None:
            underlay_color = QColor(*layer_style.outline_rgba)
            underlay_color.setAlpha(
                max(0, min(255, int(round(underlay_color.alpha() * distance_alpha * layer_opacity))))
            )
            underlay_pen = QPen(underlay_color)
            underlay_pen.setWidthF(float(_pen_width) + 2.0)
            underlay_pen.setCosmetic(True)
            underlay_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            underlay_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(underlay_pen)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawEllipse(QPointF(0.0, 0.0), major_radius, _minor_radius)
        pen = QPen(outline_color)
        pen.setWidthF(float(_pen_width))
        pen.setCosmetic(True)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        painter.setPen(pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawEllipse(QPointF(0.0, 0.0), major_radius, _minor_radius)
        painter.restore()
    painter.restore()


def draw_water_overlay_polylines(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    water_polylines: list[WaterOverlayPolyline] | None,
    *,
    opacity: float = 0.65,
    line_width_scale: float = 1.0,
    layer_style: OverlayLayerStyle | None = None,
    is_in_fov_func: Callable[..., bool] = is_in_fov,
    altaz_to_normalized_xy_func: Callable[[float, float, tuple[float, float]], tuple[float, float]] = altaz_to_normalized_xy,
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], tuple[float, float]] = normalized_to_screen_xy,
) -> None:
    """Draw simplified water footprint rings as clipped screen polylines."""
    if not water_polylines or opacity <= 0.0:
        return
    alpha_scale = 1.0 if layer_style is None else float(layer_style.alpha_scale)
    layer_opacity = max(0.0, min(1.0, float(opacity) * alpha_scale))
    view_center, edge_fov_deg, content_fov_deg = _viewer_projection_params(viewer)
    color_rgb = WATER_OVERLAY_POINT_COLOR_RGB if layer_style is None else layer_style.rgb
    painter.save()
    for polyline in water_polylines:
        screen_points: list[QPointF] = []
        for point in polyline.points:
            if not is_in_fov_func(
                float(point.alt_deg),
                float(point.az_deg),
                view_center,
                fov_deg=float(content_fov_deg),
            ):
                if len(screen_points) >= 2:
                    pen = QPen(QColor(*color_rgb, int(round(255.0 * layer_opacity))))
                    pen.setWidthF(max(1.0, 1.35 * float(line_width_scale)))
                    pen.setCosmetic(True)
                    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
                    painter.setPen(pen)
                    painter.setBrush(Qt.BrushStyle.NoBrush)
                    painter.drawPolyline(QPolygonF(screen_points))
                screen_points = []
                continue
            nx, ny = altaz_to_normalized_xy_func(
                float(point.alt_deg),
                float(point.az_deg),
                view_center,
                edge_fov_deg=float(edge_fov_deg),
            )
            px, py = normalized_to_screen_xy_func(nx, ny, geometry)
            screen_points.append(QPointF(float(px), float(py)))
        if len(screen_points) >= 2:
            pen = QPen(QColor(*color_rgb, int(round(255.0 * layer_opacity))))
            pen.setWidthF(max(1.0, 1.35 * float(line_width_scale)))
            pen.setCosmetic(True)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(pen)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPolyline(QPolygonF(screen_points))
    painter.restore()


def _visible_water_overlay_dots(
    water_dots: list[WaterOverlayPoint],
    *,
    view_center: tuple[float, float],
    content_fov_deg: float,
    is_in_fov_func: Callable[..., bool],
) -> list[WaterOverlayPoint]:
    visible: list[WaterOverlayPoint] = []
    for point in water_dots:
        alt = float(point.alt_deg)
        az = float(point.az_deg)
        if is_in_fov_func(alt, az, view_center, fov_deg=content_fov_deg):
            visible.append(point)
    return visible


def _thin_water_overlay_dots_pairwise(
    water_dots: list[WaterOverlayPoint],
) -> list[WaterOverlayPoint]:
    grouped: dict[tuple[int, int], list[WaterOverlayPoint]] = {}
    fallback: list[WaterOverlayPoint] = []
    for point in water_dots:
        azimuth_index = getattr(point, "scan_azimuth_index", None)
        distance_index = getattr(point, "scan_distance_index", None)
        if not isinstance(azimuth_index, int) or not isinstance(distance_index, int):
            fallback.append(point)
            continue
        group_key = (int(azimuth_index), int(distance_index) // 2)
        grouped.setdefault(group_key, []).append(point)
    if not grouped:
        return fallback
    ordered_grouped: list[WaterOverlayPoint] = []
    for (azimuth_index, pair_index), points in sorted(grouped.items()):
        preferred_parity = azimuth_index % 2
        def _scan_distance_parity(point: WaterOverlayPoint) -> int:
            distance_index = getattr(point, "scan_distance_index", None)
            if not isinstance(distance_index, int):
                return -1
            return distance_index % 2

        chosen = next(
            (point for point in points if _scan_distance_parity(point) == preferred_parity),
            None,
        )
        if chosen is None:
            chosen = min(
                points,
                key=lambda item: int(getattr(item, "scan_distance_index", 0) or 0),
            )
        ordered_grouped.append(chosen)
    if fallback:
        ordered_grouped.extend(fallback)
    return ordered_grouped
