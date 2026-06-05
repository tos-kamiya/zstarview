import math
from dataclasses import dataclass
from typing import Callable, List, Tuple

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import (
    PALETTE_ASTERISM_RGB,
    TERRAIN_HORIZON_LINE_COLOR,
    URBAN_OUTLINE_LAYER_LINE_COLOR,
)
from ..types import ScreenGeometry, UrbanOutlinePolyline, ViewerData
from ..water_overlay import WaterOverlayPoint
from .geometry import normalized_to_screen_xy
from .guides import _clip_polyline_to_radius, split_by_gaps

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
TERRAIN_SECONDARY_RIDGE_VISIBLE_COLOR_RGB = PALETTE_ASTERISM_RGB
TERRAIN_SECONDARY_RIDGE_OCCLUDED_COLOR_RGB = TERRAIN_HORIZON_LINE_COLOR
TERRAIN_SECONDARY_RIDGE_OCCLUSION_BIN_DEG = 1.0
TERRAIN_SECONDARY_RIDGE_OCCLUSION_EPSILON_DEG = 0.05
TERRAIN_SECONDARY_RIDGE_SEAM_BRIDGE_SCREEN_GAP = 0.25
TERRAIN_SECONDARY_RIDGE_SEAM_BRIDGE_AZ_GAP_DEG = 4.0
TERRAIN_SECONDARY_RIDGE_GLOW_OUTER_WIDTH_SCALE = 2.05
TERRAIN_SECONDARY_RIDGE_GLOW_MID_WIDTH_SCALE = 1.35
TERRAIN_SECONDARY_RIDGE_GLOW_OUTER_ALPHA_SCALE = 0.06
TERRAIN_SECONDARY_RIDGE_GLOW_MID_ALPHA_SCALE = 0.18
WATER_OVERLAY_POINT_COLOR_RGB = (122, 218, 240)
# Keep sea tones aligned with dev-samples/basic-color-palette.html:
# HSV(191°, 49.2%, 94.1%) -> RGB(122, 218, 240).
WATER_OVERLAY_SEA_COLOR_RGB = (122, 218, 240)
WATER_OVERLAY_SEA_125_COLOR_RGB = (122, 218, 240)
WATER_OVERLAY_SEA_250_COLOR_RGB = (54, 200, 184)
WATER_OVERLAY_SEA_500_COLOR_RGB = (255, 170, 64)
WATER_OVERLAY_LAKE_COLOR_RGB = (104, 196, 168)
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


def _azimuth_bin_key(az_deg: float) -> int:
    bin_size_deg = TERRAIN_SECONDARY_RIDGE_OCCLUSION_BIN_DEG
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


def _circular_azimuth_delta_deg(start_az_deg: float, end_az_deg: float) -> float:
    delta = (float(end_az_deg) - float(start_az_deg)) % 360.0
    return min(delta, 360.0 - delta)


def _straddles_seam_azimuth(start_az_deg: float, end_az_deg: float, seam_az_deg: float) -> bool:
    seam = float(seam_az_deg) % 360.0
    start_rel = (float(start_az_deg) - seam) % 360.0
    end_rel = (float(end_az_deg) - seam) % 360.0
    return abs(start_rel - end_rel) > 180.0


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


def _terrain_secondary_ridge_glow_pass_specs(
    visible_width: float,
    base_alpha: float,
) -> tuple[tuple[float, float], ...]:
    core_width = max(0.0, float(visible_width))
    return (
        (
            max(core_width, core_width * TERRAIN_SECONDARY_RIDGE_GLOW_OUTER_WIDTH_SCALE),
            _dampen_alpha_for_narrow_width(
                float(base_alpha) * TERRAIN_SECONDARY_RIDGE_GLOW_OUTER_ALPHA_SCALE,
                core_width,
            ),
        ),
        (
            max(core_width, core_width * TERRAIN_SECONDARY_RIDGE_GLOW_MID_WIDTH_SCALE),
            _dampen_alpha_for_narrow_width(
                float(base_alpha) * TERRAIN_SECONDARY_RIDGE_GLOW_MID_ALPHA_SCALE,
                core_width,
            ),
        ),
        (
            core_width,
            _dampen_alpha_for_narrow_width(float(base_alpha), core_width),
        ),
    )


def _draw_terrain_secondary_ridge_glow(
    painter: QPainter,
    start_point: QPointF,
    end_point: QPointF,
    *,
    visible_width: float,
    visible_alpha: float,
    line_width_scale: float,
) -> None:
    for width_scale, alpha_scale in _terrain_secondary_ridge_glow_pass_specs(
        visible_width,
        visible_alpha,
    ):
        painter.setPen(
            _solid_pen(
                TERRAIN_SECONDARY_RIDGE_VISIBLE_COLOR_RGB,
                alpha_scale,
                width_scale * float(line_width_scale),
            )
        )
        painter.drawLine(start_point, end_point)


def _draw_terrain_profile_layer(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    terrain_profile_distances_m: list[float] | None,
    *,
    spec: TerrainHorizonRenderSpec,
    is_in_fov_func: Callable[..., bool],
    altaz_to_normalized_xy_func: Callable[[float, float, Tuple[float, float]], Tuple[float, float]],
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], Tuple[float, float]],
    split_by_gaps_func: Callable[[List[Tuple[float, float]]], List[List[Tuple[float, float]]]],
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


def draw_terrain_secondary_ridges(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    terrain_secondary_ridges_layers: list[list[tuple[float, float]]] | None,
    terrain_secondary_ridges_distances_m_layers: list[list[float]] | None,
    *,
    opacity: float = 0.25,
    line_width_scale: float = 1.0,
    is_in_fov_func: Callable[..., bool] = is_in_fov,
    altaz_to_normalized_xy_func: Callable[[float, float, Tuple[float, float]], Tuple[float, float]] = altaz_to_normalized_xy,
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], Tuple[float, float]] = normalized_to_screen_xy,
    split_by_gaps_func: Callable[[List[Tuple[float, float]]], List[List[Tuple[float, float]]]] = split_by_gaps,
) -> None:
    """Draw fixed-width ridge bands grouped by distance interval."""
    if opacity <= 0.0:
        return

    if not terrain_secondary_ridges_layers:
        return

    ridge_opacity = terrain_secondary_ridge_line_alpha(opacity)
    if ridge_opacity <= 0.0:
        return

    if terrain_secondary_ridges_distances_m_layers is not None and len(terrain_secondary_ridges_distances_m_layers) != len(terrain_secondary_ridges_layers):
        terrain_secondary_ridges_distances_m_layers = None
    layer_count = len(terrain_secondary_ridges_layers)
    overlay_scale = 1.7
    overlay_alpha_scale = 0.2
    max_visible_alt_by_bin: dict[int, float] = {}
    view_center = tuple(float(value) for value in viewer.view_center)
    edge_fov_deg = float(viewer.edge_fov_deg)
    content_fov_deg = float(viewer.content_fov_deg)
    seam_az_deg = (float(view_center[1]) + 180.0) % 360.0
    alpha_distance_km = TERRAIN_DISTANCE_BAND_REFERENCE_ALPHA_DISTANCE_KM

    def _can_bridge_seam(
        start_point: QPointF,
        start_altaz: tuple[float, float],
        end_point: QPointF,
        end_altaz: tuple[float, float],
    ) -> bool:
        bridge_az_gap = _circular_azimuth_delta_deg(start_altaz[1], end_altaz[1])
        bridge_straddles_seam = _straddles_seam_azimuth(
            float(start_altaz[1]),
            float(end_altaz[1]),
            seam_az_deg,
        )
        bridge_screen_gap = math.hypot(
            float(end_point.x()) - float(start_point.x()),
            float(end_point.y()) - float(start_point.y()),
        )
        return (
            bridge_az_gap <= TERRAIN_SECONDARY_RIDGE_SEAM_BRIDGE_AZ_GAP_DEG
            and bridge_straddles_seam
            and bridge_screen_gap <= TERRAIN_SECONDARY_RIDGE_SEAM_BRIDGE_SCREEN_GAP
        )

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
            opacity=opacity,
        )
        underlay_alpha = _distance_band_underlay_alpha(
            distance_km=alpha_distance_km,
            band_count=layer_count,
            opacity=opacity,
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
        visible_bridge_start: tuple[QPointF, tuple[float, float]] | None = None
        visible_bridge_end: tuple[QPointF, tuple[float, float]] | None = None
        for frag in point_fragments:
            if len(frag) < 2:
                point_offset += len(frag)
                continue
            frag_points = [QPointF(*normalized_to_screen_xy_func(nx, ny, geometry)) for nx, ny in frag]
            frag_altaz = visible_altaz[point_offset:point_offset + len(frag)]
            point_offset += len(frag)
            if visible_bridge_start is None:
                visible_bridge_start = (frag_points[0], frag_altaz[0])
            visible_bridge_end = (frag_points[-1], frag_altaz[-1])
            poly = QPolygonF(frag_points)
            if underlay_alpha > 0.0 and underlay_width > base_width:
                painter.setPen(
                    _solid_pen(
                        TERRAIN_SECONDARY_RIDGE_OCCLUDED_COLOR_RGB,
                        underlay_alpha,
                        float(underlay_width) * float(line_width_scale),
                    )
                )
                painter.drawPolyline(poly)
                painter.setPen(
                    _solid_pen(
                        TERRAIN_SECONDARY_RIDGE_OCCLUDED_COLOR_RGB,
                        band_alpha,
                        float(base_width) * float(line_width_scale),
                    )
                )
                painter.drawPolyline(poly)
            for start_idx, (start_point, end_point) in enumerate(zip(frag_points, frag_points[1:])):
                start_alt, start_az = frag_altaz[start_idx]
                end_alt, end_az = frag_altaz[start_idx + 1]
                segment_mid_alt = 0.5 * (float(start_alt) + float(end_alt))
                segment_mid_az = _circular_midpoint_azimuth_deg(float(start_az), float(end_az))
                bin_key = _azimuth_bin_key(segment_mid_az)
                is_occluded = segment_mid_alt <= (
                    max_visible_alt_by_bin.get(bin_key, float("-inf"))
                    - TERRAIN_SECONDARY_RIDGE_OCCLUSION_EPSILON_DEG
                )
                if is_occluded:
                    continue
                max_visible_alt_by_bin[bin_key] = max(
                    max_visible_alt_by_bin.get(bin_key, float("-inf")),
                    segment_mid_alt,
                )
                _draw_terrain_secondary_ridge_glow(
                    painter,
                    start_point,
                    end_point,
                    visible_width=max(0.0, base_width * overlay_scale),
                    visible_alpha=band_alpha * overlay_alpha_scale,
                    line_width_scale=line_width_scale,
                )

        if visible_bridge_start is not None and visible_bridge_end is not None:
            bridge_start_point, bridge_start_altaz = visible_bridge_start
            bridge_end_point, bridge_end_altaz = visible_bridge_end
            if _can_bridge_seam(
                bridge_start_point,
                bridge_start_altaz,
                bridge_end_point,
                bridge_end_altaz,
            ):
                bridge_mid_alt = 0.5 * (float(bridge_start_altaz[0]) + float(bridge_end_altaz[0]))
                bridge_mid_az = _circular_midpoint_azimuth_deg(
                    float(bridge_start_altaz[1]),
                    float(bridge_end_altaz[1]),
                )
                bridge_bin_key = _azimuth_bin_key(bridge_mid_az)
                if bridge_mid_alt > (
                    max_visible_alt_by_bin.get(bridge_bin_key, float("-inf"))
                    - TERRAIN_SECONDARY_RIDGE_OCCLUSION_EPSILON_DEG
                ):
                    _draw_terrain_secondary_ridge_glow(
                        painter,
                        bridge_start_point,
                        bridge_end_point,
                        visible_width=max(0.0, base_width * overlay_scale),
                        visible_alpha=band_alpha * overlay_alpha_scale,
                        line_width_scale=line_width_scale,
                    )


def draw_urban_outlines(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    urban_outlines: list[UrbanOutlinePolyline] | None,
    *,
    opacity: float = 0.2,
    line_width_scale: float = 1.0,
    is_in_fov_func: Callable[..., bool] = is_in_fov,
    altaz_to_normalized_xy_func: Callable[[float, float, Tuple[float, float]], Tuple[float, float]] = altaz_to_normalized_xy,
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], Tuple[float, float]] = normalized_to_screen_xy,
    split_by_gaps_func: Callable[[List[Tuple[float, float]]], List[List[Tuple[float, float]]]] = split_by_gaps,
) -> None:
    """Draw sampled building-top outlines directly on the sky dome."""
    if not urban_outlines:
        return
    if float(opacity) <= 0.0:
        return

    view_center, edge_fov_deg, content_fov_deg = _viewer_projection_params(viewer)
    painter.save()
    width_scale = float(line_width_scale)
    for outline_entry in urban_outlines:
        outline = list(outline_entry.points)
        distance_km = float(getattr(outline_entry, "distance_km", float("inf")))
        if len(outline) < 2:
            continue
        height_scale = _urban_outline_height_width_scale(float(getattr(outline_entry, "height_m", 0.0)))
        thickened_width_scale = width_scale * height_scale
        foreground_width = _urban_outline_foreground_width(distance_km, width_scale=width_scale)
        foreground_color = QColor(*URBAN_OUTLINE_LAYER_LINE_COLOR)
        foreground_color.setAlpha(
            int(
                round(
                    255.0
                    * _dampen_alpha_for_narrow_width(
                        _urban_outline_foreground_alpha(opacity),
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
            underlay_color = QColor(*URBAN_OUTLINE_LAYER_LINE_COLOR)
            underlay_color.setAlpha(
                int(
                    round(
                        255.0
                        * _dampen_alpha_for_narrow_width(
                            0.25 * _urban_outline_underlay_alpha(opacity),
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
            mid_underlay_color = QColor(*URBAN_OUTLINE_LAYER_LINE_COLOR)
            mid_underlay_color.setAlpha(
                int(
                    round(
                        255.0
                        * _dampen_alpha_for_narrow_width(
                            0.50 * _urban_outline_underlay_alpha(opacity),
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
            outer_underlay_color = QColor(*URBAN_OUTLINE_LAYER_LINE_COLOR)
            outer_underlay_color.setAlpha(
                int(
                    round(
                        255.0
                        * _dampen_alpha_for_narrow_width(
                            0.90 * _urban_outline_underlay_alpha(opacity),
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

        def _draw_points(points: list[tuple[float, float]]) -> None:
            if len(points) < 2:
                return
            if len(points) == 2:
                clip_radius = float(content_fov_deg) / max(1.0e-6, float(edge_fov_deg))
                fragments = _clip_polyline_to_radius(points, clip_radius)
            else:
                fragments = split_by_gaps_func(points)
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
    fast_mode: bool = False,
    pairwise_thinning: bool = True,
    is_in_fov_func: Callable[..., bool] = is_in_fov,
    altaz_to_normalized_xy_func: Callable[[float, float, tuple[float, float]], tuple[float, float]] = altaz_to_normalized_xy,
    normalized_to_screen_xy_func: Callable[[float, float, ScreenGeometry], tuple[float, float]] = normalized_to_screen_xy,
) -> None:
    """Draw sampled water surface points as small filled circles."""
    layer_opacity = max(0.0, min(1.0, float(opacity)))
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
        outline_color = QColor(
            *_water_overlay_point_color_rgb(point),
            point_alpha,
        )
        painter.save()
        painter.translate(float(px), float(py))
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
