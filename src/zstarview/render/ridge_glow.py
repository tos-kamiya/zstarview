from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QBrush, QColor, QImage, QPainter, QPen, QPixmap, QPolygonF

from ..astro import altaz_to_normalized_xy
from ..types import ScreenGeometry, ViewerData
from .background import dimalt_ring_brightness_score_from_rgba
from .geometry import normalized_to_screen_xy
from .guides import split_by_gaps
from .night_lights import NIGHT_LIGHTS_GLOW_RGB, NightLightGlowProfile, night_light_strength_factor

RIDGE_GLOW_MIN_BRIGHTNESS = 0.02
RIDGE_GLOW_FILL_RGB = (244, 246, 248)
RIDGE_GLOW_SKY_RGB = (255, 170, 48)
RIDGE_GLOW_SKY_ALPHA_FLOOR = 0.02
RIDGE_GLOW_SKY_ALPHA_SCALE = 0.85


@dataclass(frozen=True)
class RidgeGlowLayerSpec:
    upper_alt_offset_deg: float
    alpha_scale: float
    window_size: int
    focus_scale: float


RIDGE_GLOW_SKY_LAYER_SPECS = [
    RidgeGlowLayerSpec(upper_alt_offset_deg=0.3, alpha_scale=1.0, window_size=2, focus_scale=1.8),
    RidgeGlowLayerSpec(upper_alt_offset_deg=0.7, alpha_scale=0.5, window_size=4, focus_scale=1.6),
    RidgeGlowLayerSpec(upper_alt_offset_deg=1.6, alpha_scale=0.25, window_size=8, focus_scale=1.4),
    RidgeGlowLayerSpec(upper_alt_offset_deg=3.5, alpha_scale=0.125, window_size=16, focus_scale=1.2),
    RidgeGlowLayerSpec(upper_alt_offset_deg=7.4, alpha_scale=0.12, window_size=32, focus_scale=1.1),
]


def _ridge_glow_directional_altitudes(
    raw_altitudes: list[float],
    window_size: int,
    *,
    center_weight_scale: float = 1.0,
) -> list[float]:
    if len(raw_altitudes) < 2:
        return list(raw_altitudes)
    window_size = max(1, int(window_size))
    if window_size == 1:
        return list(raw_altitudes)

    raw_altitudes_array = np.asarray([float(value) for value in raw_altitudes], dtype=np.float64)
    radius = max(1, int(window_size) // 2)
    smoothed: list[float] = []
    focus_scale = max(0.1, float(center_weight_scale))
    for index in range(raw_altitudes_array.size):
        start = max(0, index - radius)
        end = min(raw_altitudes_array.size, index + radius + 1)
        window = raw_altitudes_array[start:end]
        offsets = np.arange(start, end, dtype=np.float64) - float(index)
        weights = np.power(np.clip((float(radius) + 1.0) - np.abs(offsets), 0.0, None), focus_scale)
        weights = np.clip(weights, 0.0, None)
        weight_sum = float(np.sum(weights))
        if weight_sum <= 0.0:
            smoothed.append(float(raw_altitudes_array[index]))
            continue
        smoothed.append(float(np.sum(window * (weights / weight_sum))))
    return smoothed


def _seam_relative_azimuth_deg(azimuth_deg: float, seam_az_deg: float) -> float:
    return (float(azimuth_deg) - float(seam_az_deg)) % 360.0


def _layer_distance_km(profile: NightLightGlowProfile, layer_index: int) -> float:
    band_profiles = tuple(getattr(profile, "band_profiles", ()))
    assert band_profiles, "profile.band_profiles must not be empty"
    assert 0 <= int(layer_index) < len(band_profiles), "layer_index out of range"
    return float(band_profiles[int(layer_index)].max_distance_km)


def _band_lower_edge_altitudes(
    current_layer_altitudes: np.ndarray,
    current_layer_azimuths_deg: np.ndarray,
    previous_layer_altaz: list[tuple[float, float]] | None,
) -> np.ndarray:
    current = np.asarray(current_layer_altitudes, dtype=np.float64)
    current_azimuths = np.asarray(current_layer_azimuths_deg, dtype=np.float64)
    if current.shape != current_azimuths.shape:
        raise ValueError("current_layer_azimuths_deg must match current_layer_altitudes")
    if not previous_layer_altaz:
        return current.copy()

    previous_samples = sorted(
        [
            (
                float(alt_deg),
                float(az_deg) % 360.0,
            )
            for alt_deg, az_deg in previous_layer_altaz
        ],
        key=lambda item: item[1],
    )
    previous_azimuths = np.asarray([az for _, az in previous_samples], dtype=np.float64)
    previous_altitudes = np.asarray([alt for alt, _ in previous_samples], dtype=np.float64)
    if previous_azimuths.size == 0:
        return current.copy()
    if previous_azimuths.size == 1:
        return np.full_like(current, float(previous_altitudes[0]))

    previous_azimuths_ext = np.concatenate(
        [
            previous_azimuths[-1:] - 360.0,
            previous_azimuths,
            previous_azimuths[:1] + 360.0,
        ]
    )
    previous_altitudes_ext = np.concatenate(
        [
            previous_altitudes[-1:],
            previous_altitudes,
            previous_altitudes[:1],
        ]
    )
    return np.interp(np.asarray(current_azimuths, dtype=np.float64) % 360.0, previous_azimuths_ext, previous_altitudes_ext)


def _ridge_glow_hatch_brush(color_rgb: tuple[int, int, int], opacity_scale: float) -> QBrush:
    tile_size = 8
    tile = QImage(tile_size, tile_size, QImage.Format.Format_ARGB32_Premultiplied)
    tile.fill(Qt.GlobalColor.transparent)

    painter = QPainter(tile)
    try:
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        base_color = QColor(*color_rgb)
        base_color.setAlpha(max(0, min(255, int(round(32 * float(opacity_scale))))))
        painter.fillRect(tile.rect(), base_color)

        hatch_color = QColor(*color_rgb)
        hatch_color.setAlpha(max(0, min(255, int(round(255 * float(opacity_scale))))))
        hatch_pen = QPen(hatch_color, 1)
        hatch_pen.setCosmetic(True)
        painter.setPen(hatch_pen)
        painter.drawLine(0, tile_size - 1, tile_size - 1, 0)
    finally:
        painter.end()

    return QBrush(QPixmap.fromImage(tile))


def _draw_ridge_glow_fragments(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    seam_az_deg: float,
    projected_points: list[tuple[float, float]],
    projected_draw_az: list[float],
    projected_alts: list[float],
    strengths: list[float],
    color_rgb: tuple[int, int, int],
    opacity_scale: float,
    band_thickness_deg: float,
) -> None:
    layer_specs = tuple(RIDGE_GLOW_SKY_LAYER_SPECS)
    glow_rgb = tuple(int(value) for value in color_rgb)
    point_index = 0
    for fragment in split_by_gaps(projected_points):
        if len(fragment) < 2:
            point_index += len(fragment)
            continue

        frag_strengths = strengths[point_index : point_index + len(fragment)]
        frag_strength = max(float(value) for value in frag_strengths)
        if frag_strength < RIDGE_GLOW_MIN_BRIGHTNESS:
            point_index += len(fragment)
            continue

        lower_points: list[QPointF] = []
        frag_altaz = zip(
            projected_draw_az[point_index : point_index + len(fragment)],
            projected_alts[point_index : point_index + len(fragment)],
        )
        lower_altitudes: list[float] = []
        lower_azimuths: list[float] = []
        for seam_az_deg_value, lower_alt in frag_altaz:
            az = (float(seam_az_deg_value) + seam_az_deg) % 360.0
            try:
                lower_nx, lower_ny = altaz_to_normalized_xy(
                    float(lower_alt),
                    az,
                    view_center,
                    edge_fov_deg=float(edge_fov_deg),
                )
            except Exception:
                continue
            lower_points.append(QPointF(*normalized_to_screen_xy(lower_nx, lower_ny, geometry)))
            lower_altitudes.append(float(lower_alt))
            lower_azimuths.append(float(seam_az_deg_value))

        if len(lower_points) < 2 or len(lower_altitudes) < 2:
            point_index += len(fragment)
            continue

        fragment_alpha = max(
            float(RIDGE_GLOW_SKY_ALPHA_FLOOR),
            min(1.0, float(opacity_scale) * float(RIDGE_GLOW_SKY_ALPHA_SCALE)),
        )
        painter.setPen(Qt.PenStyle.NoPen)
        layer_count = len(layer_specs)
        boundary_points = [lower_points]
        for step_index in range(layer_count):
            layer_spec = layer_specs[step_index]
            upper_altitudes = _ridge_glow_directional_altitudes(
                [
                    float(lower_alt) + float(layer_spec.upper_alt_offset_deg)
                    for lower_alt in lower_altitudes
                ],
                int(layer_spec.window_size),
                center_weight_scale=float(layer_spec.focus_scale),
            )
            upper_boundary = [
                QPointF(
                    *normalized_to_screen_xy(
                        *altaz_to_normalized_xy(
                            float(smoothed_alt),
                            (float(seam_az_deg_value) + seam_az_deg) % 360.0,
                            view_center,
                            edge_fov_deg=float(edge_fov_deg),
                        ),
                        geometry,
                    )
                )
                for seam_az_deg_value, smoothed_alt in zip(lower_azimuths, upper_altitudes)
            ]
            if len(upper_boundary) < 2:
                continue
            lower_boundary = boundary_points[-1]
            band_polygon = QPolygonF(lower_boundary + list(reversed(upper_boundary)))
            if band_polygon.isEmpty():
                continue
            if step_index == layer_count - 1:
                painter.setBrush(
                    _ridge_glow_hatch_brush(
                        glow_rgb,
                        fragment_alpha * float(layer_spec.alpha_scale),
                    )
                )
            else:
                color = QColor(*glow_rgb)
                color.setAlphaF(
                    max(
                        0.0,
                        min(1.0, fragment_alpha * float(layer_spec.alpha_scale)),
                    )
                )
                painter.setBrush(QBrush(color))
            painter.drawPolygon(band_polygon)
            boundary_points.append(upper_boundary)

        point_index += len(fragment)


def estimate_ridge_glow_color_for_sky_image(
    sky_img: QImage,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    *,
    edge_fov_deg: float,
    alt_deg: float = 0.0,
) -> tuple[int, int, int] | None:
    """Return a representative glow color from the sky image and fixed warm tone."""
    if sky_img.isNull() or geometry.radius < 1:
        return None

    width = sky_img.width()
    height = sky_img.height()
    if width <= 0 or height <= 0:
        return None

    cx = float(geometry.center[0])
    cy = float(geometry.center[1])
    radius = float(geometry.radius)
    samples: list[tuple[float, tuple[int, int, int]]] = []
    az_deg = 0.0
    while az_deg <= 360.0 + 1.0e-6:
        nx, ny = altaz_to_normalized_xy(
            float(alt_deg),
            float(az_deg),
            view_center,
            edge_fov_deg=edge_fov_deg,
        )
        x = int(round(cx + (nx * radius)))
        y = int(round(cy + (ny * radius)))
        if 0 <= x < width and 0 <= y < height:
            color = sky_img.pixelColor(x, y)
            if color.alpha() > 0:
                samples.append(
                    (
                        dimalt_ring_brightness_score_from_rgba(
                            color.red(),
                            color.green(),
                            color.blue(),
                            color.alpha(),
                        ),
                        (color.red(), color.green(), color.blue()),
                    )
                )
        az_deg += 30.0

    if not samples:
        return None

    samples.sort(key=lambda item: item[0])
    mid = len(samples) // 2
    if len(samples) % 2 == 1:
        _, rgb = samples[mid]
    else:
        _, rgb = samples[mid - 1]

    red, green, blue = rgb
    sky_mix = 0.3
    warm_mix = 0.7
    base_red, base_green, base_blue = NIGHT_LIGHTS_GLOW_RGB
    return (
        int(round((float(red) * sky_mix) + (float(base_red) * warm_mix))),
        int(round((float(green) * sky_mix) + (float(base_green) * warm_mix))),
        int(round((float(blue) * sky_mix) + (float(base_blue) * warm_mix))),
    )


def draw_ridge_glow_normal(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    profile: NightLightGlowProfile | None,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
    viewer_data: ViewerData | None = None,
    view_center: tuple[float, float] | None = None,
    opacity: float = 1.0,
    ridge_glow_opacity: float | None = None,
    ridge_glow_color_rgb: tuple[int, int, int] | None = None,
    sun_alt_deg: float | None = None,
    edge_fov_deg: float | None = None,
    content_fov_deg: float | None = None,
) -> None:
    if viewer_data is not None:
        view_center = tuple(float(value) for value in viewer_data.view_center)
        edge_fov_deg = float(viewer_data.edge_fov_deg)
        content_fov_deg = float(viewer_data.content_fov_deg)
    if view_center is None:
        raise ValueError("view_center is required when viewer_data is not provided")
    if edge_fov_deg is None:
        raise ValueError("edge_fov_deg is required when viewer_data is not provided")
    if content_fov_deg is None:
        raise ValueError("content_fov_deg is required when viewer_data is not provided")

    layer_opacity = max(0.0, min(1.0, float(opacity)))
    ridge_glow_layer_opacity = (
        layer_opacity if ridge_glow_opacity is None else max(0.0, min(1.0, float(ridge_glow_opacity)))
    )
    sun_factor = 1.0 if sun_alt_deg is None else night_light_strength_factor(sun_alt_deg)
    if (
        profile is None
        or not profile.samples
        or (layer_opacity <= 0.0 and ridge_glow_layer_opacity <= 0.0)
        or sun_factor <= 0.0
    ):
        return
    if terrain_profile_altaz is None and not terrain_secondary_ridges_altaz_layers:
        return

    samples = [sample for sample in profile.samples if float(sample.strength) > 0.0]
    if not samples:
        return

    seam_az_deg = (float(view_center[1]) + 180.0) % 360.0
    fill_rgb = tuple(int(value) for value in RIDGE_GLOW_FILL_RGB)
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)

    layer_altaz_sets: list[list[tuple[float, float]]] = []
    layer_profile_samples: list[tuple] = []
    if terrain_secondary_ridges_altaz_layers:
        band_profile_samples = tuple(band_profile.samples for band_profile in getattr(profile, "band_profiles", ()))
        layer_altaz_sets.extend(list(layer) for layer in terrain_secondary_ridges_altaz_layers)
        for band_index, _layer in enumerate(terrain_secondary_ridges_altaz_layers):
            if band_profile_samples:
                layer_profile_samples.append(
                    band_profile_samples[min(band_index, len(band_profile_samples) - 1)]
                )
            else:
                layer_profile_samples.append(profile.samples)

    if not layer_altaz_sets and not terrain_profile_altaz:
        painter.restore()
        return

    if layer_altaz_sets:
        for layer_index, layer_altaz in enumerate(layer_altaz_sets):
            if not layer_altaz or len(layer_altaz) < 2:
                continue
            selected_samples = layer_profile_samples[min(layer_index, len(layer_profile_samples) - 1)]
            ordered = sorted(
                (sample for sample in selected_samples if float(sample.strength) > 0.0),
                key=lambda sample: _seam_relative_azimuth_deg(sample.azimuth_deg, seam_az_deg),
            )
            if len(ordered) < 2:
                continue

            night_az = np.asarray(
                [_seam_relative_azimuth_deg(sample.azimuth_deg, seam_az_deg) for sample in ordered],
                dtype=np.float64,
            )
            night_strengths = np.asarray([float(sample.strength) for sample in ordered], dtype=np.float64)
            if night_az.size < 2:
                continue

            night_az_ext = np.concatenate([night_az[-1:] - 360.0, night_az, night_az[:1] + 360.0])
            night_strengths_ext = np.concatenate([night_strengths[-1:], night_strengths, night_strengths[:1]])
            layer_samples = sorted(
                (
                    (float(alt_deg), float(az_deg) % 360.0)
                    for alt_deg, az_deg in layer_altaz
                ),
                key=lambda item: _seam_relative_azimuth_deg(item[1], seam_az_deg),
            )
            layer_az = np.asarray(
                [_seam_relative_azimuth_deg(az, seam_az_deg) for _, az in layer_samples],
                dtype=np.float64,
            )
            layer_horizon_alt = np.asarray([alt for alt, _ in layer_samples], dtype=np.float64)
            if layer_az.size < 2:
                continue

            previous_layer_altaz = (
                [(0.0, az_deg) for _, az_deg in layer_samples]
                if layer_index == 0
                else layer_altaz_sets[layer_index - 1]
            )
            band_lower_edge_alts = _band_lower_edge_altitudes(layer_horizon_alt, layer_az, previous_layer_altaz)
            layer_strengths = np.interp(layer_az, night_az_ext, night_strengths_ext)
            draw_az = layer_az
            projected_points: list[tuple[float, float]] = []
            projected_draw_az: list[float] = []
            projected_current_alts: list[float] = []
            projected_lower_alts: list[float] = []
            strengths: list[float] = []
            for seam_az_deg_value, current_alt, lower_alt, strength in zip(
                draw_az.tolist(),
                layer_horizon_alt.tolist(),
                band_lower_edge_alts.tolist(),
                layer_strengths.tolist(),
            ):
                try:
                    nx, ny = altaz_to_normalized_xy(
                        float(lower_alt),
                        (float(seam_az_deg_value) + seam_az_deg) % 360.0,
                        view_center,
                        edge_fov_deg=float(edge_fov_deg),
                    )
                except Exception:
                    continue
                projected_points.append((float(nx), float(ny)))
                projected_draw_az.append(float(seam_az_deg_value))
                projected_current_alts.append(float(current_alt))
                projected_lower_alts.append(float(lower_alt))
                strengths.append(float(strength))

            if len(projected_points) < 2:
                continue

            layer_distance_km = _layer_distance_km(profile, layer_index)
            distance_scale = max(0.35, min(1.0, float(np.sqrt(0.5 / max(0.5, layer_distance_km)))))
            width_scale = max(0.35, min(1.0, float(np.sqrt(0.5 / max(0.5, layer_distance_km)))))
            band_thickness_deg = float(profile.band_half_width_deg) * 2.0 * width_scale

            street_alpha = min(
                1.0,
                1.0 * layer_opacity * sun_factor * distance_scale,
            )
            if street_alpha <= 0.0:
                continue
            color = QColor(*fill_rgb)
            color.setAlphaF(street_alpha)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(color)

            point_index = 0
            for fragment in split_by_gaps(projected_points):
                if len(fragment) < 2:
                    point_index += len(fragment)
                    continue
                frag_strengths = strengths[point_index : point_index + len(fragment)]
                frag_strength = max(float(value) for value in frag_strengths)
                if frag_strength < RIDGE_GLOW_MIN_BRIGHTNESS:
                    point_index += len(fragment)
                    continue

                lower_points: list[QPointF] = []
                upper_points: list[QPointF] = []
                for seam_az_deg_value, current_alt, lower_alt, _strength in zip(
                    projected_draw_az[point_index : point_index + len(fragment)],
                    projected_current_alts[point_index : point_index + len(fragment)],
                    projected_lower_alts[point_index : point_index + len(fragment)],
                    frag_strengths,
                ):
                    az = (float(seam_az_deg_value) + seam_az_deg) % 360.0
                    upper_alt = float(current_alt) + float(band_thickness_deg)
                    if upper_alt <= float(lower_alt):
                        lower_points = []
                        upper_points = []
                        break
                    try:
                        lower_nx, lower_ny = altaz_to_normalized_xy(
                            lower_alt,
                            az,
                            view_center,
                            edge_fov_deg=float(edge_fov_deg),
                        )
                        upper_nx, upper_ny = altaz_to_normalized_xy(
                            upper_alt,
                            az,
                            view_center,
                            edge_fov_deg=float(edge_fov_deg),
                        )
                    except Exception:
                        continue
                    lower_points.append(QPointF(*normalized_to_screen_xy(lower_nx, lower_ny, geometry)))
                    upper_points.append(QPointF(*normalized_to_screen_xy(upper_nx, upper_ny, geometry)))

                if len(lower_points) < 2 or len(upper_points) < 2:
                    point_index += len(fragment)
                    continue

                band_polygon = QPolygonF(lower_points + list(reversed(upper_points)))
                if band_polygon.isEmpty():
                    point_index += len(fragment)
                    continue
                band_polygon.append(lower_points[0])
                painter.drawPolygon(band_polygon)
                point_index += len(fragment)

    if terrain_profile_altaz:
        main_az_strengths = sorted(
            (sample for sample in profile.samples if float(sample.strength) > 0.0),
            key=lambda sample: _seam_relative_azimuth_deg(sample.azimuth_deg, seam_az_deg),
        )
        if len(main_az_strengths) >= 2:
            main_az = np.asarray(
                [_seam_relative_azimuth_deg(sample.azimuth_deg, seam_az_deg) for sample in main_az_strengths],
                dtype=np.float64,
            )
            main_strengths = np.asarray([float(sample.strength) for sample in main_az_strengths], dtype=np.float64)
            main_az_ext = np.concatenate([main_az[-1:] - 360.0, main_az, main_az[:1] + 360.0])
            main_strengths_ext = np.concatenate([main_strengths[-1:], main_strengths, main_strengths[:1]])
            main_samples = sorted(
                (
                    (float(alt_deg), float(az_deg) % 360.0)
                    for alt_deg, az_deg in terrain_profile_altaz
                ),
                key=lambda item: _seam_relative_azimuth_deg(item[1], seam_az_deg),
            )
            main_azimuths = np.asarray(
                [_seam_relative_azimuth_deg(az, seam_az_deg) for _, az in main_samples],
                dtype=np.float64,
            )
            main_alts = np.asarray([alt for alt, _ in main_samples], dtype=np.float64)
            if main_azimuths.size >= 2:
                main_strength_interp = np.interp(main_azimuths, main_az_ext, main_strengths_ext)
                projected_points: list[tuple[float, float]] = []
                projected_draw_az: list[float] = []
                projected_alts: list[float] = []
                strengths: list[float] = []
                for seam_az_deg_value, alt_value, strength_value in zip(
                    main_azimuths.tolist(),
                    main_alts.tolist(),
                    main_strength_interp.tolist(),
                ):
                    try:
                        nx, ny = altaz_to_normalized_xy(
                            float(alt_value),
                            (float(seam_az_deg_value) + seam_az_deg) % 360.0,
                            view_center,
                            edge_fov_deg=float(edge_fov_deg),
                        )
                    except Exception:
                        continue
                    projected_points.append((float(nx), float(ny)))
                    projected_draw_az.append(float(seam_az_deg_value))
                    projected_alts.append(float(alt_value))
                    strengths.append(float(strength_value))

                if len(projected_points) >= 2:
                    _draw_ridge_glow_fragments(
                        painter,
                        geometry=geometry,
                        view_center=view_center,
                        edge_fov_deg=edge_fov_deg,
                        seam_az_deg=seam_az_deg,
                        projected_points=projected_points,
                        projected_draw_az=projected_draw_az,
                        projected_alts=projected_alts,
                        strengths=strengths,
                        color_rgb=(
                            RIDGE_GLOW_SKY_RGB
                            if ridge_glow_color_rgb is None
                            else ridge_glow_color_rgb
                        ),
                        opacity_scale=ridge_glow_layer_opacity * sun_factor,
                        band_thickness_deg=float(profile.band_half_width_deg) * 1.4,
                    )

    painter.restore()
