from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QBrush, QColor, QPainter, QPolygonF

from ..astro import altaz_to_normalized_xy
from ..night_lights import NightLightGlowProfile, night_light_strength_factor
from ..types import ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy
from .guides import split_by_gaps

NIGHT_LIGHTS_MIN_BRIGHTNESS = 0.02
NIGHT_LIGHTS_GLOW_RGB = (244, 246, 248)
NIGHT_LIGHTS_SKY_GLOW_RGB = NIGHT_LIGHTS_GLOW_RGB
NIGHT_LIGHTS_SKY_GLOW_ALPHA_FLOOR = 0.01
NIGHT_LIGHTS_STREET_LIGHT_GLOW_ALPHA_BASE = 1.0
NIGHT_LIGHTS_SKY_GLOW_ALPHA_SCALE = 0.8
NIGHT_LIGHTS_SKY_GLOW_SEGMENT_COUNT = 3
NIGHT_LIGHTS_SKY_GLOW_WINDOW_EXPONENT = 1.0
NIGHT_LIGHTS_SKY_GLOW_ALPHA_RATIO = 0.42
NIGHT_LIGHTS_SKY_GLOW_WIDTH_WEIGHTS = (10.0, 20.0, 40.0)
NIGHT_LIGHTS_SKY_GLOW_WINDOW_SIZES = tuple(
    2 ** (index + 3) for index in range(NIGHT_LIGHTS_SKY_GLOW_SEGMENT_COUNT)
)
NIGHT_LIGHTS_DISTANCE_NEAR_KM = 0.5
NIGHT_LIGHTS_DISTANCE_FAR_KM = 128.0
NIGHT_LIGHTS_MIN_WIDTH_SCALE = 0.35


def _night_light_distance_scale(distance_km: float) -> float:
    d = max(
        float(NIGHT_LIGHTS_DISTANCE_NEAR_KM),
        min(float(NIGHT_LIGHTS_DISTANCE_FAR_KM), float(distance_km)),
    )
    return float(np.sqrt(float(NIGHT_LIGHTS_DISTANCE_NEAR_KM) / d))


def _night_light_band_width_scale(distance_km: float) -> float:
    scale = _night_light_distance_scale(distance_km)
    return max(float(NIGHT_LIGHTS_MIN_WIDTH_SCALE), min(1.0, scale))


def _layer_distance_km(
    profile: NightLightGlowProfile,
    layer_index: int,
) -> float:
    band_profiles = tuple(getattr(profile, "band_profiles", ()))
    assert band_profiles, "profile.band_profiles must not be empty"
    assert 0 <= int(layer_index) < len(band_profiles), "layer_index out of range"
    return float(band_profiles[int(layer_index)].max_distance_km)


def _band_lower_edge_altitudes(
    current_layer_altitudes: np.ndarray,
    previous_layer_altitudes: np.ndarray | None,
) -> np.ndarray:
    current = np.asarray(current_layer_altitudes, dtype=np.float64)
    if previous_layer_altitudes is None:
        previous = np.zeros_like(current)
    else:
        previous = np.asarray(previous_layer_altitudes, dtype=np.float64)
        if previous.shape != current.shape:
            raise ValueError("previous_layer_altitudes must match current_layer_altitudes")
    return 0.9 * previous + 0.1 * current


def _seam_relative_azimuth_deg(azimuth_deg: float, seam_az_deg: float) -> float:
    return (float(azimuth_deg) - float(seam_az_deg)) % 360.0


def _sky_glow_step_boundaries() -> np.ndarray:
    step_count = int(NIGHT_LIGHTS_SKY_GLOW_SEGMENT_COUNT)
    if step_count < 1:
        raise ValueError("NIGHT_LIGHTS_SKY_GLOW_SEGMENT_COUNT must be positive")
    widths = np.asarray(NIGHT_LIGHTS_SKY_GLOW_WIDTH_WEIGHTS, dtype=np.float64)
    if widths.size != step_count:
        raise ValueError("NIGHT_LIGHTS_SKY_GLOW_WIDTH_WEIGHTS must match segment count")
    widths = widths / np.sum(widths)
    boundaries = np.concatenate(([0.0], np.cumsum(widths)))
    boundaries[-1] = 1.0
    return boundaries


def _sky_glow_step_alpha_scales() -> np.ndarray:
    step_count = int(NIGHT_LIGHTS_SKY_GLOW_SEGMENT_COUNT)
    if step_count < 1:
        raise ValueError("NIGHT_LIGHTS_SKY_GLOW_SEGMENT_COUNT must be positive")
    ratio = max(0.0, min(1.0, float(NIGHT_LIGHTS_SKY_GLOW_ALPHA_RATIO)))
    return np.asarray(
        [ratio**index for index in range(step_count)],
        dtype=np.float64,
    )


def _sky_glow_directional_altitudes(raw_altitudes: list[float], window_size: int) -> list[float]:
    if len(raw_altitudes) < 2:
        return list(raw_altitudes)
    window_size = max(1, int(window_size))
    if window_size == 1:
        return list(raw_altitudes)

    raw_altitudes_array = np.asarray([float(value) for value in raw_altitudes], dtype=np.float64)
    max_down_step_map = {
        4: 0.4,
        8: 0.2,
        16: 0.1,
        32: 0.05,
    }
    max_down_step_deg = float(max_down_step_map.get(window_size, 0.05 * (32.0 / float(window_size))))
    max_down_step_deg = max(0.0, max_down_step_deg)

    forward = np.asarray(raw_altitudes_array, dtype=np.float64).copy()
    for index in range(1, forward.size):
        current = float(raw_altitudes_array[index])
        previous = float(forward[index - 1])
        if current >= previous:
            forward[index] = current
        else:
            forward[index] = max(current, previous - max_down_step_deg)

    backward = np.asarray(raw_altitudes_array, dtype=np.float64).copy()
    for index in range(backward.size - 2, -1, -1):
        current = float(raw_altitudes_array[index])
        next_value = float(backward[index + 1])
        if current >= next_value:
            backward[index] = current
        else:
            backward[index] = max(current, next_value - max_down_step_deg)

    combined = np.maximum(forward, backward)
    return [float(value) for value in combined.tolist()]


def _draw_sky_glow_fragments(
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
    step_alpha_scales = _sky_glow_step_alpha_scales()
    point_index = 0
    for fragment in split_by_gaps(projected_points):
        if len(fragment) < 2:
            point_index += len(fragment)
            continue

        frag_strengths = strengths[point_index:point_index + len(fragment)]
        frag_strength = max(float(value) for value in frag_strengths)
        if frag_strength < NIGHT_LIGHTS_MIN_BRIGHTNESS:
            point_index += len(fragment)
            continue

        lower_points: list[QPointF] = []
        upper_raw_altitudes: list[float] = []
        frag_altaz = zip(
            projected_draw_az[point_index:point_index + len(fragment)],
            projected_alts[point_index:point_index + len(fragment)],
        )
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
            upper_raw_altitudes.append(float(lower_alt) + float(band_thickness_deg))

        if len(lower_points) < 2 or len(upper_raw_altitudes) < 2:
            point_index += len(fragment)
            continue

        upper_altitude_bands = [
            _sky_glow_directional_altitudes(upper_raw_altitudes, window_size)
            for window_size in NIGHT_LIGHTS_SKY_GLOW_WINDOW_SIZES
        ]
        upper_boundaries = [
            [
                QPointF(
                    *normalized_to_screen_xy(*altaz_to_normalized_xy(
                        float(smoothed_alt),
                        (float(seam_az_deg_value) + seam_az_deg) % 360.0,
                        view_center,
                        edge_fov_deg=float(edge_fov_deg),
                    ), geometry)
                )
                for seam_az_deg_value, smoothed_alt in zip(
                    projected_draw_az[point_index:point_index + len(fragment)],
                    smoothed_altitudes,
                )
            ]
            for smoothed_altitudes in upper_altitude_bands
        ]
        boundary_points = [lower_points, *upper_boundaries]
        if any(len(points) < 2 for points in boundary_points):
            point_index += len(fragment)
            continue

        fragment_alpha = max(
            float(NIGHT_LIGHTS_SKY_GLOW_ALPHA_FLOOR),
            min(1.0, float(opacity_scale) * float(NIGHT_LIGHTS_SKY_GLOW_ALPHA_SCALE)),
        )
        painter.setPen(Qt.PenStyle.NoPen)
        for step_index in range(len(upper_boundaries) - 1, -1, -1):
            lower_boundary = boundary_points[step_index]
            upper_boundary = boundary_points[step_index + 1]
            band_polygon = QPolygonF(lower_boundary + list(reversed(upper_boundary)))
            if band_polygon.isEmpty():
                continue
            color = QColor(*color_rgb)
            color.setAlphaF(
                max(
                    0.0,
                    min(
                        1.0,
                        fragment_alpha * float(step_alpha_scales[step_index]),
                    ),
                )
            )
            painter.setBrush(QBrush(color))
            painter.drawPolygon(band_polygon)

        point_index += len(fragment)


def _draw_night_light_glow_impl(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    profile: NightLightGlowProfile | None,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
    view_center: tuple[float, float],
    opacity: float = 1.0,
    sky_glow_opacity: float | None = None,
    sun_alt_deg: float | None = None,
    edge_fov_deg: float,
) -> None:
    layer_opacity = max(0.0, min(1.0, float(opacity)))
    sky_glow_layer_opacity = layer_opacity if sky_glow_opacity is None else max(0.0, min(1.0, float(sky_glow_opacity)))
    sun_factor = 1.0 if sun_alt_deg is None else night_light_strength_factor(sun_alt_deg)
    if (
        profile is None
        or not profile.samples
        or (layer_opacity <= 0.0 and sky_glow_layer_opacity <= 0.0)
        or sun_factor <= 0.0
    ):
        return
    samples = [sample for sample in profile.samples if float(sample.strength) > 0.0]
    if not samples:
        return

    seam_az_deg = (float(view_center[1]) + 180.0) % 360.0
    fill_rgb = tuple(int(value) for value in NIGHT_LIGHTS_GLOW_RGB)
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)

    layer_altaz_sets: list[list[tuple[float, float]]] = []
    layer_profile_samples: list[tuple] = []
    if terrain_secondary_ridges_altaz_layers:
        band_profile_samples = tuple(
            band_profile.samples
            for band_profile in getattr(profile, "band_profiles", ())
        )
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
        # Secondary ridge layers are zero-based here: layer 0 is the first secondary ridge.
        for layer_index, layer_altaz in enumerate(layer_altaz_sets):
            if not layer_altaz or len(layer_altaz) < 2:
                continue
            selected_samples = layer_profile_samples[min(layer_index, len(layer_profile_samples) - 1)]
            ordered = sorted(
                (
                    sample
                    for sample in selected_samples
                    if float(sample.strength) > 0.0
                ),
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

            layer_strengths = np.interp(layer_az, night_az_ext, night_strengths_ext)
            band_lower_edge_alts = _band_lower_edge_altitudes(layer_horizon_alt, layer_horizon_alt)
            draw_az = layer_az
            projected_points: list[tuple[float, float]] = []
            projected_draw_az: list[float] = []
            projected_lower_alts: list[float] = []
            strengths: list[float] = []
            for seam_az_deg_value, lower_alt, strength in zip(
                draw_az.tolist(),
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
                projected_lower_alts.append(float(lower_alt))
                strengths.append(float(strength))

            if len(projected_points) < 2:
                continue

            layer_distance_km = _layer_distance_km(profile, layer_index)
            distance_scale = _night_light_distance_scale(layer_distance_km)
            width_scale = _night_light_band_width_scale(layer_distance_km)
            band_thickness_deg = float(profile.band_half_width_deg) * 2.0 * width_scale

            street_alpha = min(
                1.0,
                float(NIGHT_LIGHTS_STREET_LIGHT_GLOW_ALPHA_BASE)
                * layer_opacity
                * sun_factor
                * distance_scale,
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
                frag_strengths = strengths[point_index:point_index + len(fragment)]
                frag_strength = max(float(value) for value in frag_strengths)
                if frag_strength < NIGHT_LIGHTS_MIN_BRIGHTNESS:
                    point_index += len(fragment)
                    continue

                lower_points: list[QPointF] = []
                upper_points: list[QPointF] = []
                for seam_az_deg_value, lower_alt, _strength in zip(
                    projected_draw_az[point_index:point_index + len(fragment)],
                    projected_lower_alts[point_index:point_index + len(fragment)],
                    frag_strengths,
                ):
                    az = (float(seam_az_deg_value) + seam_az_deg) % 360.0
                    upper_alt = float(lower_alt) + float(band_thickness_deg)
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
            (
                sample
                for sample in profile.samples
                if float(sample.strength) > 0.0
            ),
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
                    _draw_sky_glow_fragments(
                        painter,
                        geometry=geometry,
                        view_center=view_center,
                        edge_fov_deg=edge_fov_deg,
                        seam_az_deg=seam_az_deg,
                        projected_points=projected_points,
                        projected_draw_az=projected_draw_az,
                        projected_alts=projected_alts,
                        strengths=strengths,
                        color_rgb=NIGHT_LIGHTS_SKY_GLOW_RGB,
                        opacity_scale=sky_glow_layer_opacity * sun_factor,
                        band_thickness_deg=float(profile.band_half_width_deg) * 1.4,
                    )

    painter.restore()


def draw_night_light_glow_normal(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    profile: NightLightGlowProfile | None,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
    viewer_data: ViewerData | None = None,
    view_center: tuple[float, float] | None = None,
    opacity: float = 1.0,
    sky_glow_opacity: float | None = None,
    sun_alt_deg: float | None = None,
    edge_fov_deg: float | None = None,
    content_fov_deg: float | None = None,
) -> None:
    """Draw the full normal-mode night light glow."""
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
    _draw_night_light_glow_impl(
        painter,
        geometry=geometry,
        profile=profile,
        terrain_profile_altaz=terrain_profile_altaz,
        terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
        view_center=view_center,
        opacity=opacity,
        sky_glow_opacity=sky_glow_opacity,
        sun_alt_deg=sun_alt_deg,
        edge_fov_deg=edge_fov_deg,
    )


def draw_night_light_glow(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    profile: NightLightGlowProfile | None,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
    viewer_data: ViewerData | None = None,
    view_center: tuple[float, float] | None = None,
    opacity: float = 1.0,
    sky_glow_opacity: float | None = None,
    sun_alt_deg: float | None = None,
    edge_fov_deg: float | None = None,
    content_fov_deg: float | None = None,
    fast_mode: bool = False,
) -> None:
    """Compatibility wrapper kept for existing callers."""
    if fast_mode:
        return
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
    draw_night_light_glow_normal(
        painter,
        geometry=geometry,
        profile=profile,
        terrain_profile_altaz=terrain_profile_altaz,
        terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
        viewer_data=viewer_data,
        view_center=view_center,
        opacity=opacity,
        sky_glow_opacity=sky_glow_opacity,
        sun_alt_deg=sun_alt_deg,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
    )
