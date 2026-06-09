from __future__ import annotations

import math

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPolygonF

from ..astro import altaz_to_normalized_xy
from ..night_lights import NightLightGlowProfile, night_light_strength_factor
from ..types import ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy
from .guides import split_by_gaps

NIGHT_LIGHTS_BASE_ALPHA_SCALE = 1.0
NIGHT_LIGHTS_BAND_SPECS: tuple[tuple[float, float], ...] = (
    (0.1, 3.2),
    (0.3, 1.7),
    (0.50, 0.9),
    (1.0, 0.5),
)
NIGHT_LIGHTS_MIN_BRIGHTNESS = 0.02
NIGHT_LIGHTS_GLOW_RGB = (244, 246, 248)
NIGHT_LIGHTS_DRAW_AZIMUTH_STEP_DEG = 0.5
NIGHT_LIGHTS_DISTANCE_NEAR_KM = 0.5
NIGHT_LIGHTS_DISTANCE_FAR_KM = 128.0


def _band_width_px(
    *,
    center_alt_deg: float,
    band_half_width_deg: float,
    view_center: tuple[float, float],
    geometry: ScreenGeometry,
    edge_fov_deg: float,
) -> float:
    try:
        nx1, ny1 = altaz_to_normalized_xy(
            float(center_alt_deg),
            float(view_center[1]),
            view_center,
            edge_fov_deg=float(edge_fov_deg),
        )
        nx2, ny2 = altaz_to_normalized_xy(
            float(center_alt_deg) + float(band_half_width_deg),
            float(view_center[1]),
            view_center,
            edge_fov_deg=float(edge_fov_deg),
        )
    except Exception:
        return max(1.0, float(geometry.radius) * 0.02)
    x1, y1 = normalized_to_screen_xy(nx1, ny1, geometry)
    x2, y2 = normalized_to_screen_xy(nx2, ny2, geometry)
    return max(1.0, math.hypot(float(x2) - float(x1), float(y2) - float(y1)) * 2.0)


def _build_band_polygon(
    points: list[QPointF],
    *,
    half_width_px: float,
) -> QPolygonF:
    if len(points) < 2:
        return QPolygonF()

    half_width = max(0.5, float(half_width_px))
    left_points: list[QPointF] = []
    right_points: list[QPointF] = []
    last_index = len(points) - 1
    for index, point in enumerate(points):
        if index == 0:
            tangent_x = float(points[1].x()) - float(point.x())
            tangent_y = float(points[1].y()) - float(point.y())
        elif index == last_index:
            tangent_x = float(point.x()) - float(points[index - 1].x())
            tangent_y = float(point.y()) - float(points[index - 1].y())
        else:
            tangent_x = float(points[index + 1].x()) - float(points[index - 1].x())
            tangent_y = float(points[index + 1].y()) - float(points[index - 1].y())
        length = math.hypot(tangent_x, tangent_y)
        if length <= 1.0e-9:
            normal_x, normal_y = 0.0, 1.0
        else:
            normal_x = -tangent_y / length
            normal_y = tangent_x / length
        offset_x = normal_x * half_width
        offset_y = normal_y * half_width
        left_points.append(QPointF(float(point.x()) + offset_x, float(point.y()) + offset_y))
        right_points.append(QPointF(float(point.x()) - offset_x, float(point.y()) - offset_y))

    polygon = QPolygonF(left_points + list(reversed(right_points)))
    if not polygon.isEmpty():
        polygon.append(left_points[0])
    return polygon


def _seam_relative_azimuth_deg(azimuth_deg: float, seam_az_deg: float) -> float:
    return (float(azimuth_deg) - float(seam_az_deg)) % 360.0


def _draw_night_light_glow_impl(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    profile: NightLightGlowProfile | None,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
    view_center: tuple[float, float],
    opacity: float = 1.0,
    sun_alt_deg: float | None = None,
    edge_fov_deg: float,
) -> None:
    layer_opacity = max(0.0, min(1.0, float(opacity)))
    sun_factor = 1.0 if sun_alt_deg is None else night_light_strength_factor(sun_alt_deg)
    if profile is None or not profile.samples or layer_opacity <= 0.0 or sun_factor <= 0.0:
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
    flat_layer_altaz = [
        (0.0, float(sample.azimuth_deg) % 360.0)
        for sample in profile.samples
    ]
    if terrain_secondary_ridges_altaz_layers:
        if terrain_profile_altaz:
            layer_altaz_sets.append(list(terrain_profile_altaz))
        else:
            layer_altaz_sets.append(flat_layer_altaz)
        layer_profile_samples.append(profile.samples)
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
    elif terrain_profile_altaz:
        layer_altaz_sets.append(list(terrain_profile_altaz))
        layer_profile_samples.append(profile.samples)
    else:
        layer_altaz_sets.append(flat_layer_altaz)
        layer_profile_samples.append(profile.samples)

    if not layer_altaz_sets:
        painter.restore()
        return

    band_half_width_deg = float(profile.band_half_width_deg)
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
        draw_az = layer_az
        projected_points: list[tuple[float, float]] = []
        strengths: list[float] = []
        center_alts: list[float] = []
        for seam_az_deg_value, horizon_alt, strength in zip(draw_az.tolist(), layer_horizon_alt.tolist(), layer_strengths.tolist()):
            center_alt = float(horizon_alt) + float(profile.band_center_offset_deg)
            try:
                nx, ny = altaz_to_normalized_xy(
                    center_alt,
                    (float(seam_az_deg_value) + seam_az_deg) % 360.0,
                    view_center,
                    edge_fov_deg=float(edge_fov_deg),
                )
            except Exception:
                continue
            projected_points.append((float(nx), float(ny)))
            strengths.append(float(strength))
            center_alts.append(center_alt)

        if len(projected_points) < 2:
            continue

        width_px = _band_width_px(
            center_alt_deg=float(sum(center_alts) / len(center_alts)),
            band_half_width_deg=band_half_width_deg,
            view_center=view_center,
            geometry=geometry,
            edge_fov_deg=edge_fov_deg,
        )
        if width_px <= 0.0:
            continue

        point_index = 0
        for fragment in split_by_gaps(projected_points):
            if len(fragment) < 2:
                point_index += len(fragment)
                continue
            frag_strengths = strengths[point_index:point_index + len(fragment)]
            frag_points = [QPointF(*normalized_to_screen_xy(nx, ny, geometry)) for nx, ny in fragment]
            frag_strength = max(float(value) for value in frag_strengths)
            if frag_strength < NIGHT_LIGHTS_MIN_BRIGHTNESS:
                point_index += len(fragment)
                continue

            alpha = min(
                1.0,
                frag_strength
                * NIGHT_LIGHTS_BASE_ALPHA_SCALE
                * layer_opacity
                * sun_factor
            )
            color = QColor(*fill_rgb)

            for alpha_scale, width_scale in NIGHT_LIGHTS_BAND_SPECS:
                color.setAlphaF(max(0.0, min(1.0, alpha * alpha_scale)))
                painter.setPen(Qt.PenStyle.NoPen)
                painter.setBrush(color)
                band_polygon = _build_band_polygon(
                    frag_points,
                    half_width_px=width_px * float(width_scale) * 0.5,
                )
                if band_polygon.isEmpty():
                    continue
                painter.drawPolygon(band_polygon)
            point_index += len(fragment)

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
        sun_alt_deg=sun_alt_deg,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
    )
