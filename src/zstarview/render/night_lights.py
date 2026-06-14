from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPolygonF

from ..astro import altaz_to_normalized_xy
from ..night_lights import NightLightGlowProfile, night_light_strength_factor
from ..types import ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy
from .guides import split_by_gaps

NIGHT_LIGHTS_MIN_BRIGHTNESS = 0.02


NIGHT_LIGHTS_GLOW_RGB = (244, 246, 248)
NIGHT_LIGHTS_STREET_LIGHT_GLOW_ALPHA_BASE = 1.0
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


def _seam_relative_azimuth_deg(azimuth_deg: float, seam_az_deg: float) -> float:
    return (float(azimuth_deg) - float(seam_az_deg)) % 360.0


def _draw_night_light_glow_impl(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    profile: NightLightGlowProfile | None,
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
    view_center: tuple[float, float],
    opacity: float = 1.0,
    sun_alt_deg: float | None = None,
    edge_fov_deg: float,
) -> None:
    layer_opacity = max(0.0, min(1.0, float(opacity)))
    sun_factor = 1.0 if sun_alt_deg is None else night_light_strength_factor(sun_alt_deg)
    if (
        profile is None
        or not profile.samples
        or layer_opacity <= 0.0
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

    layer_altaz_sets = [list(layer) for layer in terrain_secondary_ridges_altaz_layers or ()]
    if not layer_altaz_sets:
        painter.restore()
        return

    band_profile_samples = tuple(
        band_profile.samples
        for band_profile in getattr(profile, "band_profiles", ())
    )
    if not band_profile_samples:
        painter.restore()
        return

    # When the terrain pipeline exposes the extra 250m ridge layer, the night-light
    # bands use adjacent ridge pairs: the first street-light band sits between the
    # 250m and 500m ridges, and later bands continue from there.
    use_shifted_ridge_pairs = len(layer_altaz_sets) == len(band_profile_samples) + 1
    ridge_start_index = 1 if use_shifted_ridge_pairs else 0
    band_count = min(len(band_profile_samples), len(layer_altaz_sets) - ridge_start_index)
    if band_count <= 0:
        painter.restore()
        return

    for band_index in range(band_count):
        layer_index = ridge_start_index + band_index
        layer_altaz = layer_altaz_sets[layer_index]
        if not layer_altaz or len(layer_altaz) < 2:
            continue
        selected_samples = tuple(
            sample
            for sample in band_profile_samples[band_index]
            if float(sample.strength) > 0.0
        )
        if len(selected_samples) < 2:
            continue

        ordered = sorted(
            selected_samples,
            key=lambda sample: _seam_relative_azimuth_deg(sample.azimuth_deg, seam_az_deg),
        )
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
            layer_altaz_sets[layer_index - 1]
            if layer_index > 0
            else [(0.0, az_deg) for _, az_deg in layer_samples]
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

        layer_distance_km = _layer_distance_km(profile, band_index)
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
            middle_points: list[QPointF] = []
            upper_points: list[QPointF] = []
            for seam_az_deg_value, current_alt, lower_alt, _strength in zip(
                projected_draw_az[point_index:point_index + len(fragment)],
                projected_current_alts[point_index:point_index + len(fragment)],
                projected_lower_alts[point_index:point_index + len(fragment)],
                frag_strengths,
            ):
                az = (float(seam_az_deg_value) + seam_az_deg) % 360.0
                upper_alt = float(current_alt) + float(band_thickness_deg)
                if upper_alt <= float(lower_alt):
                    lower_points = []
                    break
                mid_alt = float(lower_alt) + ((upper_alt - float(lower_alt)) * 0.5)
                try:
                    lower_nx, lower_ny = altaz_to_normalized_xy(
                        lower_alt,
                        az,
                        view_center,
                        edge_fov_deg=float(edge_fov_deg),
                    )
                    mid_nx, mid_ny = altaz_to_normalized_xy(
                        mid_alt,
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
                middle_points.append(QPointF(*normalized_to_screen_xy(mid_nx, mid_ny, geometry)))
                upper_points.append(QPointF(*normalized_to_screen_xy(upper_nx, upper_ny, geometry)))

            if (
                len(lower_points) < 2
                or len(middle_points) < 2
                or len(upper_points) < 2
            ):
                point_index += len(fragment)
                continue

            lower_half_polygon = QPolygonF(lower_points + list(reversed(middle_points)))
            upper_half_polygon = QPolygonF(middle_points + list(reversed(upper_points)))
            if lower_half_polygon.isEmpty() or upper_half_polygon.isEmpty():
                point_index += len(fragment)
                continue
            painter.setBrush(
                QColor(fill_rgb[0], fill_rgb[1], fill_rgb[2], int(round(max(0.0, min(1.0, street_alpha)) * 255.0)))
            )
            painter.drawPolygon(lower_half_polygon)
            painter.setBrush(
                QColor(
                    fill_rgb[0],
                    fill_rgb[1],
                    fill_rgb[2],
                    int(round(max(0.0, min(1.0, street_alpha * 0.5)) * 255.0)),
                )
            )
            painter.drawPolygon(upper_half_polygon)
            point_index += len(fragment)

    painter.restore()


def draw_night_light_glow_normal(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    profile: NightLightGlowProfile | None,
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
    viewer_data: ViewerData | None = None,
    view_center: tuple[float, float] | None = None,
    opacity: float = 1.0,
    sun_alt_deg: float | None = None,
    edge_fov_deg: float | None = None,
    fast_mode: bool = False,
) -> None:
    """Draw the street-light glow."""
    if fast_mode:
        return
    if viewer_data is not None:
        view_center = tuple(float(value) for value in viewer_data.view_center)
        edge_fov_deg = float(viewer_data.edge_fov_deg)
    if view_center is None:
        raise ValueError("view_center is required when viewer_data is not provided")
    if edge_fov_deg is None:
        raise ValueError("edge_fov_deg is required when viewer_data is not provided")
    _draw_night_light_glow_impl(
        painter,
        geometry=geometry,
        profile=profile,
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
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
    viewer_data: ViewerData | None = None,
    view_center: tuple[float, float] | None = None,
    opacity: float = 1.0,
    sun_alt_deg: float | None = None,
    edge_fov_deg: float | None = None,
    fast_mode: bool = False,
) -> None:
    """Compatibility wrapper kept for existing callers."""
    if fast_mode:
        return
    if viewer_data is not None:
        view_center = tuple(float(value) for value in viewer_data.view_center)
        edge_fov_deg = float(viewer_data.edge_fov_deg)
    if view_center is None:
        raise ValueError("view_center is required when viewer_data is not provided")
    if edge_fov_deg is None:
        raise ValueError("edge_fov_deg is required when viewer_data is not provided")
    draw_night_light_glow_normal(
        painter,
        geometry=geometry,
        profile=profile,
        terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
        viewer_data=viewer_data,
        view_center=view_center,
        opacity=opacity,
        sun_alt_deg=sun_alt_deg,
        edge_fov_deg=edge_fov_deg,
    )
