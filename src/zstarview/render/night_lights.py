from __future__ import annotations

import math

import numpy as np
from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QPainter, QPen

from ..astro import altaz_to_normalized_xy
from ..night_lights import NightLightGlowProfile
from ..paths import ThemeStyle
from ..types import ScreenGeometry
from .geometry import normalized_to_screen_xy
from .guides import split_by_gaps

NIGHT_LIGHTS_CORE_ALPHA_SCALE = 1.0
NIGHT_LIGHTS_MID_ALPHA_SCALE = 0.44
NIGHT_LIGHTS_OUTER_ALPHA_SCALE = 0.16
NIGHT_LIGHTS_MIN_BRIGHTNESS = 0.02
NIGHT_LIGHTS_GLOW_RGB = (240, 200, 140)
NIGHT_LIGHTS_DRAW_AZIMUTH_STEP_DEG = 0.5


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


def draw_night_light_glow(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewport_rect: QRectF,
    profile: NightLightGlowProfile | None,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    view_center: tuple[float, float],
    theme: ThemeStyle,
    edge_fov_deg: float,
    content_fov_deg: float,
) -> None:
    del viewport_rect, theme, content_fov_deg

    if profile is None or not profile.samples:
        return
    samples = [sample for sample in profile.samples if float(sample.strength) > 0.0]
    if not samples:
        return

    seam_az_deg = (float(view_center[1]) + 180.0) % 360.0
    ordered = sorted(samples, key=lambda sample: (float(sample.azimuth_deg) - seam_az_deg) % 360.0)

    terrain_samples: list[tuple[float, float]] = []
    if terrain_profile_altaz:
        terrain_samples = sorted(
            ((float(alt_deg), float(az_deg) % 360.0) for alt_deg, az_deg in terrain_profile_altaz),
            key=lambda item: item[1],
        )
    if terrain_samples:
        terrain_az = np.asarray([az for _, az in terrain_samples], dtype=np.float64)
        terrain_alt = np.asarray([alt for alt, _ in terrain_samples], dtype=np.float64)
    else:
        terrain_az = np.asarray([float(sample.azimuth_deg) % 360.0 for sample in ordered], dtype=np.float64)
        terrain_alt = np.asarray([float(sample.horizon_alt_deg) for sample in ordered], dtype=np.float64)

    night_az = np.asarray([float(sample.azimuth_deg) % 360.0 for sample in ordered], dtype=np.float64)
    night_strengths = np.asarray([float(sample.strength) for sample in ordered], dtype=np.float64)
    if night_az.size < 2 or terrain_az.size < 2:
        return

    night_az_sorted = np.argsort(night_az)
    night_az = night_az[night_az_sorted]
    night_strengths = night_strengths[night_az_sorted]
    night_az_ext = np.concatenate([night_az[-1:] - 360.0, night_az, night_az[:1] + 360.0])
    night_strengths_ext = np.concatenate([night_strengths[-1:], night_strengths, night_strengths[:1]])
    terrain_az_sorted = np.argsort(terrain_az)
    terrain_az = terrain_az[terrain_az_sorted]
    terrain_alt = terrain_alt[terrain_az_sorted]
    terrain_az_ext = np.concatenate([terrain_az[-1:] - 360.0, terrain_az, terrain_az[:1] + 360.0])
    terrain_alt_ext = np.concatenate([terrain_alt[-1:], terrain_alt, terrain_alt[:1]])

    draw_az = np.arange(0.0, 360.0, NIGHT_LIGHTS_DRAW_AZIMUTH_STEP_DEG, dtype=np.float64)
    interpolated_strengths = np.interp(draw_az, night_az_ext, night_strengths_ext)
    horizon_alts = np.interp(draw_az, terrain_az_ext, terrain_alt_ext)

    projected_points: list[tuple[float, float]] = []
    strengths: list[float] = []
    center_alts: list[float] = []
    for az_deg, horizon_alt, strength in zip(draw_az.tolist(), horizon_alts.tolist(), interpolated_strengths.tolist()):
        center_alt = float(horizon_alt) + float(profile.band_center_offset_deg)
        try:
            nx, ny = altaz_to_normalized_xy(
                center_alt,
                float(az_deg),
                view_center,
                edge_fov_deg=float(edge_fov_deg),
            )
        except Exception:
            continue
        projected_points.append((float(nx), float(ny)))
        strengths.append(float(strength))
        center_alts.append(center_alt)

    if len(projected_points) < 2:
        return

    band_half_width_deg = float(profile.band_half_width_deg)
    width_px = _band_width_px(
        center_alt_deg=float(sum(center_alts) / len(center_alts)),
        band_half_width_deg=band_half_width_deg,
        view_center=view_center,
        geometry=geometry,
        edge_fov_deg=edge_fov_deg,
    )
    if width_px <= 0.0:
        return

    fill_rgb = tuple(int(value) for value in NIGHT_LIGHTS_GLOW_RGB)
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)

    point_index = 0
    for fragment in split_by_gaps(projected_points):
        if len(fragment) < 2:
            point_index += len(fragment)
            continue
        frag_strengths = strengths[point_index:point_index + len(fragment)]
        frag_points = [QPointF(*normalized_to_screen_xy(nx, ny, geometry)) for nx, ny in fragment]
        for start, end, start_strength, end_strength in zip(
            frag_points,
            frag_points[1:],
            frag_strengths,
            frag_strengths[1:],
        ):
            strength = max(float(start_strength), float(end_strength))
            if strength < NIGHT_LIGHTS_MIN_BRIGHTNESS:
                continue
            alpha = min(1.0, strength * NIGHT_LIGHTS_CORE_ALPHA_SCALE)
            color = QColor(*fill_rgb)

            color.setAlphaF(max(0.0, min(1.0, alpha * NIGHT_LIGHTS_OUTER_ALPHA_SCALE)))
            outer_pen = QPen(color, width_px * 3.0, Qt.PenStyle.SolidLine)
            outer_pen.setCosmetic(True)
            outer_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            outer_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(outer_pen)
            painter.drawLine(start, end)

            color.setAlphaF(max(0.0, min(1.0, alpha * NIGHT_LIGHTS_MID_ALPHA_SCALE)))
            mid_pen = QPen(color, width_px * 1.9, Qt.PenStyle.SolidLine)
            mid_pen.setCosmetic(True)
            mid_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            mid_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(mid_pen)
            painter.drawLine(start, end)

            color.setAlphaF(max(0.0, min(1.0, alpha * NIGHT_LIGHTS_CORE_ALPHA_SCALE)))
            core_pen = QPen(color, width_px * 1.0, Qt.PenStyle.SolidLine)
            core_pen.setCosmetic(True)
            core_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            core_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(core_pen)
            painter.drawLine(start, end)
        point_index += len(fragment)

    painter.restore()
