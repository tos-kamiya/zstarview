"""See-through Earth guide drawn below the horizon.

Each pixel below the horizon is traced as a ray through the Earth to the far
side of the globe. The surface point it reaches is tinted with a land or ocean
color, and its brightness follows the local solar altitude at that point, so
the sunlit part of the far side of the Earth appears under the observer's feet.
On the night side, city lights and thin white coastline outlines are drawn.
"""

from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import cast

import numpy as np
from PIL import Image, ImageDraw
from PySide6.QtCore import QLineF, QRectF, Qt
from PySide6.QtGui import QColor, QImage, QPainter, QPen

from ..astro import is_in_fov_vectorized
from ..night_lights_constants import NIGHT_LIGHTS_GLOW_RGB
from ..paths import (
    EARTH_GUIDE_DEFAULT_OPACITY,
    EARTH_GUIDE_LAND_FILE,
    EARTH_GUIDE_LAND_LOD_FILES,
    EARTH_GUIDE_LINE_COLOR,
    EARTH_GUIDE_OCEAN_COLOR,
    OverlayLayerStyle,
)
from ..types import ScreenGeometry, ViewerData
from .geometry import (
    _altaz_to_normalized_xy_vectorized,
    _normalized_to_screen_xy_vectorized,
)
from .ground_mask import interpolate_terrain_horizon_altitude, inverse_project_disc
from .qt_image import np_rgba_to_qimage

logger = logging.getLogger(__name__)
EARTH_GUIDE_NIGHT_LIGHTS_PATH = Path(__file__).resolve().parent.parent / "data" / "earth_guide_night_lights_2025.png"
EARTH_GUIDE_NIGHT_LIGHT_RGB = NIGHT_LIGHTS_GLOW_RGB
EARTH_GUIDE_NIGHT_LIGHT_ALPHA_SCALE = 1.75
EARTH_GUIDE_NIGHT_LIGHT_START_DEG = -8.0
EARTH_GUIDE_NIGHT_LIGHT_FULL_DEG = -14.0
EARTH_MEAN_RADIUS_M = 6_371_008.8

EARTH_GUIDE_LOD_COUNT = 4
EARTH_GUIDE_MASK_LOD = EARTH_GUIDE_LOD_COUNT - 1
EARTH_GUIDE_MASK_RES_DEG = 0.25
EARTH_GUIDE_MASK_SUPERSAMPLE = 2

EARTH_GUIDE_HORIZON_MARGIN_DEG = 1.0
# Soft edge against the horizon/terrain limit.
EARTH_GUIDE_HORIZON_FADE_DEG = 1.0
# Land near the observer is squeezed into thin streaks just below the horizon,
# so nothing is drawn within NEAR_CUTOFF_KM (surface distance). Outlines fade
# in until NEAR_CUTOFF_KM + NEAR_FADE_KM. The fill is cut sharply at the
# cutoff (FILL_NEAR_FADE_KM = 0): a fade would read as dimmer daylight.
EARTH_GUIDE_NEAR_CUTOFF_KM = 500.0
EARTH_GUIDE_NEAR_FADE_KM = 300.0
EARTH_GUIDE_FILL_NEAR_FADE_KM = 0.0

# Ocean alpha relative to land; lower this if the oceans look too busy.
EARTH_GUIDE_OCEAN_ALPHA_SCALE = 1.0

# Brightness at local solar altitude 0 deg; rises linearly to 1.0 at zenith.
# Kept high so the whole day side reads as lit, not only the subsolar area.
EARTH_GUIDE_DAY_MIN_LEVEL = 0.8
# Twilight fades from DAY_MIN_LEVEL at 0 deg to zero at this solar altitude.
# It overlaps the night outline fade-in, so the fill and outlines cross-fade.
EARTH_GUIDE_TWILIGHT_LIMIT_DEG = -18.0
# Brightness on the night side.
EARTH_GUIDE_NIGHT_LEVEL = 0.0

# Thin coastline outlines on the night side. They fade in from
# NIGHT_OUTLINE_START_DEG to full strength at NIGHT_OUTLINE_FULL_DEG.
EARTH_GUIDE_NIGHT_OUTLINE_RGB = (255, 255, 255)
EARTH_GUIDE_NIGHT_OUTLINE_START_DEG = -8.0
EARTH_GUIDE_NIGHT_OUTLINE_FULL_DEG = -14.0
# Outline alpha relative to the surface peak alpha.
EARTH_GUIDE_NIGHT_OUTLINE_ALPHA_SCALE = 0.8
EARTH_GUIDE_NIGHT_OUTLINE_WIDTH_PX = 1.0
EARTH_GUIDE_NIGHT_OUTLINE_STEP_DEG = 0.5
EARTH_GUIDE_NIGHT_OUTLINE_ALPHA_LEVELS = 8

EARTH_GUIDE_RENDER_SCALE = 0.5


@dataclass(frozen=True)
class EarthGuideRing:
    source_name: str
    label_name: str | None
    points_lonlat_deg: np.ndarray
    points_xyz: np.ndarray
    approx_area_deg2: float | None = None


def _lonlat_to_unit_xyz(lon_deg: float, lat_deg: float) -> tuple[float, float, float]:
    lon_rad = math.radians(lon_deg)
    lat_rad = math.radians(lat_deg)
    cos_lat = math.cos(lat_rad)
    return (
        cos_lat * math.cos(lon_rad),
        cos_lat * math.sin(lon_rad),
        math.sin(lat_rad),
    )


def _observer_basis(
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return observer position and east/north/up unit vectors (Earth-fixed)."""
    lat_rad = math.radians(observer_lat_deg)
    lon_rad = math.radians(observer_lon_deg)
    up = np.array(
        (
            math.cos(lat_rad) * math.cos(lon_rad),
            math.cos(lat_rad) * math.sin(lon_rad),
            math.sin(lat_rad),
        ),
        dtype=np.float64,
    )
    east = np.array((-math.sin(lon_rad), math.cos(lon_rad), 0.0), dtype=np.float64)
    north = np.array(
        (
            -math.sin(lat_rad) * math.cos(lon_rad),
            -math.sin(lat_rad) * math.sin(lon_rad),
            math.cos(lat_rad),
        ),
        dtype=np.float64,
    )
    radius_scale = 1.0 + (max(0.0, float(observer_height_m)) / EARTH_MEAN_RADIUS_M)
    origin = up * radius_scale
    return origin, east, north, up


def _observer_horizon_dip_deg(observer_height_m: float) -> float:
    height_m = max(0.0, float(observer_height_m))
    if height_m <= 0.0:
        return 0.0
    ratio = EARTH_MEAN_RADIUS_M / (EARTH_MEAN_RADIUS_M + height_m)
    return math.degrees(math.acos(max(-1.0, min(1.0, ratio))))


def visible_altitude_limit_deg(
    az_deg: np.ndarray,
    observer_height_m: float,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> np.ndarray:
    """Return the highest altitude where the Earth guide may appear."""
    horizon_limit = -(
        _observer_horizon_dip_deg(observer_height_m) + EARTH_GUIDE_HORIZON_MARGIN_DEG
    )
    az = np.asarray(az_deg, dtype=np.float32)
    if not terrain_profile_altaz:
        return np.full(az.shape, horizon_limit, dtype=np.float32)
    terrain_limit = interpolate_terrain_horizon_altitude(az, terrain_profile_altaz)
    return np.minimum(terrain_limit, np.float32(horizon_limit))


def altaz_to_earth_fixed(
    alt_deg: np.ndarray | float,
    az_deg: np.ndarray | float,
    *,
    east: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
) -> np.ndarray:
    """Convert observer alt/az directions to Earth-fixed unit vectors."""
    alt = np.radians(np.asarray(alt_deg, dtype=np.float64))
    az = np.radians(np.asarray(az_deg, dtype=np.float64))
    cos_alt = np.cos(alt)
    e = (cos_alt * np.sin(az))[..., None]
    n = (cos_alt * np.cos(az))[..., None]
    u = np.sin(alt)[..., None]
    return e * east + n * north + u * up


def earth_guide_sunlight_level(sin_sun_alt: np.ndarray) -> np.ndarray:
    """Return the surface brightness factor for the local solar altitude.

    The day side runs from ``EARTH_GUIDE_DAY_MIN_LEVEL`` at sunrise/sunset to
    1.0 with the Sun at the zenith. Twilight fades smoothly to the night level
    at ``EARTH_GUIDE_TWILIGHT_LIMIT_DEG``.
    """
    sin_h = np.clip(np.asarray(sin_sun_alt, dtype=np.float64), -1.0, 1.0)
    day = EARTH_GUIDE_DAY_MIN_LEVEL + (1.0 - EARTH_GUIDE_DAY_MIN_LEVEL) * sin_h
    h_deg = np.degrees(np.arcsin(sin_h))
    limit = float(EARTH_GUIDE_TWILIGHT_LIMIT_DEG)
    t = np.clip((h_deg - limit) / max(1.0e-6, -limit), 0.0, 1.0)
    smooth = t * t * (3.0 - 2.0 * t)
    twilight = EARTH_GUIDE_NIGHT_LEVEL + (
        EARTH_GUIDE_DAY_MIN_LEVEL - EARTH_GUIDE_NIGHT_LEVEL
    ) * smooth
    return np.where(sin_h >= 0.0, day, twilight)


@lru_cache(maxsize=EARTH_GUIDE_LOD_COUNT)
def load_earth_guide_rings(lod: int = 0) -> tuple[EarthGuideRing, ...]:
    lod_index = max(0, min(EARTH_GUIDE_LOD_COUNT - 1, int(lod)))
    if lod_index == 0 and not Path(EARTH_GUIDE_LAND_LOD_FILES[0]).is_file():
        path = Path(EARTH_GUIDE_LAND_FILE)
    else:
        path = Path(EARTH_GUIDE_LAND_LOD_FILES[lod_index])
    payload = json.loads(path.read_text(encoding="utf-8"))
    raw_rings = payload.get("rings", [])
    rings: list[EarthGuideRing] = []
    for item in raw_rings:
        if not isinstance(item, dict):
            continue
        raw_points = item.get("points_lonlat_deg")
        if not isinstance(raw_points, list):
            continue
        lonlat: list[tuple[float, float]] = []
        xyz: list[tuple[float, float, float]] = []
        for pair in raw_points:
            if not isinstance(pair, list) or len(pair) < 2:
                continue
            lon_deg = float(pair[0])
            lat_deg = float(pair[1])
            lonlat.append((lon_deg, lat_deg))
            xyz.append(_lonlat_to_unit_xyz(lon_deg, lat_deg))
        if len(lonlat) < 3:
            continue
        approx_area_deg2: float | None = None
        raw_area = item.get("approx_area_deg2")
        if isinstance(raw_area, (int, float)):
            approx_area_deg2 = float(raw_area)
        rings.append(
            EarthGuideRing(
                source_name=str(item.get("source_name", "")),
                label_name=(str(item["label_name"]) if item.get("label_name") else None),
                points_lonlat_deg=np.asarray(lonlat, dtype=np.float64),
                points_xyz=np.asarray(xyz, dtype=np.float64),
                approx_area_deg2=approx_area_deg2,
            )
        )
    return tuple(rings)


def build_land_mask(
    rings: tuple[EarthGuideRing, ...],
    *,
    res_deg: float = EARTH_GUIDE_MASK_RES_DEG,
    supersample: int = EARTH_GUIDE_MASK_SUPERSAMPLE,
) -> np.ndarray:
    """Rasterize land rings into a lon/lat grid of land fractions (0..1).

    Row 0 is the northernmost band and column 0 starts at longitude -180.
    Rings are filled as simple polygons in the lon/lat plane; the bundled
    Natural Earth rings are already split at the antimeridian.
    """
    ss = max(1, int(supersample))
    cell = float(res_deg) / ss
    width = int(round(360.0 / cell))
    height = int(round(180.0 / cell))
    image = Image.new("L", (width, height), 0)
    draw = ImageDraw.Draw(image)
    for ring in rings:
        pts = ring.points_lonlat_deg
        if len(pts) < 3:
            continue
        xy = [
            ((float(lon) + 180.0) / cell, (90.0 - float(lat)) / cell)
            for lon, lat in pts
        ]
        draw.polygon(xy, fill=255)
    fine = np.asarray(image, dtype=np.float32) / 255.0
    if ss > 1:
        fine = fine.reshape(height // ss, ss, width // ss, ss).mean(axis=(1, 3))
    return fine.astype(np.float32)


@lru_cache(maxsize=1)
def load_land_mask() -> np.ndarray:
    """Return the cached land-fraction grid used by the Earth guide."""
    return build_land_mask(load_earth_guide_rings(EARTH_GUIDE_MASK_LOD))


@lru_cache(maxsize=1)
def load_earth_guide_night_lights() -> np.ndarray | None:
    """Load the bundled, display-scaled global night-light grid once."""
    try:
        with Image.open(EARTH_GUIDE_NIGHT_LIGHTS_PATH) as image:
            if image.size not in ((3600, 1800), (7200, 3600)) or image.mode != "L":
                raise ValueError("unexpected night-light image format")
            return np.asarray(image, dtype=np.uint8).copy()
    except (OSError, ValueError) as exc:
        logger.warning("Cannot load Earth guide night lights: %s", exc)
        return None


def sample_land_fraction(
    mask: np.ndarray, lon_deg: np.ndarray, lat_deg: np.ndarray
) -> np.ndarray:
    """Bilinearly sample a land mask at the given longitudes and latitudes."""
    rows, cols = mask.shape
    res_lon = 360.0 / cols
    res_lat = 180.0 / rows
    x = (np.asarray(lon_deg, dtype=np.float64) + 180.0) / res_lon - 0.5
    y = (90.0 - np.asarray(lat_deg, dtype=np.float64)) / res_lat - 0.5
    y = np.clip(y, 0.0, rows - 1.0)
    x0 = np.floor(x)
    y0 = np.clip(np.floor(y), 0, rows - 1)
    fx = x - x0
    fy = y - y0
    x0i = x0.astype(np.int64) % cols
    x1i = (x0i + 1) % cols
    y0i = y0.astype(np.int64)
    y1i = np.minimum(y0i + 1, rows - 1)
    top = mask[y0i, x0i] * (1.0 - fx) + mask[y0i, x1i] * fx
    bottom = mask[y1i, x0i] * (1.0 - fx) + mask[y1i, x1i] * fx
    return np.asarray(top * (1.0 - fy) + bottom * fy, dtype=np.float32)


def far_surface_points(
    dirs: np.ndarray,
    *,
    origin: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Trace rays from the observer to where they leave the unit sphere.

    Returns ``(points, valid)``. A ray that misses the globe (possible just
    below the horizon for an elevated observer) is marked invalid.
    """
    b = dirs @ origin
    c = float(origin @ origin) - 1.0
    disc = b * b - c
    valid = disc >= 0.0
    t = -b + np.sqrt(np.maximum(disc, 0.0))
    points = origin + t[..., None] * dirs
    norms = np.linalg.norm(points, axis=-1, keepdims=True)
    points = points / np.maximum(norms, 1.0e-12)
    return points, valid


def earth_guide_peak_alpha(opacity: float, visibility_boost: float = 1.0) -> float:
    """Return the surface alpha under full sunlight."""
    opacity = max(0.0, float(opacity))
    if opacity <= 0.0:
        return 0.0
    alpha = opacity * max(1.0, float(visibility_boost))
    return max(0.0, min(1.0, alpha))


def earth_guide_line_alpha(opacity: float) -> float:
    """Return the alpha curve shared with other guide lines (never-rises)."""
    opacity = max(0.0, min(1.0, float(opacity)))
    return max(opacity, min(1.0, 0.32 + (opacity * 0.78)))


def _smoothstep01(t: np.ndarray) -> np.ndarray:
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def earth_guide_near_fade(
    points: np.ndarray,
    up: np.ndarray,
    *,
    fade_km: float = EARTH_GUIDE_NEAR_FADE_KM,
) -> np.ndarray:
    """Return 0..1 visibility by surface distance from the observer.

    Zero within ``EARTH_GUIDE_NEAR_CUTOFF_KM``; rises smoothly over
    ``fade_km`` beyond it, or steps straight to 1.0 when ``fade_km`` is 0.
    """
    cos_c = np.clip(points @ up, -1.0, 1.0)
    distance_km = np.arccos(cos_c) * (EARTH_MEAN_RADIUS_M / 1000.0)
    if fade_km <= 0.0:
        return (distance_km >= EARTH_GUIDE_NEAR_CUTOFF_KM).astype(np.float64)
    return _smoothstep01(
        (distance_km - EARTH_GUIDE_NEAR_CUTOFF_KM) / float(fade_km)
    )


def render_earth_guide_rgba(
    width: int,
    height: int,
    *,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    sun_altaz: tuple[float, float] | None,
    terrain_profile_altaz: list[tuple[float, float]] | None = None,
    earth_guide_opacity: float = EARTH_GUIDE_DEFAULT_OPACITY,
    visibility_boost: float = 1.0,
    land_rgb: tuple[int, int, int] = EARTH_GUIDE_LINE_COLOR,
    ocean_rgb: tuple[int, int, int] = EARTH_GUIDE_OCEAN_COLOR,
    alpha_scale: float = 1.0,
) -> np.ndarray | None:
    """Render the Earth guide into a straight-alpha RGBA array.

    Returns ``None`` when nothing would be drawn.
    """
    if width <= 0 or height <= 0:
        return None
    peak_alpha = earth_guide_peak_alpha(earth_guide_opacity, visibility_boost)
    peak_alpha *= max(0.0, float(alpha_scale))
    if peak_alpha <= 0.0:
        return None

    alt, az, inside = inverse_project_disc(
        width,
        height,
        geometry,
        viewer_data.view_center,
        edge_fov_deg=float(viewer_data.edge_fov_deg),
        content_fov_deg=float(viewer_data.content_fov_deg),
    )
    if alt.size == 0:
        return None

    observer_height_m = float(viewer_data.observer_height_m)
    limit = visible_altitude_limit_deg(az, observer_height_m, terrain_profile_altaz)
    below = alt < limit
    if not np.any(below):
        return None
    alt_b = alt[below]
    az_b = az[below]
    horizon_fade = _smoothstep01(
        (limit[below] - alt_b) / max(1.0e-6, EARTH_GUIDE_HORIZON_FADE_DEG)
    )

    origin, east, north, up = _observer_basis(
        float(viewer_data.lat_deg),
        float(viewer_data.lon_deg),
        observer_height_m,
    )
    dirs = altaz_to_earth_fixed(alt_b, az_b, east=east, north=north, up=up)
    points, valid = far_surface_points(dirs, origin=origin)

    if sun_altaz is None:
        light = np.full(points.shape[0], EARTH_GUIDE_DAY_MIN_LEVEL, dtype=np.float64)
        sun_strength = None
    else:
        sun_dir = altaz_to_earth_fixed(
            float(sun_altaz[0]), float(sun_altaz[1]), east=east, north=north, up=up
        )
        sun_strength = points @ sun_dir
        light = earth_guide_sunlight_level(sun_strength)

    lon = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    lat = np.degrees(np.arcsin(np.clip(points[:, 2], -1.0, 1.0)))
    land = sample_land_fraction(load_land_mask(), lon, lat).astype(np.float64)
    ocean = 1.0 - land

    land_a = land
    ocean_a = ocean * EARTH_GUIDE_OCEAN_ALPHA_SCALE
    coverage = land_a + ocean_a
    horizon_fade = horizon_fade * earth_guide_near_fade(
        points, up, fade_km=EARTH_GUIDE_FILL_NEAR_FADE_KM
    )
    alpha = peak_alpha * light * horizon_fade * coverage * valid
    safe = np.maximum(coverage, 1.0e-6)[:, None]
    rgb = (
        land_a[:, None] * np.asarray(land_rgb, dtype=np.float64)
        + ocean_a[:, None] * np.asarray(ocean_rgb, dtype=np.float64)
    ) / safe

    if sun_strength is not None:
        night_weight = earth_guide_night_light_weight(sun_strength)
        night_pixels = (night_weight > 0.0) & valid & (horizon_fade > 0.0)
        if np.any(night_pixels):
            night_grid = load_earth_guide_night_lights()
            if night_grid is not None:
                intensity = sample_land_fraction(night_grid, lon[night_pixels], lat[night_pixels]) / 255.0
                night_alpha = np.clip(
                    peak_alpha * EARTH_GUIDE_NIGHT_LIGHT_ALPHA_SCALE
                    * horizon_fade[night_pixels] * night_weight[night_pixels] * intensity,
                    0.0, 1.0,
                )
                old_alpha = alpha[night_pixels]
                combined_alpha = night_alpha + old_alpha * (1.0 - night_alpha)
                rgb[night_pixels] = (
                    np.asarray(EARTH_GUIDE_NIGHT_LIGHT_RGB) * night_alpha[:, None]
                    + rgb[night_pixels] * (old_alpha * (1.0 - night_alpha))[:, None]
                ) / np.maximum(combined_alpha[:, None], 1.0e-6)
                alpha[night_pixels] = combined_alpha

    out = np.zeros((height, width, 4), dtype=np.uint8)
    inside_idx = np.flatnonzero(inside.ravel())
    target = inside_idx[below]
    flat = out.reshape(-1, 4)
    flat[target, :3] = np.clip(np.rint(rgb), 0, 255).astype(np.uint8)
    flat[target, 3] = np.clip(np.rint(alpha * 255.0), 0, 255).astype(np.uint8)
    if not np.any(flat[:, 3]):
        return None
    return out


def earth_guide_night_light_weight(sin_sun_alt: np.ndarray) -> np.ndarray:
    """Fade city lights in between solar altitudes -8 and -14 degrees."""
    sin_h = np.clip(np.asarray(sin_sun_alt, dtype=np.float64), -1.0, 1.0)
    h_deg = np.degrees(np.arcsin(sin_h))
    t = np.clip(
        (EARTH_GUIDE_NIGHT_LIGHT_START_DEG - h_deg)
        / (EARTH_GUIDE_NIGHT_LIGHT_START_DEG - EARTH_GUIDE_NIGHT_LIGHT_FULL_DEG),
        0.0, 1.0,
    )
    return t * t * (3.0 - 2.0 * t)


def earth_guide_night_outline_weight(sin_sun_alt: np.ndarray) -> np.ndarray:
    """Return the night-side outline strength (0..1) for the local solar altitude.

    Zero above ``EARTH_GUIDE_NIGHT_OUTLINE_START_DEG``, rising smoothly to 1.0
    at ``EARTH_GUIDE_NIGHT_OUTLINE_FULL_DEG``.
    """
    sin_h = np.clip(np.asarray(sin_sun_alt, dtype=np.float64), -1.0, 1.0)
    h_deg = np.degrees(np.arcsin(sin_h))
    start = float(EARTH_GUIDE_NIGHT_OUTLINE_START_DEG)
    full = float(EARTH_GUIDE_NIGHT_OUTLINE_FULL_DEG)
    t = np.clip((start - h_deg) / max(1.0e-6, start - full), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _is_seam_segment(lonlat0: np.ndarray, lonlat1: np.ndarray) -> bool:
    """Return True for ring edges that are data seams, not coastlines.

    The bundled rings are split along the antimeridian and Antarctica is
    closed along the South Pole; those edges must not be drawn as coastline.
    """
    on_meridian = abs(abs(float(lonlat0[0])) - 180.0) < 1.0e-6 and abs(
        abs(float(lonlat1[0])) - 180.0
    ) < 1.0e-6
    on_pole = float(lonlat0[1]) <= -90.0 + 1.0e-6 and float(lonlat1[1]) <= -90.0 + 1.0e-6
    return on_meridian or on_pole


@lru_cache(maxsize=EARTH_GUIDE_LOD_COUNT)
def load_coastline_segments(
    lod: int = EARTH_GUIDE_MASK_LOD,
) -> tuple[np.ndarray, np.ndarray]:
    """Return densified coastline segments as Earth-fixed unit vectors.

    Returns ``(start_xyz, end_xyz)``, each of shape (M, 3). Ring edges are
    split so that no segment spans more than ``EARTH_GUIDE_NIGHT_OUTLINE_STEP_DEG``.
    """
    step_rad = math.radians(max(1.0e-3, EARTH_GUIDE_NIGHT_OUTLINE_STEP_DEG))
    starts: list[np.ndarray] = []
    ends: list[np.ndarray] = []
    for ring in load_earth_guide_rings(lod):
        lonlat = ring.points_lonlat_deg
        xyz = ring.points_xyz
        count = len(xyz)
        for index in range(count):
            nxt = (index + 1) % count
            if _is_seam_segment(lonlat[index], lonlat[nxt]):
                continue
            a = xyz[index]
            b = xyz[nxt]
            angle = math.acos(max(-1.0, min(1.0, float(a @ b))))
            pieces = max(1, int(math.ceil(angle / step_rad)))
            t = np.linspace(0.0, 1.0, pieces + 1)[:, None]
            pts = (1.0 - t) * a + t * b
            pts /= np.maximum(np.linalg.norm(pts, axis=1, keepdims=True), 1.0e-12)
            starts.append(pts[:-1])
            ends.append(pts[1:])
    if not starts:
        empty = np.empty((0, 3), dtype=np.float64)
        return empty, empty
    return np.concatenate(starts), np.concatenate(ends)


def _project_surface_points_to_screen(
    points: np.ndarray,
    *,
    origin: np.ndarray,
    east: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Project far-side surface points to screen; return (xy, visibility 0..1)."""
    rays = points - origin
    rays /= np.maximum(np.linalg.norm(rays, axis=1, keepdims=True), 1.0e-12)
    alt = np.degrees(np.arcsin(np.clip(rays @ up, -1.0, 1.0)))
    az = np.degrees(np.arctan2(rays @ east, rays @ north)) % 360.0
    limit = visible_altitude_limit_deg(
        az, float(viewer_data.observer_height_m), terrain_profile_altaz
    ).astype(np.float64)
    visibility = _smoothstep01(
        (limit - alt) / max(1.0e-6, EARTH_GUIDE_HORIZON_FADE_DEG)
    ) * earth_guide_near_fade(points, up)
    visibility = visibility * is_in_fov_vectorized(
        alt, az, viewer_data.view_center, fov_deg=float(viewer_data.content_fov_deg)
    )
    nx, ny = _altaz_to_normalized_xy_vectorized(
        alt, az, viewer_data.view_center, edge_fov_deg=float(viewer_data.edge_fov_deg)
    )
    sx, sy = _normalized_to_screen_xy_vectorized(nx, ny, geometry)
    return np.column_stack((sx, sy)), visibility


def compute_night_outline_segments(
    *,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    sun_altaz: tuple[float, float] | None,
    terrain_profile_altaz: list[tuple[float, float]] | None = None,
    lod: int = EARTH_GUIDE_MASK_LOD,
    night_only: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return coastline segments in screen space.

    Returns ``(start_xy, end_xy, weight)``; weight combines night strength,
    horizon fade and field-of-view visibility. With ``night_only`` (the
    default) only the night side is returned, and nothing without a Sun
    position. With ``night_only=False`` all visible coastlines are returned
    regardless of sunlight; fast mode uses this in place of the fill.
    """
    empty = (np.empty((0, 2)), np.empty((0, 2)), np.empty((0,)))
    if night_only and sun_altaz is None:
        return empty
    seg_a, seg_b = load_coastline_segments(lod)
    if len(seg_a) == 0:
        return empty
    origin, east, north, up = _observer_basis(
        float(viewer_data.lat_deg),
        float(viewer_data.lon_deg),
        float(viewer_data.observer_height_m),
    )
    if night_only:
        assert sun_altaz is not None
        sun_dir = altaz_to_earth_fixed(
            float(sun_altaz[0]), float(sun_altaz[1]), east=east, north=north, up=up
        )
        night_a = earth_guide_night_outline_weight(seg_a @ sun_dir)
        night_b = earth_guide_night_outline_weight(seg_b @ sun_dir)
    else:
        night_a = np.ones(len(seg_a), dtype=np.float64)
        night_b = night_a
    candidate = (night_a > 0.0) & (night_b > 0.0)
    if not np.any(candidate):
        return empty
    seg_a = seg_a[candidate]
    seg_b = seg_b[candidate]
    xy_a, vis_a = _project_surface_points_to_screen(
        seg_a,
        origin=origin,
        east=east,
        north=north,
        up=up,
        geometry=geometry,
        viewer_data=viewer_data,
        terrain_profile_altaz=terrain_profile_altaz,
    )
    xy_b, vis_b = _project_surface_points_to_screen(
        seg_b,
        origin=origin,
        east=east,
        north=north,
        up=up,
        geometry=geometry,
        viewer_data=viewer_data,
        terrain_profile_altaz=terrain_profile_altaz,
    )
    weight = np.minimum(night_a[candidate] * vis_a, night_b[candidate] * vis_b)
    # Drop segments that jump across the screen (e.g. through the observer).
    span = np.hypot(*(xy_b - xy_a).T)
    weight = np.where(span <= max(8.0, 0.25 * float(geometry.radius)), weight, 0.0)
    keep = weight > 0.0
    return xy_a[keep], xy_b[keep], weight[keep]


def _draw_night_outlines(
    painter: QPainter,
    start_xy: np.ndarray,
    end_xy: np.ndarray,
    weight: np.ndarray,
    *,
    peak_alpha: float,
) -> None:
    if len(weight) == 0 or peak_alpha <= 0.0:
        return
    levels = max(1, int(EARTH_GUIDE_NIGHT_OUTLINE_ALPHA_LEVELS))
    bucket = np.clip(np.ceil(weight * levels).astype(np.int64), 1, levels)
    painter.save()
    try:
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        for level in range(1, levels + 1):
            idx = np.flatnonzero(bucket == level)
            if idx.size == 0:
                continue
            color = QColor(*EARTH_GUIDE_NIGHT_OUTLINE_RGB)
            color.setAlphaF(max(0.0, min(1.0, peak_alpha * level / levels)))
            pen = QPen(color, EARTH_GUIDE_NIGHT_OUTLINE_WIDTH_PX)
            pen.setCosmetic(True)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(pen)
            painter.drawLines(
                [
                    QLineF(
                        float(start_xy[i, 0]),
                        float(start_xy[i, 1]),
                        float(end_xy[i, 0]),
                        float(end_xy[i, 1]),
                    )
                    for i in idx
                ]
            )
    finally:
        painter.restore()


def _draw_surface_fill(
    painter: QPainter,
    *,
    target_size: tuple[int, int],
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    sun_altaz: tuple[float, float] | None,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    earth_guide_opacity: float,
    visibility_boost: float,
    land_rgb: tuple[int, int, int],
    alpha_scale: float,
) -> None:
    full_w, full_h = target_size
    scale = EARTH_GUIDE_RENDER_SCALE
    w = max(1, int(math.ceil(full_w * scale)))
    h = max(1, int(math.ceil(full_h * scale)))
    # Low-res pixel j covers full-res pixels [j/scale, (j+1)/scale); sample at
    # its center, expressed in the full-res pixel-index convention.
    offset = 0.5 - (0.5 * scale)
    scaled_geometry = ScreenGeometry(
        center=cast(
            tuple[int, int],
            (
                (float(geometry.center[0]) * scale) - offset,
                (float(geometry.center[1]) * scale) - offset,
            ),
        ),
        radius=cast(int, float(geometry.radius) * scale),
    )
    rgba = render_earth_guide_rgba(
        w,
        h,
        geometry=scaled_geometry,
        viewer_data=viewer_data,
        sun_altaz=sun_altaz,
        terrain_profile_altaz=terrain_profile_altaz,
        earth_guide_opacity=earth_guide_opacity,
        visibility_boost=visibility_boost,
        land_rgb=land_rgb,
        alpha_scale=alpha_scale,
    )
    if rgba is None:
        return
    image: QImage = np_rgba_to_qimage(rgba)
    painter.save()
    try:
        painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
        painter.drawImage(QRectF(0.0, 0.0, w / scale, h / scale), image)
    finally:
        painter.restore()


def draw_earth_guide(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    sun_altaz: tuple[float, float] | None = None,
    terrain_profile_altaz: list[tuple[float, float]] | None = None,
    earth_guide_opacity: float = EARTH_GUIDE_DEFAULT_OPACITY,
    visibility_boost: float = 1.0,
    fast_mode: bool = False,
    layer_style: OverlayLayerStyle | None = None,
    target_size: tuple[int, int] | None = None,
) -> None:
    """Draw the sunlit see-through Earth below the horizon.

    Normally the day side is filled (at reduced resolution, scaled up with
    smoothing) and the night side gets thin white coastline outlines. In
    ``fast_mode`` the fill is skipped and all visible coastlines are drawn as
    outlines from the coarsest LOD instead.
    """
    if target_size is None:
        device = painter.device()
        if device is None:
            return
        target_size = (int(device.width()), int(device.height()))
    if target_size[0] <= 0 or target_size[1] <= 0:
        return
    land_rgb = EARTH_GUIDE_LINE_COLOR if layer_style is None else layer_style.rgb
    alpha_scale = 1.0 if layer_style is None else float(layer_style.alpha_scale)
    if not fast_mode:
        _draw_surface_fill(
            painter,
            target_size=target_size,
            geometry=geometry,
            viewer_data=viewer_data,
            sun_altaz=sun_altaz,
            terrain_profile_altaz=terrain_profile_altaz,
            earth_guide_opacity=earth_guide_opacity,
            visibility_boost=visibility_boost,
            land_rgb=land_rgb,
            alpha_scale=alpha_scale,
        )

    outline_alpha = (
        earth_guide_peak_alpha(earth_guide_opacity, visibility_boost)
        * max(0.0, alpha_scale)
        * EARTH_GUIDE_NIGHT_OUTLINE_ALPHA_SCALE
    )
    if outline_alpha <= 0.0:
        return
    start_xy, end_xy, weight = compute_night_outline_segments(
        geometry=geometry,
        viewer_data=viewer_data,
        sun_altaz=sun_altaz,
        terrain_profile_altaz=terrain_profile_altaz,
        lod=0 if fast_mode else EARTH_GUIDE_MASK_LOD,
        night_only=not fast_mode,
    )
    _draw_night_outlines(
        painter, start_xy, end_xy, weight, peak_alpha=min(1.0, outline_alpha)
    )
