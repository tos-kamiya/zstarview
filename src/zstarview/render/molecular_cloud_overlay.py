"""Minimal runtime renderer for the locally prepared AKARI display cache."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import Galactic, SkyCoord
from astropy.time import Time
from astropy.coordinates import EarthLocation

from ..astro import build_icrs_to_altaz_matrix
from ..paths import CACHE_PATH
from ..types import ScreenGeometry, ViewProjection
from .sky_disc import _inverse_project_disc, _smoothstep

MOLECULAR_CLOUD_CACHE = (
    Path(CACHE_PATH)
    / "molecular-cloud"
    / "akari-far-infrared-all-sky"
    / "release-1"
    / "schema-1"
    / "akari-galactic-display.npz"
)
MOLECULAR_CLOUD_MAX_SUN_ALT_DEG = -4.0
MOLECULAR_CLOUD_FULL_SUN_ALT_DEG = -12.0
MOLECULAR_CLOUD_OPACITY = 0.15
MOLECULAR_CLOUD_GAMMA = 0.7
MOLECULAR_CLOUD_VALUE_KNEE = 1.0
# Set to "jwst" to preview the JWST-inspired wavelength palette.
MOLECULAR_CLOUD_PALETTE = "akari"

_JWST_BLUE = (0.25, 0.35, 1.00)
_JWST_GREEN = (0.65, 0.75, 0.20)
_JWST_RED = (1.00, 0.30, 0.10)

_ICRS_BASIS = SkyCoord(
    ra=np.array([0.0, 90.0, 0.0]) * u.deg,
    dec=np.array([0.0, 0.0, 90.0]) * u.deg,
    frame="icrs",
)
_ICRS_TO_GALACTIC = _ICRS_BASIS.transform_to(Galactic).cartesian.xyz.to_value(u.one)


def is_molecular_cloud_cache_available() -> bool:
    return MOLECULAR_CLOUD_CACHE.is_file()


def _apply_molecular_cloud_value_knee(rgb: np.ndarray) -> np.ndarray:
    """Compress HSV value while preserving hue and saturation."""
    clipped = np.clip(rgb, 0.0, 1.0)
    value = np.max(clipped, axis=1, keepdims=True)
    mapped_value = value / (1.0 + MOLECULAR_CLOUD_VALUE_KNEE * value)
    scale = np.divide(mapped_value, value, out=np.zeros_like(value), where=value > 0.0)
    return clipped * scale


@lru_cache(maxsize=1)
def _load_display_asset(path: str, mtime_ns: int) -> tuple[np.ndarray, tuple[int, ...]] | None:
    del mtime_ns
    try:
        with np.load(path) as archive:
            data = np.asarray(archive["data"], dtype=np.float32) / 65535.0
            bands = tuple(int(value) for value in np.asarray(archive["bands"]).tolist())
    except (FileNotFoundError, OSError, KeyError, ValueError):
        return None
    if data.ndim != 3 or data.shape[0] != len(bands):
        return None
    return data, bands


def _sample_galactic_asset(
    data: np.ndarray,
    bands: tuple[int, ...],
    *,
    alt_deg: np.ndarray,
    az_deg: np.ndarray,
    time_obj: Time,
    observer_lat_deg: float,
    observer_lon_deg: float,
) -> np.ndarray:
    location = EarthLocation(
        lat=float(observer_lat_deg) * u.deg,
        lon=float(observer_lon_deg) * u.deg,
    )
    alt_rad = np.radians(alt_deg)
    az_rad = np.radians(az_deg)
    altaz_vectors = np.column_stack(
        (
            np.cos(alt_rad) * np.cos(az_rad),
            np.cos(alt_rad) * np.sin(az_rad),
            np.sin(alt_rad),
        )
    )
    icrs_to_altaz = build_icrs_to_altaz_matrix(time_obj, location)
    icrs_vectors = altaz_vectors @ icrs_to_altaz
    galactic_vectors = icrs_vectors @ _ICRS_TO_GALACTIC.T
    gal_lon = np.degrees(np.arctan2(galactic_vectors[:, 1], galactic_vectors[:, 0])) % 360.0
    gal_lat = np.degrees(np.arcsin(np.clip(galactic_vectors[:, 2], -1.0, 1.0)))
    width = data.shape[2]
    height = data.shape[1]
    x = np.floor(gal_lon / 360.0 * width).astype(np.int64) % width
    y = np.clip(np.rint((90.0 - gal_lat) / 180.0 * (height - 1)), 0, height - 1).astype(np.int64)

    channels = {band: data[index, y, x] for index, band in enumerate(bands)}
    red = channels.get(160, np.zeros_like(gal_lon, dtype=np.float32))
    green = channels.get(140, red)
    blue = channels.get(90, green)
    if MOLECULAR_CLOUD_PALETTE == "jwst":
        return np.column_stack(
            (
                _JWST_BLUE[0] * blue + _JWST_GREEN[0] * green + _JWST_RED[0] * red,
                _JWST_BLUE[1] * blue + _JWST_GREEN[1] * green + _JWST_RED[1] * red,
                _JWST_BLUE[2] * blue + _JWST_GREEN[2] * green + _JWST_RED[2] * red,
            )
        )
    return np.column_stack((red, green, blue))


def render_molecular_cloud_overlay(
    *,
    width: int,
    height: int,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    content_fov_deg: float,
    sun_alt_deg: float | None,
    time_obj: Time | None,
    observer_lat_deg: float | None,
    observer_lon_deg: float | None,
    opacity: float = MOLECULAR_CLOUD_OPACITY,
) -> np.ndarray | None:
    """Return an additive RGB overlay sampled from the local AKARI asset."""
    if (
        float(opacity) <= 0.0
        or not is_molecular_cloud_cache_available()
        or sun_alt_deg is None
        or time_obj is None
        or observer_lat_deg is None
        or observer_lon_deg is None
        or float(sun_alt_deg) > MOLECULAR_CLOUD_MAX_SUN_ALT_DEG
    ):
        return None
    try:
        stat = MOLECULAR_CLOUD_CACHE.stat()
    except OSError:
        return None
    asset = _load_display_asset(str(MOLECULAR_CLOUD_CACHE), int(stat.st_mtime_ns))
    if asset is None:
        return None
    data, bands = asset
    alt, az, inside = _inverse_project_disc(
        int(width),
        int(height),
        geometry,
        ViewProjection(
            view_center=(float(view_center[0]), float(view_center[1])),
            edge_fov_deg=float(edge_fov_deg),
            content_fov_deg=float(content_fov_deg),
        ),
    )
    overlay = np.zeros((int(height), int(width), 4), dtype=np.uint8)
    if alt.size == 0:
        return overlay
    rgb = _sample_galactic_asset(
        data,
        bands,
        alt_deg=alt,
        az_deg=az,
        time_obj=time_obj,
        observer_lat_deg=float(observer_lat_deg),
        observer_lon_deg=float(observer_lon_deg),
    )
    rgb = np.power(np.clip(rgb, 0.0, 1.0), MOLECULAR_CLOUD_GAMMA)
    rgb = _apply_molecular_cloud_value_knee(rgb)
    night_amount = 1.0 - _smoothstep(
        MOLECULAR_CLOUD_FULL_SUN_ALT_DEG,
        MOLECULAR_CLOUD_MAX_SUN_ALT_DEG,
        float(sun_alt_deg),
    )
    alpha = float(np.clip(opacity, 0.0, 1.0)) * night_amount
    overlay[..., :3][inside] = np.clip(np.round(rgb * alpha * 255.0), 0, 255).astype(np.uint8)
    overlay[..., 3][inside] = 255
    return overlay
