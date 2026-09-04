from __future__ import annotations

import datetime as dt
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pyproj import Transformer

from ..clouddisc.altaz_grid import (
    ALT_AZ_GRID_ALT_BINS,
    ALT_AZ_GRID_ALT_MAX_DEG,
    ALT_AZ_GRID_ALT_MIN_DEG,
    ALT_AZ_GRID_AZ_BINS,
    ALT_AZ_GRID_AZ_MAX_DEG,
    ALT_AZ_GRID_AZ_MIN_DEG,
    CloudAltAzGrid,
)
from ..clouddisc.types import SourceKey
from ..paths import GEOSATELLITE_EQDC_LONLAT_FILE
from ..utils.geostationary_image_mapping import build_equidistant_conic_projection
from .cache import read_projection_cache, write_projection_cache

DEFAULT_OUTPUT_PATH = Path("latest_cloud.png")
DEFAULT_RENDER_RADIUS_PX = 255
EARTH_RADIUS_KM = 6371.0
DEFAULT_CLOUD_HEIGHT_KM = 5.0
DEFAULT_GRID_NPZ = Path(GEOSATELLITE_EQDC_LONLAT_FILE)
DEFAULT_MIN_ALT_DEG = 3.0
DEFAULT_PROJECTION_SAMPLE_STEP = 1
EUROPE_MIN_LAT_DEG = 32.0
EUROPE_MAX_LAT_DEG = 73.0
EUROPE_MIN_LON_DEG = -15.0
EUROPE_MAX_LON_DEG = 35.0


@dataclass(frozen=True, slots=True)
class ProjectionInverse:
    projection: object
    lonlat_to_proj: Transformer
    inverse_matrix: np.ndarray
    offset: np.ndarray

    def lonlat_to_pixel(self, lon_deg: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        x_proj, y_proj = self.lonlat_to_proj.transform(lon_deg, lat_deg)
        proj = np.stack(
            [
                np.asarray(x_proj, dtype=np.float64) - float(self.offset[0]),
                np.asarray(y_proj, dtype=np.float64) - float(self.offset[1]),
            ],
            axis=0,
        )
        pixel = self.inverse_matrix @ proj.reshape(2, -1)
        x_px = pixel[0].reshape(np.asarray(lon_deg).shape)
        y_px = pixel[1].reshape(np.asarray(lat_deg).shape)
        return x_px, y_px


def _load_projection_inverse(grid_npz: Path) -> ProjectionInverse:
    with np.load(grid_npz, allow_pickle=False) as data:
        required = {
            "longitude_of_projection_origin",
            "latitude_of_projection_origin",
            "standard_parallel_1",
            "standard_parallel_2",
            "pixel_to_proj_x_coeffs",
            "pixel_to_proj_y_coeffs",
        }
        missing = sorted(name for name in required if name not in data)
        if missing:
            raise ValueError(f"{grid_npz} is missing required fields: {', '.join(missing)}")
        projection = build_equidistant_conic_projection(
            longitude_of_projection_origin=float(np.asarray(data["longitude_of_projection_origin"], dtype=np.float64).item()),
            latitude_of_projection_origin=float(np.asarray(data["latitude_of_projection_origin"], dtype=np.float64).item()),
            standard_parallel_1=float(np.asarray(data["standard_parallel_1"], dtype=np.float64).item()),
            standard_parallel_2=float(np.asarray(data["standard_parallel_2"], dtype=np.float64).item()),
        )
        coeffs_x = np.asarray(data["pixel_to_proj_x_coeffs"], dtype=np.float64)
        coeffs_y = np.asarray(data["pixel_to_proj_y_coeffs"], dtype=np.float64)
        if coeffs_x.shape != (3,) or coeffs_y.shape != (3,):
            raise ValueError("Pixel-to-projection coefficients must be length-3 vectors.")
        matrix = np.array(
            [
                [float(coeffs_x[0]), float(coeffs_x[1])],
                [float(coeffs_y[0]), float(coeffs_y[1])],
            ],
            dtype=np.float64,
        )
        det = float(np.linalg.det(matrix))
        if abs(det) < 1e-12:
            raise ValueError("Projection matrix is singular and cannot be inverted.")
        inverse = np.linalg.inv(matrix)
        offset = np.array([float(coeffs_x[2]), float(coeffs_y[2])], dtype=np.float64)
    return ProjectionInverse(
        projection=projection,
        lonlat_to_proj=Transformer.from_crs("EPSG:4326", projection, always_xy=True),
        inverse_matrix=inverse,
        offset=offset,
    )


def _sample_bilinear(gray: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    h, w = gray.shape
    x_safe = np.where(np.isfinite(x), x, 0.0)
    y_safe = np.where(np.isfinite(y), y, 0.0)
    x0 = np.floor(x_safe).astype(np.int64)
    y0 = np.floor(y_safe).astype(np.int64)
    x1 = x0 + 1
    y1 = y0 + 1
    valid = (x0 >= 0) & (y0 >= 0) & (x1 < w) & (y1 < h)
    out = np.zeros(x.shape, dtype=np.float32)
    if not np.any(valid):
        return out
    xv0 = x0[valid]
    yv0 = y0[valid]
    xv1 = x1[valid]
    yv1 = y1[valid]
    dx = (x[valid] - xv0).astype(np.float32)
    dy = (y[valid] - yv0).astype(np.float32)
    v00 = gray[yv0, xv0]
    v10 = gray[yv0, xv1]
    v01 = gray[yv1, xv0]
    v11 = gray[yv1, xv1]
    out_valid = (
        v00 * (1.0 - dx) * (1.0 - dy)
        + v10 * dx * (1.0 - dy)
        + v01 * (1.0 - dx) * dy
        + v11 * dx * dy
    )
    out[valid] = out_valid.astype(np.float32, copy=False)
    return out


def _altaz_grid_to_source_pixel_coords(
    alt_grid: np.ndarray,
    az_grid: np.ndarray,
    *,
    observer_lat: float,
    observer_lon: float,
    grid_npz: Path,
    cloud_height_km: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map an observer-centric alt/az grid to source-image pixel coordinates."""
    projection_inverse = _load_projection_inverse(grid_npz.expanduser())

    observer_pos_ecef = np.asarray(
        [
            float(EARTH_RADIUS_KM) * math.cos(math.radians(observer_lat)) * math.cos(math.radians(observer_lon)),
            float(EARTH_RADIUS_KM) * math.cos(math.radians(observer_lat)) * math.sin(math.radians(observer_lon)),
            float(EARTH_RADIUS_KM) * math.sin(math.radians(observer_lat)),
        ],
        dtype=np.float64,
    )
    up_vec = np.asarray(
        [
            math.cos(math.radians(observer_lat)) * math.cos(math.radians(observer_lon)),
            math.cos(math.radians(observer_lat)) * math.sin(math.radians(observer_lon)),
            math.sin(math.radians(observer_lat)),
        ],
        dtype=np.float64,
    )
    east = np.asarray([-math.sin(math.radians(observer_lon)), math.cos(math.radians(observer_lon)), 0.0], dtype=np.float64)
    north = np.asarray(
        [
            -math.sin(math.radians(observer_lat)) * math.cos(math.radians(observer_lon)),
            -math.sin(math.radians(observer_lat)) * math.sin(math.radians(observer_lon)),
            math.cos(math.radians(observer_lat)),
        ],
        dtype=np.float64,
    )

    alt_rad = np.radians(np.asarray(alt_grid, dtype=np.float64))
    az_rad = np.radians(np.asarray(az_grid, dtype=np.float64))
    d = (
        np.sin(alt_rad)[..., None] * up_vec
        + np.cos(alt_rad)[..., None]
        * (np.sin(az_rad)[..., None] * east + np.cos(az_rad)[..., None] * north)
    )
    d /= np.linalg.norm(d, axis=2, keepdims=True) + 1e-12

    cloud_shell_km = float(EARTH_RADIUS_KM) + float(cloud_height_km)
    b_quad = 2.0 * np.sum(observer_pos_ecef * d, axis=2)
    c_quad = float(np.dot(observer_pos_ecef, observer_pos_ecef)) - float(cloud_shell_km * cloud_shell_km)
    discriminant = b_quad * b_quad - 4.0 * c_quad
    valid_intersection = discriminant >= 0
    sqrt_disc = np.sqrt(np.maximum(discriminant, 0.0))
    t1 = (-b_quad - sqrt_disc) / 2.0
    t2 = (-b_quad + sqrt_disc) / 2.0
    t = np.where(t1 > 1e-6, t1, np.where(t2 > 1e-6, t2, np.nan))

    points = observer_pos_ecef + d * t[..., None]
    x_int = points[..., 0]
    y_int = points[..., 1]
    z_int = points[..., 2]
    lon_grid = np.full(alt_grid.shape, np.nan, dtype=np.float32)
    lat_grid = np.full(alt_grid.shape, np.nan, dtype=np.float32)
    lon_grid[valid_intersection] = np.degrees(np.arctan2(y_int, x_int))[valid_intersection]
    hyp = np.hypot(x_int, y_int)
    lat_grid[valid_intersection] = np.degrees(np.arctan2(z_int, hyp))[valid_intersection]

    x_src, y_src = projection_inverse.lonlat_to_pixel(lon_grid, lat_grid)
    valid = valid_intersection & np.isfinite(x_src) & np.isfinite(y_src)
    return (
        np.asarray(x_src, dtype=np.float32),
        np.asarray(y_src, dtype=np.float32),
        np.asarray(valid, dtype=bool),
    )


def build_altaz_grid_from_source_gray(
    source_gray: np.ndarray,
    *,
    observer_lat: float,
    observer_lon: float,
    kind: str,
    time_utc: dt.datetime,
    grid_npz: Path = DEFAULT_GRID_NPZ,
    cloud_height_km: float = DEFAULT_CLOUD_HEIGHT_KM,
    alt_bins: int = ALT_AZ_GRID_ALT_BINS,
    az_bins: int = ALT_AZ_GRID_AZ_BINS,
) -> CloudAltAzGrid:
    """Project a geo-satellite source image directly into the observer alt/az grid."""
    gray = np.asarray(source_gray, dtype=np.uint8)
    if gray.ndim != 2:
        raise ValueError("source_gray must have shape (H, W)")

    alt_edges = np.linspace(ALT_AZ_GRID_ALT_MIN_DEG, ALT_AZ_GRID_ALT_MAX_DEG, alt_bins + 1)
    az_edges = np.linspace(ALT_AZ_GRID_AZ_MIN_DEG, ALT_AZ_GRID_AZ_MAX_DEG, az_bins + 1)
    alt_centers = (alt_edges[:-1] + alt_edges[1:]) * 0.5
    az_centers = (az_edges[:-1] + az_edges[1:]) * 0.5
    alt_grid, az_grid = np.meshgrid(alt_centers, az_centers, indexing="ij")

    x_src, y_src, valid = _altaz_grid_to_source_pixel_coords(
        alt_grid,
        az_grid,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
        grid_npz=grid_npz,
        cloud_height_km=cloud_height_km,
    )
    sampled = _sample_bilinear(gray, x_src, y_src)

    amount = np.zeros((alt_bins, az_bins), dtype=np.float32)
    if np.any(valid):
        amount[valid] = sampled[valid].astype(np.float32, copy=False) / 255.0

    missing_mask = np.zeros((alt_bins, az_bins), dtype=np.uint8)
    missing_mask[~valid] = 255

    coverage_ratio = float(np.count_nonzero(valid)) / float(valid.size) if valid.size > 0 else 1.0
    source_key = SourceKey(
        satellite="Geo-sat",
        provider=str(kind),
        timeslot_utc=time_utc,
    )
    grid_resolution_deg = min(
        (ALT_AZ_GRID_AZ_MAX_DEG - ALT_AZ_GRID_AZ_MIN_DEG) / max(1, az_bins),
        (ALT_AZ_GRID_ALT_MAX_DEG - ALT_AZ_GRID_ALT_MIN_DEG) / max(1, alt_bins),
    )
    return CloudAltAzGrid(
        amount=amount,
        missing_mask=missing_mask,
        alt_min_deg=ALT_AZ_GRID_ALT_MIN_DEG,
        alt_max_deg=ALT_AZ_GRID_ALT_MAX_DEG,
        az_min_deg=ALT_AZ_GRID_AZ_MIN_DEG,
        az_max_deg=ALT_AZ_GRID_AZ_MAX_DEG,
        observer_lat=float(observer_lat),
        observer_lon=float(observer_lon),
        satellite="Geo-sat",
        product=str(kind),
        time_utc=time_utc,
        shells_km=(),
        source_key=source_key,
        coverage_ratio=coverage_ratio,
        source_completeness_ratio=None,
        grid_resolution_deg=float(grid_resolution_deg),
    )


def _projection_cache_key(
    *,
    source_shape: tuple[int, int],
    observer_lat: float,
    observer_lon: float,
    alt: float,
    az: float,
    fov_deg: float,
    cloud_height_km: float,
    radius_px: int,
    sample_step: int,
    grid_npz: Path,
) -> str:
    grid_path = grid_npz.expanduser().resolve()
    try:
        grid_stat = grid_path.stat()
        grid_mtime_ns = int(grid_stat.st_mtime_ns)
        grid_size = int(grid_stat.st_size)
    except OSError:
        grid_mtime_ns = -1
        grid_size = -1
    payload = {
        "alt": round(float(alt), 6),
        "az": round(float(az), 6),
        "cloud_height_km": round(float(cloud_height_km), 6),
        "edge_fov_deg": round(float(fov_deg), 6),
        "grid_mtime_ns": grid_mtime_ns,
        "grid_path": str(grid_path),
        "grid_size": grid_size,
        "observer_lat": round(float(observer_lat), 6),
        "observer_lon": round(float(observer_lon), 6),
        "radius_px": int(radius_px),
        "sample_step": int(sample_step),
        "source_shape": [int(source_shape[0]), int(source_shape[1])],
    }
    return json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _build_projection_sample(
    *,
    observer_lat: float,
    observer_lon: float,
    alt: float,
    az: float,
    fov_deg: float,
    grid_npz: Path,
    cloud_height_km: float,
    radius_px: int,
    sample_step: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    projection_inverse = _load_projection_inverse(grid_npz.expanduser())
    image_size = 2 * int(radius_px)
    sample_step = int(sample_step)
    if sample_step < 1:
        raise ValueError("sample_step must be >= 1")
    coarse_size = int(math.ceil(image_size / float(sample_step)))
    pixel_coords = ((np.arange(coarse_size, dtype=np.float32) + 0.5) * float(sample_step)) - float(radius_px)
    x, y = np.meshgrid(pixel_coords, -pixel_coords)
    cloud_shell_km = float(EARTH_RADIUS_KM) + float(cloud_height_km)
    fov_limit = float(fov_deg) * math.sqrt(2.0)
    if float(fov_deg) > fov_limit + 1e-6:
        raise ValueError(f"mask_fov_deg ({fov_deg} deg) exceeds geometric limit ({fov_limit:.2f} deg)")
    if not (0 < float(fov_deg) <= 180.0):
        raise ValueError("edge_fov_deg must be in (0, 180]")

    observer_pos_ecef = np.asarray(
        [
            float(EARTH_RADIUS_KM) * math.cos(math.radians(observer_lat)) * math.cos(math.radians(observer_lon)),
            float(EARTH_RADIUS_KM) * math.cos(math.radians(observer_lat)) * math.sin(math.radians(observer_lon)),
            float(EARTH_RADIUS_KM) * math.sin(math.radians(observer_lat)),
        ],
        dtype=np.float64,
    )
    up_vec = np.asarray(
        [
            math.cos(math.radians(observer_lat)) * math.cos(math.radians(observer_lon)),
            math.cos(math.radians(observer_lat)) * math.sin(math.radians(observer_lon)),
            math.sin(math.radians(observer_lat)),
        ],
        dtype=np.float64,
    )
    east = np.asarray([-math.sin(math.radians(observer_lon)), math.cos(math.radians(observer_lon)), 0.0], dtype=np.float64)
    north = np.asarray(
        [
            -math.sin(math.radians(observer_lat)) * math.cos(math.radians(observer_lon)),
            -math.sin(math.radians(observer_lat)) * math.sin(math.radians(observer_lon)),
            math.cos(math.radians(observer_lat)),
        ],
        dtype=np.float64,
    )
    az_rad = math.radians(az)
    alt_rad = math.radians(alt)
    center_view_dir = (
        math.sin(alt_rad) * up_vec
        + math.cos(alt_rad) * (math.sin(az_rad) * east + math.cos(az_rad) * north)
    )
    center_view_dir /= np.linalg.norm(center_view_dir) or 1.0
    ty_unnormalized = up_vec - np.dot(up_vec, center_view_dir) * center_view_dir
    ty = ty_unnormalized / (np.linalg.norm(ty_unnormalized) or 1.0)
    tx = np.cross(center_view_dir, ty)

    pixel_radius = np.hypot(x, y)
    rho_deg = (float(fov_deg) * pixel_radius / float(radius_px)).astype(np.float32)
    psi = np.arctan2(y, x)
    cos_rho = np.cos(np.deg2rad(rho_deg))
    sin_rho = np.sin(np.deg2rad(rho_deg))
    b = np.cos(psi)[..., None] * tx + np.sin(psi)[..., None] * ty
    d = cos_rho[..., None] * center_view_dir + sin_rho[..., None] * b
    d /= np.linalg.norm(d, axis=2, keepdims=True) + 1e-12
    alt_vis = np.degrees(np.arcsin(np.dot(d, up_vec)))
    visible_mask = alt_vis >= float(DEFAULT_MIN_ALT_DEG)
    fov_mask = rho_deg <= float(fov_deg) + 1e-6
    inside = fov_mask & visible_mask

    b_quad = 2.0 * np.sum(observer_pos_ecef * d, axis=2)
    c_quad = float(np.dot(observer_pos_ecef, observer_pos_ecef)) - float(cloud_shell_km * cloud_shell_km)
    discriminant = b_quad * b_quad - 4.0 * c_quad
    valid_intersection = discriminant >= 0
    sqrt_disc = np.sqrt(np.maximum(discriminant, 0.0))
    t1 = (-b_quad - sqrt_disc) / 2.0
    t2 = (-b_quad + sqrt_disc) / 2.0
    t = np.where(t1 > 1e-6, t1, np.where(t2 > 1e-6, t2, np.nan))
    P = observer_pos_ecef + d * t[..., None]
    x_int, y_int, z_int = P[..., 0], P[..., 1], P[..., 2]
    lon_grid = np.full((coarse_size, coarse_size), np.nan, dtype=np.float32)
    lat_grid = np.full((coarse_size, coarse_size), np.nan, dtype=np.float32)
    lon_grid[valid_intersection] = np.degrees(np.arctan2(y_int, x_int))[valid_intersection]
    hyp = np.hypot(x_int, y_int)
    lat_grid[valid_intersection] = np.degrees(np.arctan2(z_int, hyp))[valid_intersection]
    lon_grid[~inside] = np.nan
    lat_grid[~inside] = np.nan
    x_src, y_src = projection_inverse.lonlat_to_pixel(lon_grid, lat_grid)
    valid = inside & np.isfinite(x_src) & np.isfinite(y_src)
    return (
        np.asarray(x_src, dtype=np.float32),
        np.asarray(y_src, dtype=np.float32),
        np.asarray(valid, dtype=bool),
    )


def project_gray_image_to_disc(
    source_gray: np.ndarray,
    *,
    observer_lat: float,
    observer_lon: float,
    alt: float,
    az: float,
    fov_deg: float,
    grid_npz: Path = DEFAULT_GRID_NPZ,
    cloud_height_km: float = DEFAULT_CLOUD_HEIGHT_KM,
    radius_px: int = DEFAULT_RENDER_RADIUS_PX,
) -> np.ndarray:
    image_size = 2 * int(radius_px)
    source_shape = tuple(int(dim) for dim in np.asarray(source_gray).shape[:2])
    cache_key = _projection_cache_key(
        source_shape=source_shape,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
        alt=alt,
        az=az,
        fov_deg=fov_deg,
        cloud_height_km=cloud_height_km,
        radius_px=radius_px,
        sample_step=DEFAULT_PROJECTION_SAMPLE_STEP,
        grid_npz=grid_npz,
    )
    cached = read_projection_cache(cache_key=cache_key, source_shape=source_shape)
    if cached is None:
        x_src, y_src, valid = _build_projection_sample(
            observer_lat=observer_lat,
            observer_lon=observer_lon,
            alt=alt,
            az=az,
            fov_deg=fov_deg,
            grid_npz=grid_npz,
            cloud_height_km=cloud_height_km,
            radius_px=radius_px,
            sample_step=DEFAULT_PROJECTION_SAMPLE_STEP,
        )
        write_projection_cache(
            cache_key=cache_key,
            source_shape=source_shape,
            x_src=x_src,
            y_src=y_src,
            valid_mask=valid,
        )
    else:
        x_src, y_src, valid = cached
    sampled = _sample_bilinear(np.asarray(source_gray), x_src, y_src)
    output = np.zeros((image_size, image_size), dtype=np.uint8)
    if DEFAULT_PROJECTION_SAMPLE_STEP > 1:
        sampled = np.repeat(np.repeat(sampled, DEFAULT_PROJECTION_SAMPLE_STEP, axis=0), DEFAULT_PROJECTION_SAMPLE_STEP, axis=1)[
            :image_size,
            :image_size,
        ]
        valid = np.repeat(np.repeat(valid, DEFAULT_PROJECTION_SAMPLE_STEP, axis=0), DEFAULT_PROJECTION_SAMPLE_STEP, axis=1)[
            :image_size,
            :image_size,
        ]
    output[valid] = np.clip(np.rint(sampled[valid]), 0.0, 255.0).astype(np.uint8, copy=False)
    return output


def render_gray_image_to_cloud_rgba(gray: np.ndarray) -> np.ndarray:
    """Convert a grayscale cloud amount image into a white RGBA cloud layer."""
    arr = np.asarray(gray, dtype=np.uint8)
    if arr.ndim != 2:
        raise ValueError("gray must have shape (H, W)")
    rgba = np.zeros((arr.shape[0], arr.shape[1], 4), dtype=np.uint8)
    rgba[..., :3] = 255
    rgba[..., 3] = arr
    return rgba
