from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pyproj import Transformer

from ..clouddisc.projectors.az import az_project_lonlat_grid
from ..paths import GEOSATELLITE_EQDC_LONLAT_FILE
from ..utils.geostationary_image_mapping import build_equidistant_conic_projection

DEFAULT_OUTPUT_PATH = Path("latest_cloud.png")
DEFAULT_RENDER_RADIUS_PX = 255
EARTH_RADIUS_KM = 6371.0
DEFAULT_CLOUD_HEIGHT_KM = 5.0
DEFAULT_GRID_NPZ = Path(GEOSATELLITE_EQDC_LONLAT_FILE)
DEFAULT_MIN_ALT_DEG = 3.0
MAX_FOV_DEG = 135.0
EUROPE_MIN_LAT_DEG = 32.0
EUROPE_MAX_LAT_DEG = 73.0
EUROPE_MIN_LON_DEG = -15.0
EUROPE_MAX_LON_DEG = 35.0


def parse_observer_spec(value: str) -> tuple[float, float]:
    text = (value or "").strip()
    if not text.startswith("@"):
        raise argparse.ArgumentTypeError("Observer must be given as '@lat,lon'.")
    payload = text[1:].strip()
    parts = [part.strip() for part in payload.split(",")]
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise argparse.ArgumentTypeError("Observer must be given as '@lat,lon'.")
    try:
        lat = float(parts[0])
        lon = float(parts[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Observer must be given as '@lat,lon'.") from exc
    if not (-90.0 <= lat <= 90.0):
        raise argparse.ArgumentTypeError("Latitude must be between -90 and 90 degrees.")
    if not (-180.0 <= lon <= 180.0):
        raise argparse.ArgumentTypeError("Longitude must be between -180 and 180 degrees.")
    return float(lat), float(lon)


def parse_fov_deg(value: str) -> float:
    try:
        fov = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid FOV: {value!r}") from exc
    if not (0.0 < fov <= MAX_FOV_DEG):
        raise argparse.ArgumentTypeError(
            f"FOV must be greater than 0 and at most {MAX_FOV_DEG:.0f} degrees."
        )
    return float(fov)


def validate_observer_bounds(lat: float, lon: float) -> None:
    if not (EUROPE_MIN_LAT_DEG <= lat <= EUROPE_MAX_LAT_DEG):
        raise argparse.ArgumentTypeError(
            f"Latitude must be between {EUROPE_MIN_LAT_DEG:.0f} and {EUROPE_MAX_LAT_DEG:.0f} degrees for the Europe workflow."
        )
    if not (EUROPE_MIN_LON_DEG <= lon <= EUROPE_MAX_LON_DEG):
        raise argparse.ArgumentTypeError(
            f"Longitude must be between {EUROPE_MIN_LON_DEG:.0f} and {EUROPE_MAX_LON_DEG:.0f} degrees for the Europe workflow."
        )


@dataclass(frozen=True, slots=True)
class ProjectionInverse:
    projection: object
    inverse_matrix: np.ndarray
    offset: np.ndarray

    def lonlat_to_pixel(self, lon_deg: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        to_proj = Transformer.from_crs("EPSG:4326", self.projection, always_xy=True)
        x_proj, y_proj = to_proj.transform(lon_deg, lat_deg)
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
    return ProjectionInverse(projection=projection, inverse_matrix=inverse, offset=offset)


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
    projection_inverse = _load_projection_inverse(grid_npz.expanduser())
    image_size = 2 * int(radius_px)
    cloud_shell_km = float(EARTH_RADIUS_KM) + float(cloud_height_km)
    lon_grid, lat_grid, inside = az_project_lonlat_grid(
        observer_lat,
        observer_lon,
        alt,
        az,
        radius_px,
        float(cloud_shell_km),
        alt_min_deg=float(DEFAULT_MIN_ALT_DEG),
        mask_fov_deg=float(fov_deg),
        edge_fov_deg=float(fov_deg),
    )
    x_src, y_src = projection_inverse.lonlat_to_pixel(lon_grid, lat_grid)
    sampled = _sample_bilinear(source_gray, x_src, y_src)
    output = np.zeros((image_size, image_size), dtype=np.uint8)
    valid = inside & np.isfinite(x_src) & np.isfinite(y_src)
    output[valid] = np.clip(np.rint(sampled[valid]), 0.0, 255.0).astype(np.uint8, copy=False)
    return output


def render_gray_image_to_rgba(gray: np.ndarray) -> np.ndarray:
    arr = np.asarray(gray, dtype=np.uint8)
    if arr.ndim != 2:
        raise ValueError("gray must have shape (H, W)")
    rgba = np.zeros((arr.shape[0], arr.shape[1], 4), dtype=np.uint8)
    rgba[..., :3] = arr[..., None]
    rgba[..., 3] = 255
    return rgba


def render_gray_image_to_cloud_rgba(gray: np.ndarray) -> np.ndarray:
    """Convert a grayscale cloud amount image into a white RGBA cloud layer."""
    arr = np.asarray(gray, dtype=np.uint8)
    if arr.ndim != 2:
        raise ValueError("gray must have shape (H, W)")
    rgba = np.zeros((arr.shape[0], arr.shape[1], 4), dtype=np.uint8)
    rgba[..., :3] = 255
    rgba[..., 3] = arr
    return rgba
