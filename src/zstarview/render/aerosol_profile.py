"""Access the bundled CAMS aerosol optical-depth climatology."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import as_file, files

import numpy as np

logger = logging.getLogger(__name__)

ASSET_PATH = "aerosol/cams_aod550_climatology.npz"
DEFAULT_AOD550 = 0.12


@dataclass(frozen=True)
class AerosolClimatology:
    """Validated monthly AOD550 values on a regular latitude/longitude grid."""

    aod550: np.ndarray
    latitude: np.ndarray
    longitude: np.ndarray

    def value_at(self, latitude: float, longitude: float, month: int) -> float:
        if month < 1 or month > 12:
            raise ValueError("month must be between 1 and 12")
        if not np.isfinite(latitude) or not np.isfinite(longitude):
            raise ValueError("latitude and longitude must be finite")

        lat = float(np.clip(latitude, self.latitude.min(), self.latitude.max()))
        lon = float(longitude) % 360.0
        lat_index = int(np.abs(self.latitude - lat).argmin())
        lon_index = int(np.abs(self.longitude - lon).argmin())
        value = float(self.aod550[month - 1, lat_index, lon_index])
        if not np.isfinite(value) or value < 0.0:
            raise ValueError("bundled AOD value is invalid")
        return value


def _validate_asset(archive: np.lib.npyio.NpzFile) -> AerosolClimatology:
    required = {"aod550", "latitude", "longitude"}
    if not required.issubset(archive.files):
        raise ValueError("bundled AOD asset is missing required arrays")

    aod550 = np.asarray(archive["aod550"], dtype=np.float32)
    latitude = np.asarray(archive["latitude"], dtype=np.float32)
    longitude = np.asarray(archive["longitude"], dtype=np.float32)
    if aod550.ndim != 3 or aod550.shape[0] != 12:
        raise ValueError("bundled AOD asset must have twelve monthly layers")
    if latitude.ndim != 1 or longitude.ndim != 1:
        raise ValueError("bundled AOD coordinates must be one-dimensional")
    if aod550.shape[1:] != (latitude.size, longitude.size):
        raise ValueError("bundled AOD dimensions do not match coordinates")
    if latitude.size < 2 or longitude.size < 2:
        raise ValueError("bundled AOD grid is too small")
    if not np.all(np.isfinite(latitude)) or not np.all(np.isfinite(longitude)):
        raise ValueError("bundled AOD coordinates contain non-finite values")
    latitude_steps = np.diff(latitude)
    longitude_steps = np.diff(longitude)
    if not (
        np.all(latitude_steps > 0.0) or np.all(latitude_steps < 0.0)
    ) or not np.all(longitude_steps > 0.0):
        raise ValueError("bundled AOD coordinates must be monotonic")
    if not np.all(np.isfinite(aod550)) or np.any(aod550 < 0.0):
        raise ValueError("bundled AOD values must be finite and non-negative")
    return AerosolClimatology(aod550, latitude, longitude)


@lru_cache(maxsize=1)
def load_bundled_climatology() -> AerosolClimatology:
    """Load and validate the packaged CAMS asset once per process."""
    resource = files("zstarview.data").joinpath(ASSET_PATH)
    with as_file(resource) as asset_path, np.load(
        asset_path, allow_pickle=False
    ) as archive:
        return _validate_asset(archive)


_asset_warning_logged = False


def bundled_aod550_or_default(latitude: float, longitude: float, month: int) -> float:
    """Return the bundled AOD550 or the defensive fixed default."""
    global _asset_warning_logged
    try:
        return load_bundled_climatology().value_at(latitude, longitude, month)
    except (OSError, ValueError, KeyError, EOFError) as exc:
        if not _asset_warning_logged:
            logger.warning("Bundled CAMS AOD asset is unavailable: %s", exc)
            _asset_warning_logged = True
        return DEFAULT_AOD550
