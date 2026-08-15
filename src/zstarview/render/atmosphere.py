"""Spherical-Earth single-scattering model for the daytime sky."""

from __future__ import annotations

import math

import numpy as np

EARTH_RADIUS_KM = 6371.0
ATMOSPHERE_TOP_KM = 100.0
RAYLEIGH_SCALE_HEIGHT_KM = 8.0
AEROSOL_SCALE_HEIGHT_KM = 1.4
OZONE_SHELL_BOTTOM_KM = 15.0
OZONE_SHELL_TOP_KM = 35.0
AEROSOL_DENSITY_RATIO = 0.025
AEROSOL_REFERENCE_AOD550 = 0.15
REPRESENTATIVE_WAVELENGTHS_NM = np.array([650.0, 550.0, 450.0], dtype=np.float32)
RAYLEIGH_WAVELENGTH_EXPONENT = 4.0
AEROSOL_ANGSTROM_EXPONENT = 0.7
# Effective Chappuis-band optical depths for the representative RGB bands.
# These are normalized to a vertical crossing of the ozone shell and are a
# deliberately compact approximation of the wavelength-resolved absorption.
# Keep the absorption path available for later tuning with a deliberately
# weak contribution while comparing the twilight colour balance.
OZONE_EXTINCTION_RGB = np.array([0.0025, 0.005, 0.0001875], dtype=np.float32)

RAYLEIGH_SCATTERING_RGB = (
    REPRESENTATIVE_WAVELENGTHS_NM[-1] / REPRESENTATIVE_WAVELENGTHS_NM
) ** RAYLEIGH_WAVELENGTH_EXPONENT
RAYLEIGH_EXTINCTION_RGB = RAYLEIGH_SCATTERING_RGB.copy()
AEROSOL_SCATTERING_RGB = (
    REPRESENTATIVE_WAVELENGTHS_NM[1] / REPRESENTATIVE_WAVELENGTHS_NM
) ** AEROSOL_ANGSTROM_EXPONENT
AEROSOL_EXTINCTION_RGB = AEROSOL_SCATTERING_RGB.copy()
SUN_RADIANCE_RGB = np.ones(3, dtype=np.float32)
# Small empirical balance for the three-channel approximation.  The green
# channel otherwise remains too strong in low-sky sunset colours and makes
# orange light appear yellow.  Keep this separate from optical depths so the
# atmospheric path lengths remain physically interpretable.
SCATTERING_RGB_BALANCE = np.array([0.98, 1.0, 1.0], dtype=np.float32)
OPTICAL_DEPTH_SCALE = 0.018
AEROSOL_OPTICAL_DEPTH_SCALE = 0.018
MIE_ANISOTROPY = 0.76
DISPLAY_EXPOSURE = 2.8
# Effective blue radiance from higher-order twilight scattering. This is an
# RGB approximation, added before display conversion, rather than an emitted
# blue layer. It is intentionally strongest in the upper sky.
TWILIGHT_MULTIPLE_SCATTERING_RGB = np.array(
    [0.008, 0.020, 0.070], dtype=np.float32
)
TWILIGHT_MULTIPLE_SCATTERING_START_ALT_DEG = 3.0
TWILIGHT_MULTIPLE_SCATTERING_END_ALT_DEG = -12.0


def _ray_shell_distance(
    origin: np.ndarray,
    direction: np.ndarray,
    radius_km: float,
) -> np.ndarray:
    """Return the forward intersection distance with a spherical shell."""
    b = np.sum(origin * direction, axis=-1)
    c = np.sum(origin * origin, axis=-1) - radius_km * radius_km
    discriminant = np.maximum(0.0, b * b - c)
    return -b + np.sqrt(discriminant)


def _ray_hits_earth(
    origin: np.ndarray,
    direction: np.ndarray,
) -> np.ndarray:
    """Return whether a forward ray intersects the solid Earth."""
    b = np.sum(origin * direction, axis=-1)
    c = np.sum(origin * origin, axis=-1) - EARTH_RADIUS_KM * EARTH_RADIUS_KM
    return (b < 0.0) & (b * b >= c)


def _shell_path_length(
    origin: np.ndarray,
    direction: np.ndarray,
    max_distance: np.ndarray,
    inner_radius_km: float,
    outer_radius_km: float,
) -> np.ndarray:
    """Return the distance within a spherical shell along a finite ray."""
    b = np.sum(origin * direction, axis=-1)
    c = np.sum(origin * origin, axis=-1)

    def sphere_interval(radius_km: float) -> np.ndarray:
        discriminant = b * b - (c - radius_km * radius_km)
        valid = discriminant >= 0.0
        root = np.sqrt(np.maximum(0.0, discriminant))
        near = np.maximum(0.0, -b - root)
        far = np.minimum(max_distance, -b + root)
        return np.where(valid, np.maximum(0.0, far - near), 0.0)

    return sphere_interval(outer_radius_km) - sphere_interval(inner_radius_km)


def _twilight_multiple_scattering_radiance(
    view_alt_deg: np.ndarray,
    sun_alt_deg: float,
) -> np.ndarray:
    """Return a high-sky blue radiance approximation for civil/nautical twilight."""
    sun_t = np.clip(
        (float(sun_alt_deg) - TWILIGHT_MULTIPLE_SCATTERING_END_ALT_DEG)
        / (
            TWILIGHT_MULTIPLE_SCATTERING_START_ALT_DEG
            - TWILIGHT_MULTIPLE_SCATTERING_END_ALT_DEG
        ),
        0.0,
        1.0,
    )
    sun_weight = sun_t * sun_t * (3.0 - 2.0 * sun_t)
    view_t = np.clip(np.asarray(view_alt_deg, dtype=np.float32) / 90.0, 0.0, 1.0)
    view_weight = view_t * view_t * (3.0 - 2.0 * view_t)
    return (
        view_weight.reshape(-1, 1)
        * sun_weight
        * TWILIGHT_MULTIPLE_SCATTERING_RGB[None, :]
    )


def _height_density(position: np.ndarray, scale_height_km: float) -> np.ndarray:
    height = np.maximum(
        0.0,
        np.linalg.norm(position, axis=-1) - EARTH_RADIUS_KM,
    )
    return np.exp(-height / scale_height_km).astype(np.float32)


def _sun_is_visible(point: np.ndarray, sun_direction: np.ndarray) -> np.ndarray:
    """Return whether the path from each point toward the Sun clears Earth."""
    return ~_ray_hits_earth(point, sun_direction[None, :])


def _sun_column_densities(
    point: np.ndarray,
    sun_direction: np.ndarray,
    steps: int,
    atmosphere_top_km: float,
    aerosol_density_scale: float,
) -> tuple[np.ndarray, np.ndarray]:
    distance = _ray_shell_distance(
        point,
        sun_direction[None, :],
        EARTH_RADIUS_KM + atmosphere_top_km,
    )
    fraction = (np.arange(steps, dtype=np.float32) + 0.5) / steps
    samples = (
        point[:, None, :]
        + distance[:, None, None]
        * fraction[None, :, None]
        * sun_direction[None, None, :]
    )
    flat_samples = samples.reshape(-1, 3)
    shape = (point.shape[0], steps)
    rayleigh = _height_density(flat_samples, RAYLEIGH_SCALE_HEIGHT_KM).reshape(shape)
    aerosol = (
        AEROSOL_DENSITY_RATIO
        * aerosol_density_scale
        * _height_density(flat_samples, AEROSOL_SCALE_HEIGHT_KM).reshape(shape)
    )
    return (
        (np.mean(rayleigh, axis=1) * distance).astype(np.float32),
        (np.mean(aerosol, axis=1) * distance).astype(np.float32),
    )


def _henyey_greenstein(cos_theta: np.ndarray) -> np.ndarray:
    g = MIE_ANISOTROPY
    denominator = np.power(1.0 + g * g - 2.0 * g * cos_theta, 1.5)
    return ((1.0 - g * g) / (4.0 * math.pi * denominator)).astype(np.float32)


def _direction_vectors(
    altitude_deg: np.ndarray,
    azimuth_deg: np.ndarray,
) -> np.ndarray:
    altitude = np.radians(altitude_deg)
    azimuth = np.radians(azimuth_deg)
    return np.stack(
        [
            np.cos(altitude) * np.sin(azimuth),
            np.cos(altitude) * np.cos(azimuth),
            np.sin(altitude),
        ],
        axis=1,
    ).astype(np.float32)


def atmospheric_sky_samples(
    view_alt_deg: np.ndarray,
    view_az_deg: np.ndarray,
    sun_altaz: tuple[float, float],
    *,
    observer_height_m: float = 0.0,
    atmosphere_top_km: float = ATMOSPHERE_TOP_KM,
    view_steps: int = 32,
    sun_steps: int = 12,
    exposure: float = 1.0,
    aerosol_optical_depth: float = AEROSOL_REFERENCE_AOD550,
) -> np.ndarray:
    """Return RGB sky radiance for local altitude/azimuth directions.

    The result contains only sunlight scattered into the viewing direction.
    Direct solar radiance is intentionally not included.
    """
    view_alt = np.asarray(view_alt_deg, dtype=np.float32)
    view_az = np.asarray(view_az_deg, dtype=np.float32)
    if view_alt.shape != view_az.shape:
        raise ValueError("view_alt_deg and view_az_deg must have the same shape")
    if view_alt.size == 0:
        return np.zeros((0, 3), dtype=np.float32)
    if atmosphere_top_km <= 0.0 or view_steps < 1 or sun_steps < 1:
        raise ValueError("atmosphere_top_km, view_steps, and sun_steps must be positive")
    if not np.isfinite(aerosol_optical_depth) or aerosol_optical_depth < 0.0:
        raise ValueError("aerosol_optical_depth must be finite and non-negative")
    aerosol_density_scale = float(aerosol_optical_depth) / AEROSOL_REFERENCE_AOD550

    directions = _direction_vectors(view_alt.reshape(-1), view_az.reshape(-1))
    observer_radius = EARTH_RADIUS_KM + max(0.0, float(observer_height_m)) / 1000.0
    observer = np.array([0.0, 0.0, observer_radius], dtype=np.float32)
    sun_direction = _direction_vectors(
        np.asarray([sun_altaz[0]], dtype=np.float32),
        np.asarray([sun_altaz[1]], dtype=np.float32),
    )[0]

    max_distance = _ray_shell_distance(
        np.broadcast_to(observer, directions.shape),
        directions,
        EARTH_RADIUS_KM + atmosphere_top_km,
    )
    view_distance = max_distance / view_steps
    view_fraction = (np.arange(view_steps, dtype=np.float32) + 0.5) / view_steps
    result = np.zeros((directions.shape[0], 3), dtype=np.float32)

    for start in range(0, directions.shape[0], 4096):
        end = min(directions.shape[0], start + 4096)
        chunk_directions = directions[start:end]
        chunk_max_distance = max_distance[start:end]
        chunk_distance = view_distance[start:end]
        chunk_points = (
            observer[None, None, :]
            + chunk_max_distance[:, None, None]
            * view_fraction[None, :, None]
            * chunk_directions[:, None, :]
        )
        chunk_shape = (end - start, view_steps)
        flat_chunk_points = chunk_points.reshape(-1, 3)
        chunk_rayleigh = _height_density(
            flat_chunk_points,
            RAYLEIGH_SCALE_HEIGHT_KM,
        ).reshape(chunk_shape)
        chunk_aerosol = (
            AEROSOL_DENSITY_RATIO
            * aerosol_density_scale
            * _height_density(flat_chunk_points, AEROSOL_SCALE_HEIGHT_KM).reshape(
                chunk_shape
            )
        )
        rayleigh_view_column = np.cumsum(
            chunk_rayleigh * chunk_distance[:, None], axis=1
        ) - 0.5 * chunk_rayleigh * chunk_distance[:, None]
        aerosol_view_column = np.cumsum(
            chunk_aerosol * chunk_distance[:, None], axis=1
        ) - 0.5 * chunk_aerosol * chunk_distance[:, None]
        rayleigh_sun_column, aerosol_sun_column = _sun_column_densities(
            flat_chunk_points,
            sun_direction,
            sun_steps,
            atmosphere_top_km,
            aerosol_density_scale,
        )
        rayleigh_sun_column = rayleigh_sun_column.reshape(chunk_shape)
        aerosol_sun_column = aerosol_sun_column.reshape(chunk_shape)
        sun_visible = _sun_is_visible(flat_chunk_points, sun_direction).reshape(
            chunk_shape
        )
        view_sample_distance = chunk_distance[:, None] * view_fraction[None, :]
        view_directions = np.broadcast_to(
            chunk_directions[:, None, :], (end - start, view_steps, 3)
        )
        view_ozone_column = _shell_path_length(
            np.broadcast_to(observer, flat_chunk_points.shape),
            view_directions.reshape(-1, 3),
            view_sample_distance.reshape(-1),
            EARTH_RADIUS_KM + OZONE_SHELL_BOTTOM_KM,
            EARTH_RADIUS_KM + OZONE_SHELL_TOP_KM,
        ).reshape(chunk_shape)
        sun_ozone_path = _shell_path_length(
            flat_chunk_points,
            np.broadcast_to(sun_direction, flat_chunk_points.shape),
            _ray_shell_distance(
                flat_chunk_points,
                np.broadcast_to(sun_direction, flat_chunk_points.shape),
                EARTH_RADIUS_KM + atmosphere_top_km,
            ),
            EARTH_RADIUS_KM + OZONE_SHELL_BOTTOM_KM,
            EARTH_RADIUS_KM + OZONE_SHELL_TOP_KM,
        ).reshape(chunk_shape)
        ozone_column = (
            view_ozone_column + sun_ozone_path
        ) / (OZONE_SHELL_TOP_KM - OZONE_SHELL_BOTTOM_KM)
        cos_phase = np.sum(
            sun_direction[None, None, :] * directions[start:end, None, :],
            axis=-1,
        )
        rayleigh_phase = 3.0 / (16.0 * math.pi) * (1.0 + cos_phase**2)
        mie_phase = _henyey_greenstein(cos_phase)
        rayleigh_transmission = np.exp(
            -OPTICAL_DEPTH_SCALE
            * (rayleigh_sun_column[:, :, None] + rayleigh_view_column[:, :, None])
            * RAYLEIGH_EXTINCTION_RGB[None, None, :]
        )
        aerosol_transmission = np.exp(
            -AEROSOL_OPTICAL_DEPTH_SCALE
            * (aerosol_sun_column[:, :, None] + aerosol_view_column[:, :, None])
            * AEROSOL_EXTINCTION_RGB[None, None, :]
        )
        ozone_transmission = np.exp(
            -ozone_column[:, :, None] * OZONE_EXTINCTION_RGB[None, None, :]
        )
        transmission = rayleigh_transmission * aerosol_transmission * ozone_transmission
        scattering = (
            chunk_rayleigh[:, :, None]
            * RAYLEIGH_SCATTERING_RGB[None, None, :]
            * rayleigh_phase[:, :, None]
            + chunk_aerosol[:, :, None]
            * AEROSOL_SCATTERING_RGB[None, None, :]
            * mie_phase[:, :, None]
        )
        scattering *= SCATTERING_RGB_BALANCE[None, None, :]
        result[start:end] = np.sum(
            scattering
            * SUN_RADIANCE_RGB[None, None, :]
            * transmission
            * sun_visible[:, :, None]
            * chunk_distance[:, None, None],
            axis=1,
        )

    result += _twilight_multiple_scattering_radiance(
        view_alt.reshape(-1),
        float(sun_altaz[0]),
    )
    blocked = _ray_hits_earth(
        np.broadcast_to(observer, directions.shape),
        directions,
    )
    result[blocked] = 0.0
    result *= max(0.0, float(exposure))
    return np.clip(1.0 - np.exp(-result * DISPLAY_EXPOSURE), 0.0, 1.0).astype(
        np.float32
    )


def direct_solar_transmission_rgb(
    sun_alt_deg: float,
    *,
    observer_height_m: float = 0.0,
    atmosphere_top_km: float = ATMOSPHERE_TOP_KM,
    sun_steps: int = 32,
    aerosol_optical_depth: float = AEROSOL_REFERENCE_AOD550,
) -> np.ndarray:
    """Return RGB transmittance for direct sunlight reaching the observer.

    This is the extinction half of the sky model. It deliberately excludes
    in-scattered radiance because a solar-disc image is direct sunlight and
    must be colorized multiplicatively rather than with an additive sky glow.
    """
    sun_alt = float(sun_alt_deg)
    height = float(observer_height_m)
    top_km = float(atmosphere_top_km)
    steps = int(sun_steps)
    aod = float(aerosol_optical_depth)
    if not np.isfinite(sun_alt) or not np.isfinite(height) or height < 0.0:
        raise ValueError("sun_alt_deg and observer_height_m must be finite")
    if top_km <= 0.0 or steps < 1:
        raise ValueError("atmosphere_top_km must be positive and sun_steps >= 1")
    if not np.isfinite(aod) or aod < 0.0:
        raise ValueError("aerosol_optical_depth must be finite and non-negative")
    if sun_alt < 0.0:
        return np.zeros(3, dtype=np.float32)

    aerosol_density_scale = aod / AEROSOL_REFERENCE_AOD550
    observer_radius = EARTH_RADIUS_KM + height / 1000.0
    observer = np.array([[0.0, 0.0, observer_radius]], dtype=np.float32)
    sun_direction = _direction_vectors(
        np.asarray([sun_alt], dtype=np.float32),
        np.asarray([0.0], dtype=np.float32),
    )[0]
    if bool(_ray_hits_earth(observer, sun_direction[None, :])[0]):
        return np.zeros(3, dtype=np.float32)

    rayleigh_column, aerosol_column = _sun_column_densities(
        observer,
        sun_direction,
        steps,
        top_km,
        aerosol_density_scale,
    )
    distance = _ray_shell_distance(
        observer,
        sun_direction[None, :],
        EARTH_RADIUS_KM + top_km,
    )
    ozone_column = _shell_path_length(
        observer,
        sun_direction[None, :],
        distance,
        EARTH_RADIUS_KM + OZONE_SHELL_BOTTOM_KM,
        EARTH_RADIUS_KM + OZONE_SHELL_TOP_KM,
    ) / (OZONE_SHELL_TOP_KM - OZONE_SHELL_BOTTOM_KM)
    rayleigh_tau = (
        OPTICAL_DEPTH_SCALE * rayleigh_column[0] * RAYLEIGH_EXTINCTION_RGB
    )
    aerosol_tau = (
        AEROSOL_OPTICAL_DEPTH_SCALE * aerosol_column[0] * AEROSOL_EXTINCTION_RGB
    )
    ozone_tau = ozone_column[0] * OZONE_EXTINCTION_RGB
    return np.exp(-(rayleigh_tau + aerosol_tau + ozone_tau)).astype(np.float32)
