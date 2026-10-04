"""Shade clouds represented by native satellite pixel voxel columns."""

from __future__ import annotations

import numpy as np
from numba import njit
from pyproj import Transformer

from ..clouddisc.altaz_grid import _blend_cloud_shell_weights
from ..clouddisc.render.grayscale import _bt_to_weight, _suppress_low_cloud_weight
from ..clouddisc.sampling.estimate_bt_warm_cold import (
    estimate_bt_cold_hybrid,
    estimate_bt_warm_from_equator_band,
    estimate_bt_warm_hybrid,
)


ENVIRONMENT_LIGHT_FRACTION = 0.05

SUNLIGHT_LEVELS = np.array([0.0, 0.01, 0.03, 0.10, 0.30, 1.0])
CLOUD_WHITENESS = np.array([0.18, 0.60, 0.85, 0.95, 0.99, 1.0])


@njit
def _segments(origin, direction, shape):
    """Exact intersections with axis-aligned voxel faces; distances are km."""
    result = []
    near, far = 0.0, 1000.0
    for axis in range(3):
        if abs(direction[axis]) < 1e-12:
            if origin[axis] < 0 or origin[axis] >= shape[axis]:
                return result
        else:
            a = -origin[axis] / direction[axis]
            b = (shape[axis] - origin[axis]) / direction[axis]
            near = max(near, min(a, b))
            far = min(far, max(a, b))
    t = near
    while t < far - 1e-8:
        p = origin + direction * (t + 1e-7)
        cell = np.floor(p).astype(np.int64)
        end = far
        for axis in range(3):
            if direction[axis] > 1e-12:
                end = min(end, (cell[axis] + 1 - origin[axis]) / direction[axis])
            elif direction[axis] < -1e-12:
                end = min(end, (cell[axis] - origin[axis]) / direction[axis])
        if end <= t:
            break
        if np.all(cell >= 0) and np.all(cell < shape):
            result.append((cell[0], cell[1], cell[2], end - t, t))
        t = end
    return result


@njit
def _light_with_environment(clear_sunlight: float, optical_depth: float) -> float:
    """Closed form of repeated L = L*T + environment*(1-T) steps."""
    transmission = np.exp(-optical_depth)
    environment = clear_sunlight * ENVIRONMENT_LIGHT_FRACTION
    return clear_sunlight * transmission + environment * (1.0 - transmission)


@njit
def _render(density, origin, rays, sun, base, opacity, show_grid):
    shape = np.array(density.shape)
    lights = np.full(density.shape, -1.0)
    output = base.copy()
    transmissions = np.ones(len(rays), dtype=np.float32)
    for i in range(len(rays)):
        transmission = 1.0
        color = np.zeros(3)
        for x, y, z, distance, entry in _segments(origin, rays[i], shape):
            amount = density[x, y, z]
            if amount <= 0:
                continue
            if lights[x, y, z] < 0:
                tau = 0.0
                if sun[2] > 0:
                    center = np.array([x + 0.5, y + 0.5, z + 0.5])
                    for sx, sy, sz, length, _ in _segments(center, sun, shape):
                        tau += density[sx, sy, sz] * length * 1.8
                    # Solar altitude sets clear-sky strength; clouds attenuate it.
                    sunlight = _light_with_environment(sun[2], tau)
                    lights[x, y, z] = np.interp(
                        sunlight, SUNLIGHT_LEVELS, CLOUD_WHITENESS
                    )
                else:
                    lights[x, y, z] = 0.18
            alpha = 1.0 - np.exp(-2.4 * amount * distance * opacity)
            # Cloud amount controls extinction, not a separate color multiplier.
            # Thin cloud keeps a pale color while more background shines through.
            brightness = lights[x, y, z]
            # Diagnostic dark seams on voxel faces, only where cloud exists.
            if show_grid:
                entry_point = origin + rays[i] * entry
                fractions = entry_point - np.floor(entry_point)
                face_distances = np.sort(np.minimum(fractions, 1.0 - fractions))
                if face_distances[1] < 0.035:
                    alpha = 1.0 - (1.0 - alpha) * 0.92
                    brightness *= 0.75
            color += transmission * alpha * brightness * np.array([0.96, 0.975, 1.0])
            transmission *= 1.0 - alpha
            if transmission < 1e-5:
                break
        output[i] = color + transmission * base[i]
        transmissions[i] = transmission
    return output, transmissions


def shade_native_voxels(
    source,
    lat,
    lon,
    alt,
    az,
    sun_alt,
    sun_az,
    base,
    *,
    opacity=0.85,
    show_grid=False,
    return_transmission=False,
):
    """Use raw pixels without interpolation; affine native footprints locally.

    Native geostationary projection is linearized at the observer. The resulting
    parallelogram columns retain pixel indices and spacing, with nine 1-km slabs.
    Curvature, parallax and B16 redistribution are omitted in this experiment.
    """
    area = source.data_array.attrs["area"]
    xmin, ymin, xmax, ymax = area.area_extent
    h, w = area.shape
    dx, dy = (xmax - xmin) / w, (ymax - ymin) / h
    transformer = Transformer.from_crs("EPSG:4326", area.crs, always_xy=True)
    x, y = transformer.transform(lon, lat)
    # One-km east/north offsets determine the local satellite projection metric.
    dlon = np.degrees(1.0 / (6371.0 * np.cos(np.radians(lat))))
    dlat = np.degrees(1.0 / 6371.0)
    xe, ye = transformer.transform(lon + dlon, lat)
    xn, yn = transformer.transform(lon, lat + dlat)
    metric = np.array(
        [[(xe - x) / dx, (xn - x) / dx], [-(ye - y) / dy, -(yn - y) / dy]]
    )
    pixel = np.array([(x - xmin) / dx, (ymax - y) / dy])
    if not np.all(np.isfinite(metric)):
        raise ValueError("observer is outside the satellite projection")
    # Limit the tangent volume to 200 km in each direction.
    radius = np.ceil(np.sum(np.abs(metric), axis=1) * 200).astype(int)
    lo = np.maximum(0, np.floor(pixel).astype(int) - radius)
    hi = np.minimum([w, h], np.floor(pixel).astype(int) + radius + 1)
    raw = np.asarray(
        source.data_array.isel(y=slice(lo[1], hi[1]), x=slice(lo[0], hi[0]))
        .compute()
        .values,
        dtype=np.float32,
    )
    _, eq = estimate_bt_warm_from_equator_band(
        source.data_array,
        lon_center_deg=lon,
        delta_lon=60,
        equator_lat=0,
        warm_p=97,
        half=5,
        equator_lat_half_band_deg=5,
    )
    valid = np.isfinite(raw)
    warm = estimate_bt_warm_hybrid(raw, valid, eq, fallback_bt_warm=310)
    cold = estimate_bt_cold_hybrid(raw, valid, eq, warm)
    amount = _suppress_low_cloud_weight(_bt_to_weight(raw, warm, cold))
    weights = np.repeat(_blend_cloud_shell_weights(float(np.mean(amount))), 3) / 3
    density = np.ascontiguousarray(amount.T[..., None] * weights)
    origin = np.array([pixel[0] - lo[0], pixel[1] - lo[1], -2.5])

    def directions(altitude, azimuth):
        a, b = np.radians(altitude), np.radians(azimuth)
        east, north = np.cos(a) * np.sin(b), np.cos(a) * np.cos(b)
        horizontal = np.stack(np.broadcast_arrays(east, north), axis=-1) @ metric.T
        return np.column_stack((horizontal.reshape(-1, 2), np.sin(a).reshape(-1)))

    rays = directions(alt, az)
    sun = directions(np.array([sun_alt]), np.array([sun_az]))[0]
    result, transmission = _render(density, origin, rays, sun, base, opacity, show_grid)
    info = {
        "native_shape": [h, w],
        "pixel_window_xy": [lo.tolist(), hi.tolist()],
        "satellite_pixel_spacing_m": [dx, dy],
        "local_pixel_basis_km": np.linalg.inv(metric).tolist(),
        "vertical_edges_km": np.arange(2.5, 12.0).tolist(),
        "coverage_ratio": float(np.mean(valid)),
        "bt_warm_k": float(warm),
        "bt_cold_k": float(cold),
        "b16_redistribution": False,
        "geometry": "native pixel footprints linearized at observer; flat Earth",
        "horizontal_extent_km": 200,
    }
    if return_transmission:
        return result, info, transmission
    return result, info
