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


ENVIRONMENT_LIGHT_FRACTION = 0.03

SUNLIGHT_LEVELS = np.array([0.0, 0.01, 0.03, 0.10, 0.30, 1.0])
CLOUD_WHITENESS = np.array([0.18, 0.60, 0.85, 0.95, 0.99, 1.0])


@njit(cache=True)
def _segments_in_layer(
    origin,
    direction,
    shape,
    layer,
    pixel_origin,
    pixel_from_enu,
    layer_center_enu,
):
    """Intersect one ray with a layer's locally affine voxel prisms."""
    result = []
    grid_origin = np.empty(3)
    grid_direction = np.empty(3)
    offset = origin[:2] - layer_center_enu
    for axis in range(2):
        grid_origin[axis] = pixel_origin[axis] + (
            pixel_from_enu[axis, 0] * offset[0]
            + pixel_from_enu[axis, 1] * offset[1]
        )
        grid_direction[axis] = (
            pixel_from_enu[axis, 0] * direction[0]
            + pixel_from_enu[axis, 1] * direction[1]
        )
    grid_origin[2] = origin[2]
    grid_direction[2] = direction[2]
    lower = np.array([0.0, 0.0, float(layer)])
    upper = np.array([float(shape[0]), float(shape[1]), float(layer + 1)])
    near, far = 0.0, 1000.0
    for axis in range(3):
        if abs(grid_direction[axis]) < 1e-12:
            if grid_origin[axis] < lower[axis] or grid_origin[axis] >= upper[axis]:
                return result
        else:
            a = (lower[axis] - grid_origin[axis]) / grid_direction[axis]
            b = (upper[axis] - grid_origin[axis]) / grid_direction[axis]
            near = max(near, min(a, b))
            far = min(far, max(a, b))
    if far <= near:
        return result
    t = near
    while t < far - 1e-8:
        p = grid_origin + grid_direction * (t + 1e-7)
        cell = np.floor(p).astype(np.int64)
        end = far
        for axis in range(2):
            if grid_direction[axis] > 1e-12:
                end = min(
                    end,
                    (cell[axis] + 1 - grid_origin[axis]) / grid_direction[axis],
                )
            elif grid_direction[axis] < -1e-12:
                end = min(
                    end,
                    (cell[axis] - grid_origin[axis]) / grid_direction[axis],
                )
        if end <= t:
            break
        if (
            0 <= cell[0] < shape[0]
            and 0 <= cell[1] < shape[1]
            and cell[2] == layer
        ):
            result.append(
                (
                    cell[0],
                    cell[1],
                    cell[2],
                    end - t,
                    t,
                )
            )
        t = end
    return result


@njit(cache=True)
def _segments(origin, direction, shape, pixel_origin, pixel_from_enu, layer_centers):
    """Trace ordered intervals through layer-specific rectangular prisms."""
    result = []
    if direction[2] >= 0.0:
        for layer in range(shape[2]):
            for segment in _segments_in_layer(
                origin,
                direction,
                shape,
                layer,
                pixel_origin,
                pixel_from_enu[layer],
                layer_centers[layer],
            ):
                result.append(segment)
    else:
        for layer in range(shape[2] - 1, -1, -1):
            for segment in _segments_in_layer(
                origin,
                direction,
                shape,
                layer,
                pixel_origin,
                pixel_from_enu[layer],
                layer_centers[layer],
            ):
                result.append(segment)
    return result


@njit(cache=True)
def _light_with_environment(clear_sunlight: float, optical_depth: float) -> float:
    """Closed form of repeated L = L*T + environment*(1-T) steps."""
    transmission = np.exp(-optical_depth)
    environment = clear_sunlight * ENVIRONMENT_LIGHT_FRACTION
    return clear_sunlight * transmission + environment * (1.0 - transmission)


@njit(cache=True)
def _render(
    density,
    origin,
    rays,
    sun,
    base,
    opacity,
    show_grid,
    pixel_origin,
    pixel_from_enu,
    layer_centers,
    layer_bases,
):
    shape = np.array(density.shape)
    lights = np.full(density.shape, -1.0)
    output = base.copy()
    transmissions = np.ones(len(rays), dtype=np.float32)
    for i in range(len(rays)):
        transmission = 1.0
        color = np.zeros(3)
        for x, y, z, distance, entry in _segments(
            origin,
            rays[i],
            shape,
            pixel_origin,
            pixel_from_enu,
            layer_centers,
        ):
            amount = density[x, y, z]
            if amount <= 0:
                continue
            if lights[x, y, z] < 0:
                tau = 0.0
                if sun[2] > 0:
                    pixel_delta = np.array(
                        [x + 0.5 - pixel_origin[0], y + 0.5 - pixel_origin[1]]
                    )
                    center_enu = np.empty(2)
                    for axis in range(2):
                        center_enu[axis] = layer_centers[z, axis] + (
                            layer_bases[z, axis, 0] * pixel_delta[0]
                            + layer_bases[z, axis, 1] * pixel_delta[1]
                        )
                    center = np.array([center_enu[0], center_enu[1], z + 0.5])
                    for sx, sy, sz, length, _ in _segments(
                        center,
                        sun,
                        shape,
                        pixel_origin,
                        pixel_from_enu,
                        layer_centers,
                    ):
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
                entry_world = origin + rays[i] * entry
                entry_offset = entry_world[:2] - layer_centers[z]
                entry_point = np.array(
                    [
                        pixel_origin[0]
                        + pixel_from_enu[z, 0, 0] * entry_offset[0]
                        + pixel_from_enu[z, 0, 1] * entry_offset[1],
                        pixel_origin[1]
                        + pixel_from_enu[z, 1, 0] * entry_offset[0]
                        + pixel_from_enu[z, 1, 1] * entry_offset[1],
                        entry_world[2],
                    ]
                )
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


def _height_layer_geometry(area, lat, lon, pixel, dx, dy, fallback_basis):
    """Approximate each slab's satellite-pixel grid on an altitude shell."""
    params = area.crs.to_dict()
    if params.get("proj") != "geos":
        raise ValueError("height-aware voxel geometry requires a geostationary grid")
    semi_major = float(area.crs.ellipsoid.semi_major_metre)
    semi_minor = float(area.crs.ellipsoid.semi_minor_metre)
    satellite_height = float(params["h"])
    satellite_lon = np.radians(float(params["lon_0"]))
    satellite_radius = semi_major + satellite_height
    satellite = np.array(
        [
            satellite_radius * np.cos(satellite_lon),
            satellite_radius * np.sin(satellite_lon),
            0.0,
        ],
        dtype=np.float64,
    )
    to_geographic = Transformer.from_crs(area.crs, "EPSG:4326", always_xy=True)
    to_geocentric = Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)
    observer_ecef = np.array(to_geocentric.transform(lon, lat, 0.0))
    lat_radians, lon_radians = np.radians([lat, lon])
    east_axis = np.array([-np.sin(lon_radians), np.cos(lon_radians), 0.0])
    north_axis = np.array(
        [
            -np.sin(lat_radians) * np.cos(lon_radians),
            -np.sin(lat_radians) * np.sin(lon_radians),
            np.cos(lat_radians),
        ]
    )
    xmin, ymin, xmax, ymax = area.area_extent

    def point_on_shell(pixel_xy, height_km):
        proj_x = xmin + float(pixel_xy[0]) * dx
        proj_y = ymax - float(pixel_xy[1]) * dy
        point_lon, point_lat = to_geographic.transform(proj_x, proj_y)
        if not np.isfinite(point_lon) or not np.isfinite(point_lat):
            raise ValueError("pixel is outside the geostationary projection")
        ground = np.array(to_geocentric.transform(point_lon, point_lat, 0.0))
        ray = ground - satellite
        axis_a = semi_major + height_km * 1000.0
        axis_b = semi_minor + height_km * 1000.0
        inv_a2 = 1.0 / (axis_a * axis_a)
        inv_b2 = 1.0 / (axis_b * axis_b)
        qa = (ray[0] * ray[0] + ray[1] * ray[1]) * inv_a2 + ray[2] ** 2 * inv_b2
        qb = 2.0 * (
            (satellite[0] * ray[0] + satellite[1] * ray[1]) * inv_a2
            + satellite[2] * ray[2] * inv_b2
        )
        qc = (satellite[0] ** 2 + satellite[1] ** 2) * inv_a2
        qc += satellite[2] ** 2 * inv_b2 - 1.0
        discriminant = qb * qb - 4.0 * qa * qc
        if discriminant < 0.0:
            raise ValueError("satellite ray does not intersect the altitude shell")
        fraction = (-qb - np.sqrt(discriminant)) / (2.0 * qa)
        point = satellite + fraction * ray
        offset = point - observer_ecef
        return np.array([np.dot(offset, east_axis), np.dot(offset, north_axis)]) / 1000.0

    centers = np.empty((9, 2), dtype=np.float64)
    bases = np.empty((9, 2, 2), dtype=np.float64)

    def pixel_basis_at_height(center, height_km, axis):
        delta = np.zeros(2, dtype=np.float64)
        delta[axis] = 1.0
        for sign in (1.0, -1.0):
            try:
                neighbor = point_on_shell(pixel + sign * delta, height_km)
            except ValueError:
                continue
            return (neighbor - center) / sign
        return fallback_basis[:, axis]

    for layer in range(9):
        height_km = 3.0 + layer
        center = point_on_shell(pixel, height_km)
        centers[layer] = center
        for axis in range(2):
            bases[layer, :, axis] = pixel_basis_at_height(
                center, height_km, axis
            )
    return centers, bases


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
    height_layer_transform=False,
    return_transmission=False,
):
    """Use raw pixels without interpolation and optionally vary grids by height.

    Height transforms approximate each 1-km slab as rectangular prisms using a
    local affine fit to satellite rays intersecting an altitude-offset ellipsoid.
    The prism footprint changes between slabs; its sides do not taper within a
    slab. B16 redistribution is omitted in this experiment.
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
    surface_basis = np.linalg.inv(metric)
    if height_layer_transform:
        layer_centers, layer_bases = _height_layer_geometry(
            area, lat, lon, pixel, dx, dy, surface_basis
        )
    else:
        layer_centers = np.zeros((9, 2), dtype=np.float64)
        layer_bases = np.repeat(surface_basis[None, :, :], 9, axis=0)
    pixel_from_enu = np.linalg.inv(layer_bases)
    # Include the per-height parallax offsets while keeping a 200-km local view.
    corners = np.array(
        [[-200.0, -200.0], [-200.0, 200.0], [200.0, -200.0], [200.0, 200.0]]
    )
    layer_pixel_centers = pixel[None, :] - np.einsum(
        "lij,lj->li", pixel_from_enu, layer_centers
    )
    layer_bounds = layer_pixel_centers[:, None, :] + np.einsum(
        "lij,cj->lci", pixel_from_enu, corners
    )
    pixel_lo = np.floor(np.min(layer_bounds, axis=(0, 1))).astype(int)
    pixel_hi = np.ceil(np.max(layer_bounds, axis=(0, 1))).astype(int) + 1
    lo = np.maximum(0, pixel_lo)
    hi = np.minimum([w, h], pixel_hi)
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
    pixel_origin = pixel - lo
    origin = np.array([0.0, 0.0, -2.5])

    def directions(altitude, azimuth):
        a, b = np.radians(altitude), np.radians(azimuth)
        east, north = np.cos(a) * np.sin(b), np.cos(a) * np.cos(b)
        horizontal = np.stack(np.broadcast_arrays(east, north), axis=-1)
        return np.column_stack((horizontal.reshape(-1, 2), np.sin(a).reshape(-1)))

    rays = directions(alt, az)
    sun = directions(np.array([sun_alt]), np.array([sun_az]))[0]
    result, transmission = _render(
        density,
        origin,
        rays,
        sun,
        base,
        opacity,
        show_grid,
        pixel_origin,
        pixel_from_enu,
        layer_centers,
        layer_bases,
    )
    info = {
        "native_shape": [h, w],
        "pixel_window_xy": [lo.tolist(), hi.tolist()],
        "satellite_pixel_spacing_m": [dx, dy],
        "local_pixel_basis_km": layer_bases[0].tolist(),
        "height_layer_transform": bool(height_layer_transform),
        "height_layer_geometry": [
            {
                "height_km": 3.0 + layer,
                "center_offset_enu_km": layer_centers[layer].tolist(),
                "pixel_basis_enu_km": layer_bases[layer].tolist(),
            }
            for layer in range(9)
        ],
        "vertical_edges_km": np.arange(2.5, 12.0).tolist(),
        "coverage_ratio": float(np.mean(valid)),
        "bt_warm_k": float(warm),
        "bt_cold_k": float(cold),
        "b16_redistribution": False,
        "geometry": (
            "per-height local affine fit to satellite rays and altitude ellipsoid"
            if height_layer_transform
            else "native pixel footprints linearized at observer; flat Earth"
        ),
        "horizontal_extent_km": 200,
    }
    if return_transmission:
        return result, info, transmission
    return result, info
