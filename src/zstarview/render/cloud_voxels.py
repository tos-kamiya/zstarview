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
CLOUD_AMOUNT_SUBTRACTION = 0.09
GEO_REFINEMENT_FACTOR = 3

SUNLIGHT_LEVELS = np.array([0.0, 0.01, 0.03, 0.10, 0.30, 1.0])
CLOUD_WHITENESS = np.array([0.18, 0.60, 0.85, 0.95, 0.99, 1.0])


def _validate_cloud_appearance(
    sunlight_mix, night_color_rgb, cloud_amount_subtract
) -> tuple[float, np.ndarray, float | None]:
    """Validate and normalize appearance options shared by both cloud sources."""
    sunlight_mix = float(sunlight_mix)
    if not np.isfinite(sunlight_mix):
        raise ValueError("sunlight_mix must be finite")
    sunlight_mix = float(np.clip(sunlight_mix, 0.0, 1.0))

    if cloud_amount_subtract is not None:
        cloud_amount_subtract = float(cloud_amount_subtract)
        if not np.isfinite(cloud_amount_subtract) or not (
            0.0 <= cloud_amount_subtract <= 1.0
        ):
            raise ValueError("cloud_amount_subtract must be finite and between 0 and 1")

    night_color_rgb = np.asarray(night_color_rgb, dtype=np.float32)
    if night_color_rgb.shape != (3,) or not np.all(np.isfinite(night_color_rgb)):
        raise ValueError("night_color_rgb must contain three finite RGB values")
    night_color_rgb = np.clip(night_color_rgb, 0.0, 1.0)
    return sunlight_mix, night_color_rgb, cloud_amount_subtract


def _altaz_to_directions(altitude_deg, azimuth_deg) -> np.ndarray:
    """Convert broadcastable altitude and azimuth angles to ENU unit vectors."""
    altitude_deg, azimuth_deg = np.broadcast_arrays(altitude_deg, azimuth_deg)
    altitude = np.radians(altitude_deg)
    azimuth = np.radians(azimuth_deg)
    east = np.cos(altitude) * np.sin(azimuth)
    north = np.cos(altitude) * np.cos(azimuth)
    horizontal = np.stack((east, north), axis=-1)
    return np.column_stack((horizontal.reshape(-1, 2), np.sin(altitude).reshape(-1)))


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
            pixel_from_enu[axis, 0] * offset[0] + pixel_from_enu[axis, 1] * offset[1]
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
        if 0 <= cell[0] < shape[0] and 0 <= cell[1] < shape[1] and cell[2] == layer:
            result.append((cell[0], cell[1], cell[2], end - t))
        t = end
    return result


@njit(cache=True)
def _segments(
    origin,
    direction,
    shape,
    pixel_origin,
    pixel_from_enu,
    layer_centers,
):
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
    sunlight_mix,
    night_color_rgb,
    base,
    opacity,
    pixel_origin,
    pixel_from_enu,
    layer_centers,
    layer_bases,
):
    shape = np.array(density.shape)
    light_fractions = np.full(density.shape, -1.0)
    # Recover the satellite-pixel amount before its distribution into layers.
    # Ambient night scattering follows this amount without directional shadows.
    night_amounts = np.minimum(density.sum(axis=2), 1.0)
    output = base.copy()
    transmissions = np.ones(len(rays), dtype=np.float32)
    for i in range(len(rays)):
        transmission = 1.0
        red = 0.0
        green = 0.0
        blue = 0.0
        for x, y, z, distance in _segments(
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
            if sunlight_mix > 0.0 and light_fractions[x, y, z] < 0:
                pixel_delta_x = x + 0.5 - pixel_origin[0]
                pixel_delta_y = y + 0.5 - pixel_origin[1]
                center_x = layer_centers[z, 0] + (
                    layer_bases[z, 0, 0] * pixel_delta_x
                    + layer_bases[z, 0, 1] * pixel_delta_y
                )
                center_y = layer_centers[z, 1] + (
                    layer_bases[z, 1, 0] * pixel_delta_x
                    + layer_bases[z, 1, 1] * pixel_delta_y
                )
                center = np.array([center_x, center_y, z + 0.5])
                tau = 0.0
                for sx, sy, sz, length in _segments(
                    center,
                    sun,
                    shape,
                    pixel_origin,
                    pixel_from_enu,
                    layer_centers,
                ):
                    tau += density[sx, sy, sz] * length * 1.8
                # Directional light remains active during the dusk transition.
                light_fractions[x, y, z] = _light_with_environment(1.0, tau)
            brightness = 0.0
            if sunlight_mix > 0.0:
                brightness = np.interp(
                    abs(sun[2]) * light_fractions[x, y, z],
                    SUNLIGHT_LEVELS,
                    CLOUD_WHITENESS,
                )
            alpha = 1.0 - np.exp(-2.4 * amount * distance * opacity)
            # Cloud amount controls extinction, not a separate color multiplier.
            # Thin cloud keeps a pale color while more background shines through.
            night_scale = (1.0 - sunlight_mix) * night_amounts[x, y]
            cell_red = sunlight_mix * brightness * 0.96
            cell_green = sunlight_mix * brightness * 0.975
            cell_blue = sunlight_mix * brightness
            cell_red += night_scale * night_color_rgb[0]
            cell_green += night_scale * night_color_rgb[1]
            cell_blue += night_scale * night_color_rgb[2]
            red += transmission * alpha * cell_red
            green += transmission * alpha * cell_green
            blue += transmission * alpha * cell_blue
            transmission *= 1.0 - alpha
            if transmission < 1e-5:
                break
        output[i, 0] = red + transmission * base[i, 0]
        output[i, 1] = green + transmission * base[i, 1]
        output[i, 2] = blue + transmission * base[i, 2]
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
        return (
            np.array([np.dot(offset, east_axis), np.dot(offset, north_axis)]) / 1000.0
        )

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
            bases[layer, :, axis] = pixel_basis_at_height(center, height_km, axis)
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
    sunlight_mix,
    night_color_rgb,
    opacity=0.85,
    height_layer_transform=False,
    return_transmission=False,
    cloud_amount_subtract=None,
):
    """Use raw pixels without interpolation and optionally vary grids by height.

    Height transforms approximate each 1-km slab as rectangular prisms using a
    local affine fit to satellite rays intersecting an altitude-offset ellipsoid.
    The prism footprint changes between slabs; its sides do not taper within a
    slab. B16 redistribution is omitted in this experiment.

    sunlight_mix blends directional daylight and ambient night components.
    night_color_rgb is scaled by the pixel's cloud amount before layer allocation.
    """
    sunlight_mix, night_color_rgb, cloud_amount_subtract = _validate_cloud_appearance(
        sunlight_mix, night_color_rgb, cloud_amount_subtract
    )
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
    raw_amount = _bt_to_weight(raw, warm, cold)
    if cloud_amount_subtract is None:
        amount = _suppress_low_cloud_weight(raw_amount)
    else:
        amount = np.clip(raw_amount - cloud_amount_subtract, 0.0, 1.0)
    weights = np.repeat(_blend_cloud_shell_weights(float(np.mean(amount))), 3) / 3
    density = np.ascontiguousarray(amount.T[..., None] * weights)
    pixel_origin = pixel - lo
    origin = np.array([0.0, 0.0, -2.5])

    rays = _altaz_to_directions(alt, az)
    sun = _altaz_to_directions(np.array([sun_alt]), np.array([sun_az]))[0]
    result, transmission = _render(
        density,
        origin,
        rays,
        sun,
        sunlight_mix,
        night_color_rgb,
        base,
        opacity,
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
        "sunlight_mix": sunlight_mix,
        "cloud_amount_subtract": cloud_amount_subtract,
        "night_color_rgb": night_color_rgb.tolist(),
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


def _refine_geo_cloud_amount(amount, factor, valid):
    """Interpolate every displayed Geo-satellite pixel into a finer grid.

    Keep missing source cells empty. Valid neighbors contribute in proportion
    to their bilinear weights, and edge samples clamp to the source extent.
    """
    h, w = amount.shape
    x = (np.arange(w * factor, dtype=np.float32) + 0.5) / factor - 0.5
    y = (np.arange(h * factor, dtype=np.float32) + 0.5) / factor - 0.5
    x0, y0 = np.floor(x).astype(np.intp), np.floor(y).astype(np.intp)
    wx, wy = x - x0, y - y0
    xa, xb = np.clip(x0, 0, w - 1), np.clip(x0 + 1, 0, w - 1)
    ya, yb = np.clip(y0, 0, h - 1), np.clip(y0 + 1, 0, h - 1)
    weighted_amount = np.zeros((h * factor, w * factor), dtype=np.float32)
    weight_total = np.zeros_like(weighted_amount)
    for row, row_weight in ((ya, 1.0 - wy), (yb, wy)):
        for col, col_weight in ((xa, 1.0 - wx), (xb, wx)):
            weight = row_weight[:, None] * col_weight[None, :] * valid[np.ix_(row, col)]
            weighted_amount += amount[np.ix_(row, col)] * weight
            weight_total += weight
    fine = np.zeros_like(weighted_amount)
    np.divide(weighted_amount, weight_total, out=fine, where=weight_total > 0)
    parent_valid = np.repeat(np.repeat(valid, factor, axis=0), factor, axis=1)
    fine[~parent_valid] = 0.0
    return fine


def shade_geo_satellite_voxels(
    source,
    lat,
    lon,
    alt,
    az,
    sun_alt,
    sun_az,
    base,
    *,
    sunlight_mix,
    night_color_rgb,
    opacity=0.85,
    height_layer_transform=False,
    return_transmission=False,
    cloud_amount_subtract=None,
    refinement_factor=GEO_REFINEMENT_FACTOR,
):
    """Render display-derived Geo-satellite cloud amounts in flat voxel layers."""
    from ..geosatellite.projection import _load_projection_inverse

    if height_layer_transform:
        raise ValueError("Geo-satellite voxel layers use one fixed horizontal grid")
    if refinement_factor not in (1, 3):
        raise ValueError("Geo-satellite refinement factor must be 1 or 3")
    sunlight_mix, night_color_rgb, cloud_amount_subtract = _validate_cloud_appearance(
        sunlight_mix, night_color_rgb, cloud_amount_subtract
    )

    cloud_amount = np.asarray(source.cloud_amount, dtype=np.float32)
    valid_mask = np.asarray(source.valid_mask, dtype=bool)
    if cloud_amount.ndim != 2 or valid_mask.shape != cloud_amount.shape:
        raise ValueError("Geo-satellite voxel source arrays must be matching 2D grids")
    h, w = cloud_amount.shape
    if h < 2 or w < 2:
        raise ValueError("Geo-satellite voxel source grid is too small")

    projection = _load_projection_inverse(source.grid_npz)
    pixel_x, pixel_y = projection.lonlat_to_pixel(
        np.asarray(float(lon)), np.asarray(float(lat))
    )
    pixel = np.asarray([float(pixel_x), float(pixel_y)], dtype=np.float64)
    lat_rad = np.radians(float(lat))
    dlon = np.degrees(1.0 / (6371.0 * max(1e-6, abs(np.cos(lat_rad)))))
    dlat = np.degrees(1.0 / 6371.0)
    east_x, east_y = projection.lonlat_to_pixel(
        np.asarray(float(lon) + dlon), np.asarray(float(lat))
    )
    north_x, north_y = projection.lonlat_to_pixel(
        np.asarray(float(lon)), np.asarray(float(lat) + dlat)
    )
    pixel_from_enu = np.asarray(
        [
            [float(east_x) - pixel[0], float(north_x) - pixel[0]],
            [float(east_y) - pixel[1], float(north_y) - pixel[1]],
        ],
        dtype=np.float64,
    )
    det = float(np.linalg.det(pixel_from_enu))
    if not np.all(np.isfinite(pixel_from_enu)) or abs(det) < 1e-12:
        raise ValueError("Geo-satellite local pixel basis is singular")
    layer_basis = np.linalg.inv(pixel_from_enu)
    # One common projected grid is used for all layers; satellite parallax is
    # intentionally omitted from this display-oriented Geo-satellite path.
    pixel_from_enu_layers = np.repeat(pixel_from_enu[None, :, :], 9, axis=0)
    layer_bases = np.repeat(layer_basis[None, :, :], 9, axis=0)
    layer_centers = np.zeros((9, 2), dtype=np.float64)

    corners = np.asarray(
        [[-200.0, -200.0], [-200.0, 200.0], [200.0, -200.0], [200.0, 200.0]],
        dtype=np.float64,
    )
    corner_pixels = pixel[None, :] + corners @ pixel_from_enu.T
    lo = np.maximum(0, np.floor(np.min(corner_pixels, axis=0)).astype(int))
    hi = np.minimum(
        np.asarray([w, h]),
        np.ceil(np.max(corner_pixels, axis=0)).astype(int) + 1,
    )
    if np.any(hi <= lo):
        raise ValueError("Geo-satellite voxel view does not overlap source image")
    raw_amount = cloud_amount[lo[1] : hi[1], lo[0] : hi[0]]
    valid = valid_mask[lo[1] : hi[1], lo[0] : hi[0]] & np.isfinite(raw_amount)
    amount = np.zeros(raw_amount.shape, dtype=np.float32)
    if cloud_amount_subtract is None:
        amount[valid] = np.clip(raw_amount[valid], 0.0, 1.0)
    else:
        amount[valid] = np.clip(raw_amount[valid] - cloud_amount_subtract, 0.0, 1.0)
    scene_amount = float(np.mean(amount[valid])) if np.any(valid) else 0.0
    group_weights = _blend_cloud_shell_weights(scene_amount)
    weights = np.repeat(np.asarray(group_weights, dtype=np.float32), 3) / 3.0
    pixel_origin = pixel - lo
    if refinement_factor > 1:
        amount = _refine_geo_cloud_amount(amount, refinement_factor, valid)
        pixel_origin = pixel_origin * refinement_factor
        pixel_from_enu_layers = pixel_from_enu_layers * refinement_factor
        layer_bases = layer_bases / refinement_factor
    else:
        refinement_factor = 1
    density = np.ascontiguousarray(amount.T[..., None] * weights)
    origin = np.asarray([0.0, 0.0, -0.5], dtype=np.float64)

    rays = _altaz_to_directions(alt, az)
    sun = _altaz_to_directions(np.asarray([sun_alt]), np.asarray([sun_az]))[0]
    result, transmission = _render(
        density,
        origin,
        rays,
        sun,
        sunlight_mix,
        night_color_rgb,
        base,
        float(np.clip(opacity, 0.0, 1.0)),
        pixel_origin,
        pixel_from_enu_layers,
        layer_centers,
        layer_bases,
    )
    info = {
        "native_shape": [h, w],
        "pixel_window_xy": [lo.tolist(), hi.tolist()],
        "local_pixel_basis_km": layer_basis.tolist(),
        "refinement_factor": refinement_factor,
        "refined_pixel_basis_km": layer_bases[0].tolist(),
        "height_layer_transform": False,
        "height_layer_geometry": [
            {
                "height_km": 1.0 + layer,
                "center_offset_enu_km": [0.0, 0.0],
                "pixel_basis_enu_km": layer_basis.tolist(),
            }
            for layer in range(9)
        ],
        "vertical_edges_km": np.arange(0.5, 10.0).tolist(),
        "cloud_layer_weights": weights.tolist(),
        "scene_cloud_amount": scene_amount,
        "coverage_ratio": float(np.mean(valid)),
        "source_kind": str(source.kind),
        "source_time_utc": source.time_utc.isoformat(),
        "sunlight_mix": sunlight_mix,
        "cloud_amount_subtract": cloud_amount_subtract,
        "night_color_rgb": night_color_rgb.tolist(),
        "geometry": "fixed georeferenced horizontal grid shared by all cloud layers",
        "horizontal_extent_km": 200,
    }
    if return_transmission:
        return result, info, transmission
    return result, info
