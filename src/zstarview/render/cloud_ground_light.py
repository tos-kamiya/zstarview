"""Worker-only direct geographic sampling of vertical ground illumination."""

from __future__ import annotations

import logging

import numpy as np
import rasterio
from pyproj import Geod, Transformer
from rasterio.windows import Window

from ..cloud_night_lighting import GROUND_RADIANCE_SCALE
from ..night_light_source import _ensure_night_light_tiles, _tile_name_for_latlon

logger = logging.getLogger(__name__)


def _sample_raster(dataset, longitudes, latitudes):
    """Read bounded windows, then index the original raster pixels."""
    xs, ys = longitudes, latitudes
    if dataset.crs is not None and dataset.crs != rasterio.crs.CRS.from_epsg(4326):
        xs, ys = Transformer.from_crs(
            "EPSG:4326", dataset.crs, always_xy=True
        ).transform(xs, ys)
    rows, cols = rasterio.transform.rowcol(dataset.transform, xs, ys)
    rows, cols = np.asarray(rows), np.asarray(cols)
    valid = (rows >= 0) & (rows < dataset.height) & (cols >= 0) & (cols < dataset.width)
    result = np.zeros(len(rows), dtype=np.float32)
    if not np.any(valid):
        return result
    indices = np.flatnonzero(valid)
    # Bounded blocks also handle polar footprints without reading a huge
    # bounding rectangle across a 90-degree tile.
    blocks_across = (dataset.width + 511) // 512
    codes = (rows[valid] // 512) * blocks_across + cols[valid] // 512
    order = np.argsort(codes)
    groups = np.split(order, np.flatnonzero(np.diff(codes[order])) + 1)
    for group in groups:
        selected = indices[group]
        row0 = int(rows[selected[0]] // 512) * 512
        col0 = int(cols[selected[0]] // 512) * 512
        block = dataset.read(
            1, window=Window(
                col0, row0, min(512, dataset.width - col0), min(512, dataset.height - row0)
            ), masked=True
        )
        values = np.asarray(block[rows[selected] - row0, cols[selected] - col0].filled(0))
        result[selected] = np.where(np.isfinite(values), np.maximum(values, 0), 0)
    return result


def sample_cloud_ground_light(
    lat, lon, shape, pixel_origin, layer_centers, layer_bases
) -> np.ndarray:
    """Sample the surface directly under each physical cell center.

    Run only inside the one-shot cloud worker: downloads, raster I/O, projection
    and NumPy work share its process isolation and overall render deadline.
    Tile failures are independent and contribute zero illumination.
    """
    ground = np.zeros(shape, dtype=np.float32)
    geod = Geod(ellps="WGS84")
    xx, yy = np.meshgrid(
        np.arange(shape[0]) + 0.5 - pixel_origin[0],
        np.arange(shape[1]) + 0.5 - pixel_origin[1],
        indexing="ij",
    )
    tile_paths = {}
    failed_tiles: set[str] = set()
    for layer in range(shape[2]):
        basis = layer_bases[layer]
        east = layer_centers[layer, 0] + basis[0, 0] * xx + basis[0, 1] * yy
        north = layer_centers[layer, 1] + basis[1, 0] * xx + basis[1, 1] * yy
        distance = np.hypot(east, north).ravel() * 1000.0
        azimuth = np.degrees(np.arctan2(east, north)).ravel()
        longitudes, latitudes, _ = geod.fwd(
            [float(lon)] * len(distance), [float(lat)] * len(distance),
            azimuth.tolist(), distance.tolist(),
        )
        longitudes, latitudes = np.asarray(longitudes), np.asarray(latitudes)
        # Vectorized tile names, with the same wrapping and equator convention
        # as the existing night-light sampler.
        columns = np.clip(((longitudes + 180.0) % 360.0 / 90.0).astype(int), 0, 3)
        rows = (latitudes < 0.0).astype(int)
        codes = columns + 4 * rows
        values = np.zeros(len(distance), dtype=np.float32)
        for code in np.unique(codes):
            mask = codes == code
            first = int(np.flatnonzero(mask)[0])
            tile = _tile_name_for_latlon(float(latitudes[first]), float(longitudes[first]))
            if tile in failed_tiles:
                continue
            try:
                if tile not in tile_paths:
                    tile_paths.update(_ensure_night_light_tiles(
                        tile_names={tile}, timeout_s=3.0, download_timeout_s=5.0
                    ))
                with rasterio.open(tile_paths[tile]) as dataset:
                    values[mask] = _sample_raster(dataset, longitudes[mask], latitudes[mask])
            except Exception as exc:
                failed_tiles.add(tile)
                logger.warning("Cloud ground light tile %s unavailable: %s", tile, exc)
        ground[:, :, layer] = (
            1.0 - np.exp(-values / GROUND_RADIANCE_SCALE)
        ).reshape(xx.shape)
    return ground
