# -*- coding: utf-8 -*-
"""Helpers for direct Himawari ISatSS tile ingestion without Satpy."""

from __future__ import annotations

import datetime as dt
import math
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import xarray as xr
from pyproj import CRS, Transformer

from ..geo_area import GeoArea
from ..projectors.az import altaz_to_dir_ecef, geodetic_to_ecef
from ._s3_io import list_s3_keys

DATA_VAR = "Sectorized_CMI"
GRID_VAR = "fixedgrid_projection"
_RE_TILE = re.compile(r"-T(\d{3})_")
_HIMA_BUCKETS = ("noaa-himawari9", "noaa-himawari8")
_PREFIX_ROOT = "AHI-L2-FLDK-ISatSS"

warnings.filterwarnings(
    "ignore",
    message="You will likely lose important projection information when converting to a PROJ string",
    category=UserWarning,
    module="pyproj.crs.crs",
)


@dataclass(frozen=True)
class TemplateMeta:
    bucket: str
    tile_path: Path
    tile_height: int
    tile_width: int
    product_rows: int
    product_cols: int
    tile_count: int
    x0: float
    y0: float
    x_step: float
    y_step: float
    geos_scale: float
    crs: CRS


@dataclass(frozen=True)
class TileRecord:
    token: str
    row_offset: int
    col_offset: int
    x_min: float
    x_max: float
    y_min: float
    y_max: float


@dataclass(frozen=True)
class ObserverTileSelection:
    selected_keys: list[str]
    render_tiles: list[TileRecord]
    equator_tiles: list[TileRecord]
    far_missing_render_tiles: list[TileRecord]
    near_missing_render_tiles: list[TileRecord]
    missing_equator_tiles: list[TileRecord]


def format_prefix(when_utc: dt.datetime) -> str:
    slot = when_utc.astimezone(dt.timezone.utc)
    hm = f"{slot.hour:02d}{(slot.minute // 10) * 10:02d}"
    return f"{_PREFIX_ROOT}/{slot.year:04d}/{slot.month:02d}/{slot.day:02d}/{hm}/"


def extract_tile_token(name: str) -> str:
    match = _RE_TILE.search(name)
    if not match:
        raise ValueError(f"Could not extract tile token from {name}")
    return match.group(1)


def find_matching_keys(when_utc: dt.datetime, *, satellite: str, product: str, timeout_s: float | None = None) -> tuple[str, list[str]]:
    prefix = format_prefix(when_utc)
    for bucket in _HIMA_BUCKETS:
        keys = list_s3_keys(
            bucket=bucket,
            prefix=prefix,
            satellite=satellite,
            product=product,
            time_utc=when_utc,
            timeout_s=timeout_s,
        )
        matched = sorted(key for key in keys if "M1C13" in Path(key).name and key.endswith(".nc"))
        if matched:
            return bucket, matched
    raise FileNotFoundError(f"No Himawari ISatSS M1C13 keys found for {prefix}")


def load_template_from_tile(tile_path: Path, *, bucket: str) -> TemplateMeta:
    with xr.open_dataset(tile_path) as ds:
        tile_height = int(ds.sizes["y"])
        tile_width = int(ds.sizes["x"])
        product_rows = int(ds.attrs["product_rows"])
        product_cols = int(ds.attrs["product_columns"])
        tile_count = int(ds.attrs["number_product_tiles"])
        x_step = float(ds.x.values[1] - ds.x.values[0])
        y_step = float(ds.y.values[1] - ds.y.values[0])
        x0 = float(ds.x.values[0] - int(ds.attrs["tile_column_offset"]) * x_step)
        y0 = float(ds.y.values[0] - int(ds.attrs["tile_row_offset"]) * y_step)
        crs = CRS.from_cf(dict(ds[GRID_VAR].attrs))
        geos_scale = float(ds[GRID_VAR].attrs["perspective_point_height"]) * 1e-6
        return TemplateMeta(
            bucket=bucket,
            tile_path=tile_path,
            tile_height=tile_height,
            tile_width=tile_width,
            product_rows=product_rows,
            product_cols=product_cols,
            tile_count=tile_count,
            x0=x0,
            y0=y0,
            x_step=x_step,
            y_step=y_step,
            geos_scale=geos_scale,
            crs=crs,
        )


def generate_sparse_layout(meta: TemplateMeta) -> list[TileRecord]:
    grid_rows = meta.product_rows // meta.tile_height
    grid_cols = meta.product_cols // meta.tile_width
    if meta.product_rows % meta.tile_height or meta.product_cols % meta.tile_width:
        raise ValueError("Product size is not divisible by tile size")

    if meta.tile_count == grid_rows * grid_cols:
        row_counts = [grid_cols] * grid_rows
    elif grid_rows == 10 and grid_cols == 10 and meta.tile_count == 88:
        row_counts = [6, 8, 10, 10, 10, 10, 10, 10, 8, 6]
    else:
        raise ValueError(f"Unsupported sparse tile layout: grid={grid_rows}x{grid_cols}, tile_count={meta.tile_count}")

    records: list[TileRecord] = []
    token = 1
    for row_index, count in enumerate(row_counts):
        start_col = (grid_cols - count) // 2
        for local_col in range(count):
            row_offset = row_index * meta.tile_height
            col_offset = (start_col + local_col) * meta.tile_width
            x_a = meta.x0 + col_offset * meta.x_step
            x_b = meta.x0 + (col_offset + meta.tile_width) * meta.x_step
            y_a = meta.y0 + row_offset * meta.y_step
            y_b = meta.y0 + (row_offset + meta.tile_height) * meta.y_step
            records.append(
                TileRecord(
                    token=f"{token:03d}",
                    row_offset=row_offset,
                    col_offset=col_offset,
                    x_min=min(x_a, x_b) * meta.geos_scale,
                    x_max=max(x_a, x_b) * meta.geos_scale,
                    y_min=min(y_a, y_b) * meta.geos_scale,
                    y_max=max(y_a, y_b) * meta.geos_scale,
                )
            )
            token += 1
    if len(records) != meta.tile_count:
        raise ValueError(f"Tile layout mismatch: expected {meta.tile_count}, got {len(records)}")
    return records


def _intersect_cloud_shell(observer_xyz: np.ndarray, direction_xyz: np.ndarray, shell_radius_km: float) -> tuple[float, float, float] | None:
    b_quad = 2.0 * float(np.dot(observer_xyz, direction_xyz))
    c_quad = float(np.dot(observer_xyz, observer_xyz) - shell_radius_km * shell_radius_km)
    disc = b_quad * b_quad - 4.0 * c_quad
    if disc < 0.0:
        return None
    sqrt_disc = math.sqrt(max(disc, 0.0))
    t1 = (-b_quad - sqrt_disc) / 2.0
    t2 = (-b_quad + sqrt_disc) / 2.0
    t = t1 if t1 > 1e-6 else (t2 if t2 > 1e-6 else None)
    if t is None:
        return None
    point = observer_xyz + direction_xyz * t
    return float(point[0]), float(point[1]), float(point[2])


def _visible_boundary_lonlat(*, lat_deg: float, lon_deg: float, cloud_shell_km: float, azimuth_samples: int) -> tuple[np.ndarray, np.ndarray]:
    observer_xyz = np.asarray(geodetic_to_ecef(lat_deg, lon_deg), dtype=np.float64)
    lons: list[float] = [float(lon_deg)]
    lats: list[float] = [float(lat_deg)]
    for az_deg in np.linspace(0.0, 360.0, azimuth_samples, endpoint=False):
        direction = np.asarray(altaz_to_dir_ecef(0.0, float(az_deg), lat_deg, lon_deg), dtype=np.float64)
        hit = _intersect_cloud_shell(observer_xyz, direction, cloud_shell_km)
        if hit is None:
            continue
        x, y, z = hit
        lons.append(math.degrees(math.atan2(y, x)))
        lats.append(math.degrees(math.atan2(z, math.hypot(x, y))))
    return np.asarray(lons, dtype=np.float64), np.asarray(lats, dtype=np.float64)


def _point_in_polygon(x: float, y: float, poly_x: np.ndarray, poly_y: np.ndarray) -> bool:
    inside = False
    j = poly_x.size - 1
    for i in range(poly_x.size):
        yi, yj = poly_y[i], poly_y[j]
        xi, xj = poly_x[i], poly_x[j]
        intersects = ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / ((yj - yi) or 1e-12) + xi)
        if intersects:
            inside = not inside
        j = i
    return inside


def _tile_intersects_polygon(tile: TileRecord, poly_x: np.ndarray, poly_y: np.ndarray) -> bool:
    if poly_x.size == 0:
        return False
    if poly_x.max() < tile.x_min or poly_x.min() > tile.x_max or poly_y.max() < tile.y_min or poly_y.min() > tile.y_max:
        return False
    center = ((tile.x_min + tile.x_max) * 0.5, (tile.y_min + tile.y_max) * 0.5)
    if _point_in_polygon(center[0], center[1], poly_x, poly_y):
        return True
    corners = [
        (tile.x_min, tile.y_min),
        (tile.x_min, tile.y_max),
        (tile.x_max, tile.y_min),
        (tile.x_max, tile.y_max),
    ]
    if any(_point_in_polygon(x, y, poly_x, poly_y) for x, y in corners):
        return True
    boundary_hit = (
        (poly_x >= tile.x_min)
        & (poly_x <= tile.x_max)
        & (poly_y >= tile.y_min)
        & (poly_y <= tile.y_max)
    )
    return bool(np.any(boundary_hit))


def _expand_selection(tokens: set[str], records: Sequence[TileRecord], margin_tiles: int) -> set[str]:
    if margin_tiles <= 0:
        return tokens
    row_stride = min(abs(a.row_offset - b.row_offset) for a in records for b in records if a.row_offset != b.row_offset)
    col_stride = min(abs(a.col_offset - b.col_offset) for a in records for b in records if a.col_offset != b.col_offset)
    by_token = {record.token: record for record in records}
    selected = set(tokens)
    for _ in range(margin_tiles):
        grown = set(selected)
        for token in list(selected):
            record = by_token[token]
            for other in records:
                if abs(other.row_offset - record.row_offset) <= row_stride and abs(other.col_offset - record.col_offset) <= col_stride:
                    grown.add(other.token)
        selected = grown
    return selected


def select_needed_tiles(
    *,
    lat_deg: float,
    lon_deg: float,
    meta: TemplateMeta,
    cloud_shell_km: float,
    azimuth_samples: int,
    margin_tiles: int,
) -> tuple[list[TileRecord], np.ndarray, np.ndarray]:
    lons, lats = _visible_boundary_lonlat(
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        cloud_shell_km=cloud_shell_km,
        azimuth_samples=azimuth_samples,
    )
    to_proj = Transformer.from_crs("EPSG:4326", meta.crs, always_xy=True)
    poly_x, poly_y = to_proj.transform(lons, lats)
    finite = np.isfinite(poly_x) & np.isfinite(poly_y)
    poly_x = np.asarray(poly_x[finite], dtype=np.float64)
    poly_y = np.asarray(poly_y[finite], dtype=np.float64)
    records = generate_sparse_layout(meta)
    selected_tokens = {record.token for record in records if _tile_intersects_polygon(record, poly_x, poly_y)}
    selected_tokens = _expand_selection(selected_tokens, records, margin_tiles)
    selected = [record for record in records if record.token in selected_tokens]
    selected.sort(key=lambda item: item.token)
    return selected, poly_x, poly_y


def select_equator_tiles(
    *,
    lon_center_deg: float,
    meta: TemplateMeta,
    delta_lon: float = 60.0,
    equator_lat: float = 0.0,
    step_deg: float = 1.0,
    margin_tiles: int = 0,
) -> tuple[list[TileRecord], np.ndarray, np.ndarray]:
    lons = np.arange(lon_center_deg - delta_lon, lon_center_deg + delta_lon + step_deg, step_deg, dtype=np.float64)
    lats = np.full_like(lons, float(equator_lat), dtype=np.float64)
    to_proj = Transformer.from_crs("EPSG:4326", meta.crs, always_xy=True)
    poly_x, poly_y = to_proj.transform(lons, lats)
    finite = np.isfinite(poly_x) & np.isfinite(poly_y)
    poly_x = np.asarray(poly_x[finite], dtype=np.float64)
    poly_y = np.asarray(poly_y[finite], dtype=np.float64)
    records = generate_sparse_layout(meta)
    selected_tokens = {record.token for record in records if _tile_intersects_polygon(record, poly_x, poly_y)}
    selected_tokens = _expand_selection(selected_tokens, records, margin_tiles)
    selected = [record for record in records if record.token in selected_tokens]
    selected.sort(key=lambda item: item.token)
    return selected, poly_x, poly_y


def select_equator_band_tiles(
    *,
    lon_center_deg: float,
    meta: TemplateMeta,
    delta_lon: float = 60.0,
    equator_lat: float = 0.0,
    equator_lat_half_band_deg: float = 5.0,
    step_deg: float = 1.0,
    margin_tiles: int = 0,
) -> tuple[list[TileRecord], np.ndarray, np.ndarray]:
    lat_half_band = max(0.0, float(equator_lat_half_band_deg))
    lats = np.arange(
        float(equator_lat) - lat_half_band,
        float(equator_lat) + lat_half_band + float(step_deg),
        float(step_deg),
        dtype=np.float64,
    )
    records = generate_sparse_layout(meta)
    selected_tokens: set[str] = set()
    poly_x_parts: list[np.ndarray] = []
    poly_y_parts: list[np.ndarray] = []
    to_proj = Transformer.from_crs("EPSG:4326", meta.crs, always_xy=True)
    for lat in lats:
        lons = np.arange(
            lon_center_deg - delta_lon,
            lon_center_deg + delta_lon + step_deg,
            step_deg,
            dtype=np.float64,
        )
        lats_arr = np.full_like(lons, float(lat), dtype=np.float64)
        poly_x, poly_y = to_proj.transform(lons, lats_arr)
        finite = np.isfinite(poly_x) & np.isfinite(poly_y)
        poly_x = np.asarray(poly_x[finite], dtype=np.float64)
        poly_y = np.asarray(poly_y[finite], dtype=np.float64)
        if poly_x.size > 0:
            poly_x_parts.append(poly_x)
            poly_y_parts.append(poly_y)
        selected_tokens.update(
            record.token
            for record in records
            if _tile_intersects_polygon(record, poly_x, poly_y)
        )
    selected_tokens = _expand_selection(selected_tokens, records, margin_tiles)
    selected = [record for record in records if record.token in selected_tokens]
    selected.sort(key=lambda item: item.token)
    if poly_x_parts:
        combined_x = np.concatenate(poly_x_parts)
        combined_y = np.concatenate(poly_y_parts)
    else:
        combined_x = np.array([], dtype=np.float64)
        combined_y = np.array([], dtype=np.float64)
    return selected, combined_x, combined_y


def tile_distance_km(record: TileRecord, meta: TemplateMeta, *, observer_lat: float, observer_lon: float) -> float:
    to_proj = Transformer.from_crs("EPSG:4326", meta.crs, always_xy=True)
    obs_x, obs_y = to_proj.transform(observer_lon, observer_lat)
    center_x = (record.x_min + record.x_max) * 0.5
    center_y = (record.y_min + record.y_max) * 0.5
    return float(math.hypot(center_x - float(obs_x), center_y - float(obs_y)) / 1000.0)


def stitch_tiles_from_paths(paths: Sequence[Path], *, source_label: str | None = None) -> xr.Dataset:
    if not paths:
        raise ValueError("No tile paths provided")
    with xr.open_dataset(paths[0]) as first:
        tile_height = int(first.sizes["y"])
        tile_width = int(first.sizes["x"])
        first_projection_attrs = dict(first[GRID_VAR].attrs)
        product_rows = int(first.attrs["product_rows"])
        product_cols = int(first.attrs["product_columns"])
        channel_id = int(first.attrs["channel_id"])
        x_step = float(first.x.values[1] - first.x.values[0])
        y_step = float(first.y.values[1] - first.y.values[0])
        x0 = float(first.x.values[0] - int(first.attrs["tile_column_offset"]) * x_step)
        y0 = float(first.y.values[0] - int(first.attrs["tile_row_offset"]) * y_step)

        full = np.full((product_rows, product_cols), np.nan, dtype=np.float32)
        coverage = np.zeros((product_rows, product_cols), dtype=np.uint8)

        for path in paths:
            with xr.open_dataset(path) as ds:
                row0 = int(ds.attrs["tile_row_offset"])
                col0 = int(ds.attrs["tile_column_offset"])
                row1 = row0 + tile_height
                col1 = col0 + tile_width
                tile = np.asarray(ds[DATA_VAR].values, dtype=np.float32)
                if full[row0:row1, col0:col1].shape != tile.shape:
                    raise ValueError(f"Tile shape mismatch for {path.name}")
                if np.any(coverage[row0:row1, col0:col1]):
                    raise ValueError(f"Tile overlap detected for {path.name}")
                expected_x0 = x0 + col0 * x_step
                expected_y0 = y0 + row0 * y_step
                if not math.isclose(float(ds.x.values[0]), expected_x0, rel_tol=0.0, abs_tol=1e-6):
                    raise ValueError(f"x coordinate mismatch for {path.name}")
                if not math.isclose(float(ds.y.values[0]), expected_y0, rel_tol=0.0, abs_tol=1e-6):
                    raise ValueError(f"y coordinate mismatch for {path.name}")
                full[row0:row1, col0:col1] = tile
                coverage[row0:row1, col0:col1] = 1

        x = x0 + np.arange(product_cols, dtype=np.float64) * x_step
        y = y0 + np.arange(product_rows, dtype=np.float64) * y_step
        attrs = dict(first[DATA_VAR].attrs)
        attrs.update(
            {
                "source_tile_count": len(paths),
                "source_tile_dir": source_label or str(paths[0].parent),
                "channel_id": channel_id,
                "coverage_fraction": float(coverage.mean()),
            }
        )
        data = xr.DataArray(full, dims=("y", "x"), coords={"y": y, "x": x}, name=DATA_VAR, attrs=attrs)
        projection = xr.DataArray(0, name=GRID_VAR, attrs=first_projection_attrs)
        merged = xr.Dataset(
            data_vars={DATA_VAR: data, GRID_VAR: projection, "tile_coverage_mask": (("y", "x"), coverage)},
            coords={"y": y, "x": x},
            attrs=dict(first.attrs),
        )
        merged.attrs.update(
            {
                "stitched_from_tiles": len(paths),
                "stitched_tile_height": tile_height,
                "stitched_tile_width": tile_width,
            }
        )
        return merged


def build_area_from_dataset(ds: xr.Dataset) -> GeoArea:
    proj_var = ds[GRID_VAR]
    crs = CRS.from_cf(dict(proj_var.attrs))
    x = np.asarray(ds.x.values, dtype=np.float64)
    y = np.asarray(ds.y.values, dtype=np.float64)
    geos_scale = float(proj_var.attrs["perspective_point_height"]) * 1e-6
    x_step = float(x[1] - x[0])
    y_step = float(y[1] - y[0])
    # ISatSS x/y coordinates are scan-angle-like values in microradian units.
    # Convert them to the projected units expected by the projection helpers.
    area_extent = (
        float(x[0] * geos_scale),
        float((y[0] + y.size * y_step) * geos_scale),
        float((x[0] + x.size * x_step) * geos_scale),
        float(y[0] * geos_scale),
    )
    return GeoArea(
        area_id="himawari_isatss_m1c13",
        description="Himawari ISatSS stitched M1C13",
        proj_id="geos",
        projection=crs,
        width=int(x.size),
        height=int(y.size),
        area_extent=area_extent,
    )


def attach_area_to_dataset(ds: xr.Dataset) -> xr.Dataset:
    area = build_area_from_dataset(ds)
    out = ds.copy(deep=False)
    da = out[DATA_VAR].copy(deep=False)
    attrs = dict(da.attrs)
    attrs["area"] = area
    da.attrs = attrs
    out[DATA_VAR] = da
    out.attrs = dict(out.attrs)
    out.attrs["area_id"] = area.area_id
    return out
