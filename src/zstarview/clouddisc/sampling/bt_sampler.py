# -*- coding: utf-8 -*-
"""
Builds a sampler function for resampling geostationary satellite data.

This module provides a factory function that creates a high-performance sampler.
The sampler is used to take data on a satellite's native projection (like GOES
or Himawari) and find the corresponding data values for a new set of arbitrary
longitude/latitude points. This is essential for re-projecting the satellite
view to match the observer's perspective.
"""

import warnings
from typing import Callable

import numpy as np
import xarray as xr
from pyproj import CRS, Transformer

# Suppress a UserWarning from pyproj that is not relevant to this use case.
warnings.filterwarnings(
    "ignore",
    message="You will likely lose important projection information when converting to a PROJ string",
    category=UserWarning,
    module="pyproj.crs.crs",
)


def build_bt_sampler(da: xr.DataArray) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    """
    Builds a vectorized, NaN-aware bilinear interpolation sampler for a geostationary grid.

    This is a factory function: it takes a DataArray and returns a new `sampler`
    function. This returned function is computationally expensive but highly optimized.
    It can be called with longitude and latitude arrays to get interpolated brightness
    temperatures from the original satellite data grid.

    Args:
        da: The input xarray DataArray containing brightness temperatures. It is
            expected to have an `area` attribute that defines the
            geostationary projection.

    Returns:
        A sampler function. This function takes (longitude, latitude) numpy arrays
        and returns a numpy array of the corresponding interpolated brightness temperatures.
    """
    area = da.attrs.get("area")
    if area is None:
        raise ValueError("Input DataArray must have an 'area' attribute.")

    # --- 1. Extract projection and grid information from the area definition ---
    min_x, min_y, max_x, max_y = area.area_extent
    height, width = area.shape
    pixel_size_x = (max_x - min_x) / width
    pixel_size_y = (max_y - min_y) / height
    grid_data = np.asarray(da.compute().values, dtype=np.float32)

    # --- 2. Create the coordinate transformer ---
    # This transformer will convert from standard lon/lat (EPSG:4326) to the
    # satellite's native projection system (e.g., geostationary).
    target_crs = getattr(area, "crs", CRS.from_dict(area.proj_dict))
    transformer = Transformer.from_crs("EPSG:4326", target_crs, always_xy=True)

    # Reuse float index buffers across calls when input shape is unchanged.
    ix_f_buf: np.ndarray | None = None
    iy_f_buf: np.ndarray | None = None
    buf_shape: tuple[int, ...] | None = None

    # --- 3. Define and return the sampler function ---
    def sampler(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
        """
        Performs bilinear interpolation on the source grid for the given lon/lat points.
        """
        nonlocal ix_f_buf, iy_f_buf, buf_shape
        if lon.shape != lat.shape:
            raise ValueError("lon and lat must have the same shape")

        if buf_shape != lon.shape or ix_f_buf is None or iy_f_buf is None:
            ix_f_buf = np.empty(lon.shape, dtype=np.float32)
            iy_f_buf = np.empty(lon.shape, dtype=np.float32)
            buf_shape = lon.shape

        # Project lon/lat coordinates to the satellite's grid coordinates (x, y).
        projected_x, projected_y = transformer.transform(lon, lat)

        # Convert projected (x, y) to floating-point pixel indices (ix, iy).
        # Points that cannot be projected (e.g., on the other side of the Earth)
        # will result in non-finite values, which we handle with a mask.
        finite_mask = np.isfinite(projected_x) & np.isfinite(projected_y)
        ix_f = ix_f_buf
        iy_f = iy_f_buf
        ix_f.fill(np.nan)
        iy_f.fill(np.nan)
        ix_f[finite_mask] = (projected_x[finite_mask] - min_x) / pixel_size_x
        iy_f[finite_mask] = (max_y - projected_y[finite_mask]) / pixel_size_y  # Y-axis is inverted

        # --- NaN-aware Bilinear Interpolation ---
        output_bt = np.full(lon.shape, np.nan, dtype=np.float32)

        # Create a mask for points that fall within the grid boundaries, where a 2x2
        # box of pixels can be formed for interpolation.
        in_bounds_mask = (ix_f >= 0.0) & (ix_f < (width - 1)) & (iy_f >= 0.0) & (iy_f < (height - 1))
        valid_mask = finite_mask & in_bounds_mask
        if not np.any(valid_mask):
            return output_bt

        # Get coordinates and interpolation weights for the valid points.
        ix_valid, iy_valid = ix_f[valid_mask], iy_f[valid_mask]
        ix0, iy0 = np.floor(ix_valid).astype(np.int32), np.floor(iy_valid).astype(np.int32)
        ix1, iy1 = ix0 + 1, iy0 + 1
        wx, wy = ix_valid - ix0, iy_valid - iy0  # Weights are the fractional part

        # Get the values of the four surrounding pixels.
        v00, v10 = grid_data[iy0, ix0], grid_data[iy0, ix1]
        v01, v11 = grid_data[iy1, ix0], grid_data[iy1, ix1]

        # Handle NaNs in the source data by adjusting interpolation weights.
        # If a neighboring pixel is NaN, its weight becomes zero.
        m00, m10 = np.isfinite(v00), np.isfinite(v10)
        m01, m11 = np.isfinite(v01), np.isfinite(v11)
        w00, w10 = (1 - wx) * (1 - wy) * m00, wx * (1 - wy) * m10
        w01, w11 = (1 - wx) * wy * m01, wx * wy * m11

        # Perform the weighted sum of non-NaN neighbors.
        total_weight = w00 + w10 + w01 + w11
        interp_values = np.where(m00, v00, 0.0) * w00 + np.where(m10, v10, 0.0) * w10 + np.where(m01, v01, 0.0) * w01 + np.where(m11, v11, 0.0) * w11

        # Normalize the result and handle cases where all neighbors were NaN.
        good_weights = total_weight > 1e-6
        interp_values[good_weights] /= total_weight[good_weights]
        interp_values[~good_weights] = np.nan

        # Place the interpolated values into the final output array.
        output_bt[valid_mask] = interp_values
        return output_bt

    return sampler
