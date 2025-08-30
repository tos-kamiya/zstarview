"""
Builds a sampler function for resampling geostationary satellite data.
"""

import warnings
from typing import Callable

import numpy as np
import xarray as xr
from pyproj import CRS, Transformer

# Suppress a UserWarning from pyproj when creating a CRS from a proj string,
# as the area definition from Satpy often triggers this. The information loss
# is acceptable for this use case.
warnings.filterwarnings(
    "ignore",
    message="You will likely lose important projection information when converting to a PROJ string",
    category=UserWarning,
    module="pyproj.crs.crs",
)


def build_bt_sampler(da: xr.DataArray) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    """
    Builds a vectorized sampler function for a geostationary grid.

    This function takes a DataArray of brightness temperatures on a fixed
    geostationary projection and returns a new function. This new function
    (the "sampler") can take arrays of longitude and latitude and return the
    corresponding brightness temperature for each point, using NaN-aware
    bilinear interpolation.

    Args:
        da: The input DataArray containing brightness temperatures. It is
            expected to have an `area` attribute provided by Satpy, which
            defines the geostationary projection.

    Returns:
        A sampler function that takes (longitude, latitude) numpy arrays
        and returns a numpy array of interpolated brightness temperatures.
    """
    area = da.attrs.get("area")
    if area is None:
        raise ValueError("Input DataArray must have an 'area' attribute from Satpy.")

    # --- 1. Extract projection and grid information from the area definition ---
    min_x, min_y, max_x, max_y = area.area_extent
    height, width = area.shape
    pixel_size_x = (max_x - min_x) / width
    pixel_size_y = (max_y - min_y) / height

    # Get the underlying numpy data
    grid_data = np.asarray(da.compute().values, dtype=np.float32)

    # Create the coordinate transformer from lon/lat (EPSG:4326) to the satellite's projection
    target_crs = getattr(area, "crs", CRS.from_dict(area.proj_dict))
    transformer = Transformer.from_crs("EPSG:4326", target_crs, always_xy=True)

    def sampler(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
        """
        Performs bilinear interpolation on the source grid for the given lon/lat points.
        """
        # --- 2. Project lon/lat coordinates to the satellite's grid coordinates (x, y) ---
        projected_x, projected_y = transformer.transform(lon, lat)

        # --- 3. Convert projected (x, y) to floating-point pixel indices (ix, iy) ---
        # First, identify points that are valid for projection
        finite_mask = np.isfinite(projected_x) & np.isfinite(projected_y)

        # Calculate floating point indices.
        # Using np.full allows us to handle non-finite points gracefully.
        ix_f = np.full(lon.shape, np.nan, dtype=np.float32)
        iy_f = np.full(lon.shape, np.nan, dtype=np.float32)

        ix_f[finite_mask] = (projected_x[finite_mask] - min_x) / pixel_size_x
        # The y-axis is inverted (row 0 is at the top, corresponding to max_y)
        iy_f[finite_mask] = (max_y - projected_y[finite_mask]) / pixel_size_y

        # --- 4. Perform NaN-aware bilinear interpolation ---
        # Initialize output array with NaNs
        output_bt = np.full(lon.shape, np.nan, dtype=np.float32)

        # Identify points that fall within the grid boundaries for interpolation
        # (i.e., where we can form a 2x2 box of pixels)
        in_bounds_mask = (ix_f >= 0.0) & (ix_f < (width - 1)) & (iy_f >= 0.0) & (iy_f < (height - 1))
        valid_mask = finite_mask & in_bounds_mask

        if not np.any(valid_mask):
            return output_bt  # Return all NaNs if no points are valid

        # Get coordinates and weights for valid points only
        ix_valid = ix_f[valid_mask]
        iy_valid = iy_f[valid_mask]

        # Get the integer indices of the top-left corner of the 4-pixel box
        ix0 = np.floor(ix_valid).astype(np.int32)
        iy0 = np.floor(iy_valid).astype(np.int32)
        ix1 = ix0 + 1
        iy1 = iy0 + 1

        # Get the fractional part of the indices, which are the interpolation weights
        wx = ix_valid - ix0
        wy = iy_valid - iy0

        # Get the values of the four surrounding pixels
        v00 = grid_data[iy0, ix0]
        v10 = grid_data[iy0, ix1]
        v01 = grid_data[iy1, ix0]
        v11 = grid_data[iy1, ix1]

        # --- 5. Handle NaNs in the source data ---
        # If a surrounding pixel is NaN, its contribution to the interpolation is 0.
        # We adjust the weights accordingly.
        m00 = np.isfinite(v00)
        m10 = np.isfinite(v10)
        m01 = np.isfinite(v01)
        m11 = np.isfinite(v11)

        w00 = (1 - wx) * (1 - wy) * m00
        w10 = wx * (1 - wy) * m10
        w01 = (1 - wx) * wy * m01
        w11 = wx * wy * m11

        total_weight = w00 + w10 + w01 + w11

        # Perform the interpolation
        interp_values = np.nan_to_num(v00) * w00 + np.nan_to_num(v10) * w10 + np.nan_to_num(v01) * w01 + np.nan_to_num(v11) * w11

        # Normalize the result by the total weight
        # If total_weight is 0 (all 4 neighbors were NaN), the result will be 0.
        # We replace these with NaN.
        good_weights = total_weight > 1e-6
        interp_values[good_weights] /= total_weight[good_weights]
        interp_values[~good_weights] = np.nan

        # Place the interpolated values into the output array
        output_bt[valid_mask] = interp_values

        return output_bt

    return sampler
