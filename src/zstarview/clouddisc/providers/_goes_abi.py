# -*- coding: utf-8 -*-
"""Helpers for direct GOES ABI CMI ingestion without Satpy."""

from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np
import xarray as xr
from pyproj import CRS

from ..geo_area import GeoArea


DATA_VAR = "CMI"
GRID_VAR = "goes_imager_projection"

warnings.filterwarnings(
    "ignore",
    message="You will likely lose important projection information when converting to a PROJ string",
    category=UserWarning,
    module="pyproj.crs.crs",
)


def build_area_from_cmi_dataset(ds: xr.Dataset) -> GeoArea:
    """Build a geostationary area definition from a GOES CMI dataset."""
    proj_var = ds[GRID_VAR]
    crs = CRS.from_cf(dict(proj_var.attrs))
    scale = float(proj_var.attrs["perspective_point_height"])
    x = np.asarray(ds.x.values, dtype=np.float64)
    y = np.asarray(ds.y.values, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or x.size < 2 or y.size < 2:
        raise ValueError("Expected 1-D x/y coordinates with at least 2 points each")
    x_step = float(x[1] - x[0])
    y_step = float(y[1] - y[0])
    area_extent = (
        float(x[0] * scale),
        float((y[0] + y.size * y_step) * scale),
        float((x[0] + x.size * x_step) * scale),
        float(y[0] * scale),
    )
    return GeoArea(
        area_id="goes_abi_cmi_c13",
        description="GOES ABI CMIPF C13",
        proj_id="geos",
        projection=crs,
        width=int(x.size),
        height=int(y.size),
        area_extent=area_extent,
    )


def load_cmi_with_area(path: Path) -> xr.DataArray:
    """Load a GOES CMI file and attach Satpy-like ``area`` metadata."""
    with xr.open_dataset(path) as ds:
        if DATA_VAR not in ds.variables:
            raise ValueError(f"{path.name} does not contain {DATA_VAR}")
        if GRID_VAR not in ds.variables:
            raise ValueError(f"{path.name} does not contain {GRID_VAR}")
        da = ds[DATA_VAR].astype(np.float32).compute()
        area = build_area_from_cmi_dataset(ds)

    da.attrs = dict(da.attrs)
    da.attrs["area"] = area
    return da
