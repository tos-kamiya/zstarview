# -*- coding: utf-8 -*-
"""Helpers for direct GOES ABI CMI ingestion without Satpy."""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import numpy as np
import xarray as xr
from pyproj import CRS

from ..diagnostics import DiagnosticSink, emit_diagnostic
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


def load_cmi_with_area(path: Path, *, diagnostic_sink: DiagnosticSink | None = None) -> xr.DataArray:
    """Load a GOES CMI file and attach Satpy-like ``area`` metadata."""
    logger = logging.getLogger(__name__)
    logger.debug("GOES CMI open start: %s", path)
    emit_diagnostic(
        diagnostic_sink,
        "open_source_file",
        "start",
        "Opening GOES CMI NetCDF file",
        path=path,
        decode_times=False,
    )
    with xr.open_dataset(path, decode_times=False) as ds:
        logger.debug("GOES CMI dataset opened: %s", path)
        emit_diagnostic(
            diagnostic_sink,
            "open_source_file",
            "ok",
            "GOES CMI NetCDF file opened",
            path=path,
        )
        if DATA_VAR not in ds.variables:
            raise ValueError(f"{path.name} does not contain {DATA_VAR}")
        if GRID_VAR not in ds.variables:
            raise ValueError(f"{path.name} does not contain {GRID_VAR}")
        logger.debug("GOES CMI compute start: %s", path)
        emit_diagnostic(
            diagnostic_sink,
            "load_brightness_temperature",
            "start",
            "Loading GOES CMI brightness temperature",
            path=path,
        )
        da = ds[DATA_VAR].astype(np.float32).compute()
        logger.debug("GOES CMI compute done: %s", path)
        emit_diagnostic(
            diagnostic_sink,
            "load_brightness_temperature",
            "ok",
            "GOES CMI brightness temperature loaded",
            path=path,
            shape=tuple(int(v) for v in da.shape),
            dtype=str(da.dtype),
        )
        logger.debug("GOES CMI area build start: %s", path)
        emit_diagnostic(
            diagnostic_sink,
            "build_projection_area",
            "start",
            "Building GOES projection area",
            path=path,
        )
        area = build_area_from_cmi_dataset(ds)
        logger.debug("GOES CMI area build done: %s", path)
        emit_diagnostic(
            diagnostic_sink,
            "build_projection_area",
            "ok",
            "GOES projection area built",
            path=path,
            width=area.width,
            height=area.height,
        )

    da.attrs = dict(da.attrs)
    da.attrs["area"] = area
    return da
