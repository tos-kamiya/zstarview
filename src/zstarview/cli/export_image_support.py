"""Shared constants, types, and helpers for export-image modules."""

from __future__ import annotations

import logging
import sys
import time
from collections import Counter
from dataclasses import dataclass
from typing import TypedDict

import numpy as np

from ..terrain import DEFAULT_TERRAIN_DISTANCE_SAMPLE_STEP_M
from ..types import UrbanOutlinePolyline
from ..water_mask_interface import WaterSurfaceBandStats
from ..water_overlay import WaterOverlayPoint, sample_water_overlay_points

LOGGER_NAME = "zstarview.cli.export_image"
logger = logging.getLogger(LOGGER_NAME)

EXPORT_IMAGE_METADATA_SCHEMA = "zstarview.export-image-metadata.v1"
EXPORT_IMAGE_METADATA_TEXT_KEY = "zstarview.export-image-metadata"
OPEN_METEO_CONSENT_REQUIRED_MESSAGE = (
    "Open-Meteo Free API terms have not been accepted.\n"
    "Start zstarview or zstarview-gui to review and accept the terms before "
    "using precipitation export."
)

DEFAULT_CLOUD_ALT_MIN_DEG = 1.0
DEFAULT_CLOUD_FOV_OVERSCAN_DEG = 2.0
DEFAULT_CLOUD_BASE_SIZE = 256
DEFAULT_EXPORT_IMAGE_SKY_UPDATE_INTERVAL = 60
DEFAULT_EXPORT_IMAGE_TERRAIN_AZIMUTH_STEP_DEG = 0.5
DEFAULT_EXPORT_IMAGE_TERRAIN_SAMPLE_STEP_M = DEFAULT_TERRAIN_DISTANCE_SAMPLE_STEP_M

sample_water_overlay_points_for_observer = sample_water_overlay_points


def host():
    """Return the public export_image module so tests can patch names there."""
    return sys.modules["zstarview.cli.export_image"]


@dataclass(frozen=True)
class _UrbanOutlineFetchResult:
    outlines: list[UrbanOutlinePolyline] | None
    source: str | None

class TerrainHorizonPayload(TypedDict):
    profile_altaz: list[tuple[float, float]]
    profile_distances_m: list[float]
    secondary_ridges_altaz_layers: list[list[tuple[float, float]]]
    secondary_ridges_distances_m_layers: list[list[float]]
    sample_distances_m: np.ndarray
    sample_terrain_elevation_m: np.ndarray

def _deadline_after(timeout_seconds: float) -> float | None:
    if float(timeout_seconds) <= 0.0:
        return None
    return time.monotonic() + float(timeout_seconds)

def _remaining_timeout_seconds(deadline: float | None) -> float | None:
    if deadline is None:
        return None
    return max(0.0, deadline - time.monotonic())

def _timed_out(deadline: float | None) -> bool:
    remaining = _remaining_timeout_seconds(deadline)
    return remaining is not None and remaining <= 0.0

def _water_overlay_band_stats_text(stats: WaterSurfaceBandStats) -> str:
    return (
        f"{stats.band_name} tiles={int(stats.loaded_tile_count)} "
        f"raw={int(stats.raw_point_count)} "
        f"collapsed={int(stats.collapsed_point_count)} "
        f"visible={int(stats.visible_point_count)}"
    )

def _water_overlay_band_counts(
    points: list[WaterOverlayPoint] | tuple[WaterOverlayPoint, ...],
) -> tuple[int, int, int]:
    counts = Counter(str(point.water_category).strip().lower() for point in points)
    return (
        int(counts.get("sea-125", 0)),
        int(counts.get("sea-250", 0)),
        int(counts.get("sea-500", 0)),
    )
