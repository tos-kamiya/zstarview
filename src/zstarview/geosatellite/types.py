from __future__ import annotations

import datetime as dt
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from ..clouddisc.altaz_grid import CloudAltAzGrid

GeoSatelliteKind = Literal["infrared", "visible"]


@dataclass(frozen=True, slots=True)
class GeoSatelliteDownloadResult:
    """Metadata and payload for a downloaded Geo-satellite frame."""

    fetched_at_utc: dt.datetime
    captured_at_utc: dt.datetime | None
    kind: GeoSatelliteKind
    source_url: str
    png_bytes: bytes
    content_type: str | None = None
    cache_path: Path | None = None
    metadata_path: Path | None = None


@dataclass(frozen=True, slots=True)
class GeoSatelliteIntermediateResult:
    """Intermediate images derived from a downloaded frame."""

    raw_digest: str
    proxy_gray: np.ndarray
    inpainted_gray: np.ndarray
    manifest: dict[str, object] | None = None
    proxy_path: Path | None = None
    inpainted_path: Path | None = None
    mask_path: Path | None = None


@dataclass(frozen=True, slots=True)
class GeoSatelliteVoxelSource:
    """Display-derived cloud amounts on a fixed georeferenced image grid."""

    cloud_amount: np.ndarray
    valid_mask: np.ndarray
    grid_npz: Path
    kind: GeoSatelliteKind
    time_utc: dt.datetime
    source_key: str


@dataclass(frozen=True, slots=True)
class GeoSatellitePipelineResult:
    """Full in-memory result of the experimental Geo-satellite workflow."""

    download: GeoSatelliteDownloadResult
    intermediate: GeoSatelliteIntermediateResult
    disc_gray: np.ndarray
    altaz_grid: CloudAltAzGrid
    voxel_source: GeoSatelliteVoxelSource | None = None
