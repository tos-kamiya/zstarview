"""
Defines custom data structures and exceptions for the CloudDisc library.

This module contains the `CloudMeta` dataclass for holding metadata about
rendered cloud images and a set of custom exception classes for handling
specific errors that can occur during data fetching and processing.
"""

import datetime as dt
from dataclasses import dataclass
from pathlib import Path
from typing import Any


def round_down_utc_to_slot(
    when: dt.datetime,
    *,
    slot_minutes: int = 10,
) -> dt.datetime:
    """Round a datetime down to a UTC slot boundary."""
    if slot_minutes <= 0:
        raise ValueError("slot_minutes must be > 0")

    t = when
    if t.tzinfo is None:
        t = t.replace(tzinfo=dt.timezone.utc)
    else:
        t = t.astimezone(dt.timezone.utc)

    minute = (t.minute // slot_minutes) * slot_minutes
    return t.replace(minute=minute, second=0, microsecond=0)


@dataclass(frozen=True, slots=True)
class SourceKey:
    """Identity for cloud source data acquisition/cache lookup."""

    satellite: str
    timeslot_utc: dt.datetime
    provider: str | None = None
    sat_priority: tuple[str, ...] = ("AUTO",)

    def __post_init__(self) -> None:
        object.__setattr__(self, "timeslot_utc", round_down_utc_to_slot(self.timeslot_utc))
        object.__setattr__(self, "satellite", str(self.satellite))
        if self.provider is not None:
            object.__setattr__(self, "provider", str(self.provider))
        object.__setattr__(self, "sat_priority", tuple(self.sat_priority))


@dataclass(frozen=True, slots=True)
class RenderKey:
    """Identity for camera-dependent cloud rendering."""

    source: SourceKey
    alt_deg: float
    az_deg: float
    radius_px: int
    edge_fov_deg: float = 90.0
    mask_fov_deg: float = 93.0

    def __post_init__(self) -> None:
        az = float(self.az_deg) % 360.0
        alt = float(self.alt_deg)
        alt = max(-89.999, min(89.999, alt))
        object.__setattr__(self, "az_deg", az)
        object.__setattr__(self, "alt_deg", alt)
        object.__setattr__(self, "radius_px", max(1, int(self.radius_px)))
        object.__setattr__(self, "edge_fov_deg", float(self.edge_fov_deg))
        object.__setattr__(self, "mask_fov_deg", float(self.mask_fov_deg))


@dataclass(slots=True)
class CloudSourceData:
    """Fetched cloud source data before camera-dependent rendering."""

    source_key: SourceKey
    data_array: Any
    satellite: str
    product: str
    time_utc: dt.datetime
    src_paths: list[Path]
    sampler: Any = None
    source_expected_count: int | None = None
    source_available_count: int | None = None
    source_completeness_ratio: float | None = None
    altaz_grid: Any = None


@dataclass
class CloudMeta:
    """
    Metadata for a rendered cloud image.

    This structure holds information about the source data used to generate
    a particular cloud image.

    Attributes:
        satellite: The name of the satellite used, e.g., "G19", "G18", "HIMAWARI".
        product: The data product identifier, e.g., "CMIPF-C13", "HSD-B13".
        time_utc: The UTC timestamp of the source satellite data.
        src_paths: A list of local file paths to the source data files.
    """

    satellite: str
    product: str
    time_utc: dt.datetime
    src_paths: list[Path]


# --- Custom Exceptions ---


class CloudDiscError(RuntimeError):
    """
    Base class for all custom exceptions in the clouddisc library.

    This allows for catching all library-specific errors with a single `except` block.
    It can optionally carry metadata about the operation that failed.
    """

    def __init__(self, message: str = "", *, meta: CloudMeta | None = None):
        super().__init__(message)
        self.meta = meta


class VisibilityError(CloudDiscError):
    """Raised when no suitable satellite is visible from the given location."""



class DataNotFoundError(CloudDiscError):
    """Raised when satellite data is not found for the requested time, even after searching back."""



class DownloadError(CloudDiscError):
    """Raised when an error occurs while downloading satellite data."""

    def __init__(self, message: str = "", *, meta: CloudMeta | None = None):
        super().__init__(message, meta=meta)


class DownloadCancelledError(DownloadError):
    """Raised when a download is cancelled cooperatively."""



class TimeoutError(DownloadError):
    """Raised when a download operation times out."""



class RenderError(CloudDiscError):
    """Raised when an error occurs during the image rendering process."""
