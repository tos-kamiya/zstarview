"""
Custom types and exceptions for the CloudDisc library.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List
import datetime as dt


@dataclass
class CloudMeta:
    """
    Metadata for a rendered cloud image.

    Attributes:
        satellite: The name of the satellite used, e.g., "G16", "G18", "HIMAWARI".
        product: The data product identifier, e.g., "CMIPF-C13", "HSD-B13", "ISatSS-B13".
        time_utc: The UTC timestamp of the data, rounded to 10 minutes.
        src_paths: A list of local file paths to the source data files.
    """

    satellite: str
    product: str
    time_utc: dt.datetime
    src_paths: List[Path]


class VisibilityError(RuntimeError):
    """Raised when no suitable satellite is visible from the given location."""

    ...


class DataNotFoundError(RuntimeError):
    """Raised when satellite data is not found for the requested time."""

    ...


class DownloadError(RuntimeError):
    """Raised when an error occurs while downloading satellite data."""

    ...


class RenderError(RuntimeError):
    """Raised when an error occurs during the rendering process."""

    ...
