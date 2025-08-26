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


class CloudDiscError(RuntimeError):
    """
    Base class for all clouddisc exceptions.
    Carries optional CloudMeta-like context for uniform handling.
    """
    def __init__(self, message: str = "", *, meta: "CloudMeta | None" = None):
        super().__init__(message)
        self.meta = meta


class DataNotFoundError(CloudDiscError):
    """Raised when satellite data is not found for the requested time."""

    ...
    # Inherits meta from CloudDiscError


class DownloadError(CloudDiscError):
    """Raised when an error occurs while downloading satellite data."""

    def __init__(
        self,
        message: str = "",
        *,
        transient: bool = True,
        meta: "CloudMeta | None" = None
    ):
        super().__init__(message, meta=meta)
        self.transient = transient


class RenderError(RuntimeError):
    """Raised when an error occurs during the rendering process."""

    ...
    # Inherits meta from CloudDiscError
