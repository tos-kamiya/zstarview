# -*- coding: utf-8 -*-
"""
Defines custom data structures and exceptions for the CloudDisc library.

This module contains the `CloudMeta` dataclass for holding metadata about
rendered cloud images and a set of custom exception classes for handling
specific errors that can occur during data fetching and processing.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional
import datetime as dt


@dataclass
class CloudMeta:
    """
    Metadata for a rendered cloud image.

    This structure holds information about the source data used to generate
    a particular cloud image.

    Attributes:
        satellite: The name of the satellite used, e.g., "G16", "G18", "HIMAWARI".
        product: The data product identifier, e.g., "CMIPF-C13", "HSD-B13".
        time_utc: The UTC timestamp of the source satellite data.
        src_paths: A list of local file paths to the source data files.
    """

    satellite: str
    product: str
    time_utc: dt.datetime
    src_paths: List[Path]


# --- Custom Exceptions ---


class CloudDiscError(RuntimeError):
    """
    Base class for all custom exceptions in the clouddisc library.

    This allows for catching all library-specific errors with a single `except` block.
    It can optionally carry metadata about the operation that failed.
    """

    def __init__(self, message: str = "", *, meta: Optional[CloudMeta] = None):
        super().__init__(message)
        self.meta = meta


class VisibilityError(CloudDiscError):
    """Raised when no suitable satellite is visible from the given location."""

    pass


class DataNotFoundError(CloudDiscError):
    """Raised when satellite data is not found for the requested time, even after searching back."""

    pass


class DownloadError(CloudDiscError):
    """Raised when an error occurs while downloading satellite data."""

    def __init__(self, message: str = "", *, meta: Optional[CloudMeta] = None):
        super().__init__(message, meta=meta)


class TimeoutError(DownloadError):
    """Raised when a download operation times out."""

    pass


class RenderError(CloudDiscError):
    """Raised when an error occurs during the image rendering process."""

    pass

