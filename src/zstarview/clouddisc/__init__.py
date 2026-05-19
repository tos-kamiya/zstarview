# SPDX-FileCopyrightText: 2025-present Toshihiro Kamiya <kamiya@mbj.nifty.com>
#
# SPDX-License-Identifier: MIT
"""
CloudDisc package initializer.

This module exposes the main classes and exceptions for the clouddisc package.
"""
from .cache.cleanup import cleanup_satellite_cache
from .config import CloudDiscConfig
from .core import CloudDisc
from .types import (
    CloudDiscError,
    CloudMeta,
    CloudSourceData,
    DataNotFoundError,
    DownloadCancelledError,
    DownloadError,
    RenderError,
    RenderKey,
    SourceKey,
    TimeoutError,
    VisibilityError,
    round_down_utc_to_slot,
)

__all__ = [
    "CloudDisc",
    "CloudDiscConfig",
    "CloudMeta",
    "CloudSourceData",
    "SourceKey",
    "RenderKey",
    "round_down_utc_to_slot",
    "VisibilityError",
    "CloudDiscError",
    "DataNotFoundError",
    "DownloadError",
    "DownloadCancelledError",
    "TimeoutError",
    "RenderError",
    "cleanup_satellite_cache",
]
