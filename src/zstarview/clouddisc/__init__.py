# SPDX-FileCopyrightText: 2025-present Toshihiro Kamiya <kamiya@mbj.nifty.com>
#
# SPDX-License-Identifier: MIT
"""
CloudDisc package initializer.

This module exposes the main classes and exceptions for the clouddisc package.
"""
from .config import CloudDiscConfig
from .types import CloudMeta, CloudSourceData, SourceKey, RenderKey, round_down_utc_to_slot
from .types import VisibilityError, CloudDiscError, DataNotFoundError, DownloadError, DownloadCancelledError, TimeoutError, RenderError
from .core import CloudDisc
from .cache.cleanup import cleanup_satellite_cache

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
