"""Shared defaults for cloud worker tasks."""

from __future__ import annotations

DEFAULT_CLOUD_SHELLS_KM: tuple[float, ...] = tuple(6371.0 + height for height in range(3, 12))
