"""
Configuration for the CloudDisc library.

This module defines the `CloudDiscConfig` dataclass, which holds all the
settings for fetching, caching, and processing satellite cloud data.
"""

import os
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class CloudDiscConfig:
    """
    Global configuration for the clouddisc library.

    This object is used to control various aspects of the cloud data processing,
    from data source priority to rendering thresholds.

    Attributes:
        cache_dir: The directory to cache satellite data. If None, a default path is used.
        sat_priority: The priority order for satellite selection. "AUTO" enables automatic
                      selection based on the observer's longitude.
                      Example: ("HIMAWARI", "G18", "G19")
        bt_warm_k: The brightness temperature (in Kelvin) considered "warm" for rendering.
                   This corresponds to lower, warmer clouds or the ground.
        bt_cold_k: The brightness temperature (in Kelvin) considered "cold" for rendering.
                   This corresponds to high-altitude, cold cloud tops.
        alt_min_deg: The minimum altitude (in degrees) for a satellite to be considered
                     visible and thus a viable data source.
        search_back_minutes: The time window in minutes to search backward for satellite
                             data if the most recent data is not available.
        connect_timeout, read_timeout: Timeout values for downloading IR data from satellite.
    """

    cache_dir: Path | None = None
    sat_priority: tuple[str, ...] = field(default_factory=lambda: ("AUTO",))
    bt_warm_k: float = 310.0
    bt_cold_k: float = 190.0
    alt_min_deg: float = 0.0
    search_back_minutes: int = 120
    connect_timeout: float = 5.0
    read_timeout: float = 30.0

    def cache_root(self) -> Path:
        """
        Gets the root directory for the cache.

        If `cache_dir` is explicitly set, it is used. Otherwise, a default directory
        is created at `~/.cache/zstarview/clouddisc` on Linux/macOS.
        The directory is created if it does not exist.

        Returns:
            The path to the cache root directory.
        """
        if self.cache_dir:
            root = Path(self.cache_dir).expanduser()
        else:
            # Define a standard, user-specific cache location to avoid cluttering the project directory.
            default_cache_path = (
                Path(os.path.expanduser("~")) / ".cache" / "zstarview" / "clouddisc"
            )
            root = default_cache_path

        # Ensure the directory exists before returning it.
        root.mkdir(parents=True, exist_ok=True)
        return root
