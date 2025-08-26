"""
Configuration for the CloudDisc library.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Tuple, Optional
import os


@dataclass(frozen=True)
class CloudDiscConfig:
    """
    Global configuration for the clouddisc library.

    Attributes:
        cache_dir: The directory to cache satellite data. If None, a default path is used.
        sat_priority: The priority order for satellite selection. "AUTO" for automatic selection.
                      Example: ("HIMAWARI", "G18", "G16")
        bt_warm_k: The brightness temperature (in Kelvin) to be considered "warm".
        bt_cold_k: The brightness temperature (in Kelvin) to be considered "cold".
        alt_min_deg: The minimum altitude (in degrees) for a celestial object to be visible.
        search_back_minutes: How many minutes to search back for satellite data if the latest is not available.
    """

    cache_dir: Optional[Path] = None
    sat_priority: Tuple[str, ...] = field(default_factory=lambda: ("AUTO",))
    bt_warm_k: float = 310.0
    bt_cold_k: float = 190.0
    alt_min_deg: float = -2.0
    search_back_minutes: int = 120

    def cache_root(self) -> Path:
        """
        Gets the root directory for the cache.

        If `cache_dir` is set, it is used. Otherwise, a default directory
        `~/.cache/zstarview/clouddisc` is used.
        The directory is created if it does not exist.

        Returns:
            The path to the cache root directory.
        """
        if self.cache_dir:
            root = Path(self.cache_dir).expanduser()
        else:
            # Default: ~/.cache/zstarview/clouddisc
            root = Path(os.path.expanduser("~")) / ".cache" / "zstarview" / "clouddisc"
        root.mkdir(parents=True, exist_ok=True)
        return root
