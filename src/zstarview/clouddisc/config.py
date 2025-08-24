from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Optional
import os

@dataclass(frozen=True)
class CloudDiscConfig:
    """ライブラリ全体設定"""
    cache_dir: Optional[Path] = None
    sat_priority: Tuple[str, ...] = ("AUTO",)  # or ("HIMAWARI","G18","G16")
    bt_warm_k: float = 310.0
    bt_cold_k: float = 190.0
    gamma: float = 1.6
    alt_min_deg: float = -2.0
    search_back_minutes: int = 120
    edge_antialias: bool = False

    def cache_root(self) -> Path:
        if self.cache_dir:
            root = Path(self.cache_dir).expanduser()
        else:
            # デフォルト: ~/.cache/zstarview/clouddisc
            root = Path(os.path.expanduser("~")) / ".cache" / "zstarview" / "clouddisc"
        root.mkdir(parents=True, exist_ok=True)
        return root
