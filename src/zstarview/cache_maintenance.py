from __future__ import annotations

import logging
import shutil
from pathlib import Path

from .paths import CACHE_PATH, OVERTURE_DERIVED_ROOT_DIR, OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR


logger = logging.getLogger(__name__)


def clear_long_lived_cache() -> tuple[Path, ...]:
    removed: list[Path] = []
    targets = (
        Path(CACHE_PATH) / "copernicus-dem",
        Path(OVERTURE_DERIVED_ROOT_DIR),
        Path(OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR),
    )
    for path in targets:
        if not path.exists():
            logger.info("Long-lived cache already absent: %s", path)
            continue
        shutil.rmtree(path)
        removed.append(path)
        logger.info("Cleared long-lived cache: %s", path)
    return tuple(removed)
