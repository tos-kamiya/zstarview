import logging
import os
import sys
from pathlib import Path

from .paths import LOG_PATH

logger = logging.getLogger(__name__)


def _resolve_log_level() -> int:
    raw = os.getenv("LOG_LEVEL", "").strip()
    if not raw:
        return int(logging.INFO)
    if raw.isdigit():
        return int(raw)
    level = logging.getLevelName(raw.upper())
    if isinstance(level, int):
        return int(level)
    return int(logging.INFO)


def setup_root_logger() -> logging.Logger:
    """Configure and return the root logger for the application."""
    log_level = _resolve_log_level()
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        stream=sys.stderr,
    )
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)

    log_dir = Path(LOG_PATH)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "app.log"

    try:
        file_handler = logging.FileHandler(log_path, encoding="utf-8")
    except OSError as exc:
        logger.warning("Failed to open log file %s: %s", log_path, exc)
        return root_logger

    file_handler.setLevel(log_level)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
    )
    root_logger.addHandler(file_handler)

    logger.info("Logging to file: %s", log_path)
    return root_logger
