import json
import logging
from pathlib import Path
from typing import Any

from appdirs import user_config_dir

from .paths import APP_AUTHOR, APP_ID

logger = logging.getLogger(__name__)


_config_dir = Path(user_config_dir(APP_ID, APP_AUTHOR))
_config_dir.mkdir(parents=True, exist_ok=True)
_config_file = _config_dir / "config.json"
OPEN_METEO_NONCOMMERCIAL_TERMS_VERSION = 1


def _load_config() -> dict[str, Any]:
    try:
        data = json.loads(_config_file.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def _save_config(data: dict[str, Any]) -> None:
    _config_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        _config_file.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
    except OSError as exc:
        logger.warning("Failed to save config file %s: %s", _config_file, exc)


def load_last_city() -> str | dict[str, Any] | None:
    """Load the last used location payload from the config file."""
    city = _load_config().get("city")
    if isinstance(city, str):
        return city
    return city if isinstance(city, dict) else None


def save_last_city(city: str | dict[str, Any]) -> None:
    """Save the last used location payload to the config file."""
    data = _load_config()
    data["city"] = city
    _save_config(data)


def load_last_window_geometry() -> tuple[int, int, int, int] | None:
    """Load the last saved window geometry as (x, y, width, height)."""
    geom = _load_config().get("window_geometry")
    if not isinstance(geom, dict):
        return None
    try:
        x = int(geom["x"])
        y = int(geom["y"])
        width = int(geom["width"])
        height = int(geom["height"])
    except (KeyError, TypeError, ValueError):
        return None
    return (x, y, width, height)


def save_last_window_geometry(x: int, y: int, width: int, height: int) -> None:
    """Save window geometry."""
    data = _load_config()
    data["window_geometry"] = {
        "x": int(x),
        "y": int(y),
        "width": int(width),
        "height": int(height),
    }
    _save_config(data)


def open_meteo_noncommercial_terms_accepted() -> bool:
    value = _load_config().get("open_meteo_noncommercial_terms_version")
    return value == OPEN_METEO_NONCOMMERCIAL_TERMS_VERSION


def accept_open_meteo_noncommercial_terms() -> None:
    data = _load_config()
    data["open_meteo_noncommercial_terms_version"] = (
        OPEN_METEO_NONCOMMERCIAL_TERMS_VERSION
    )
    _save_config(data)
