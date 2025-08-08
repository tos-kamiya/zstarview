import json
from pathlib import Path

from appdirs import user_config_dir

from .paths import APP_ID, APP_AUTHOR


_config_dir = Path(user_config_dir(APP_ID, APP_AUTHOR))
_config_dir.mkdir(parents=True, exist_ok=True)
_config_file = _config_dir / "config.json"


def load_last_city() -> str | None:
    try:
        data = json.loads(_config_file.read_text(encoding="utf-8"))
        return data.get("city")
    except FileNotFoundError:
        return None
    except json.JSONDecodeError:
        return None


def save_last_city(city: str) -> None:
    _config_file.write_text(json.dumps({"city": city}, ensure_ascii=False), encoding="utf-8")
