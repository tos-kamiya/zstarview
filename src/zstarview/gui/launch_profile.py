from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from appdirs import user_config_dir

from ..cli.args import parse_args
from ..paths import APP_AUTHOR, APP_ID

GUI_LAUNCH_PROFILE_FILENAME = "gui-launch-profile.json"


def _profile_file() -> Path:
    return Path(user_config_dir(APP_ID, APP_AUTHOR)) / GUI_LAUNCH_PROFILE_FILENAME


def _jsonable(value: Any) -> Any:
    if isinstance(value, set):
        return sorted((_jsonable(item) for item in value), key=repr)
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    return value


def _default_profile() -> dict[str, Any]:
    return dict(vars(parse_args([])))


def load_gui_launch_profile() -> dict[str, Any]:
    path = _profile_file()
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}
    if not isinstance(raw, dict):
        return {}
    return {str(key): value for key, value in raw.items()}


def save_gui_launch_profile(profile: Mapping[str, Any]) -> None:
    path = _profile_file()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_jsonable(dict(profile)), ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def reset_gui_launch_profile() -> None:
    save_gui_launch_profile(_default_profile())


def default_gui_launch_profile() -> dict[str, Any]:
    return _default_profile()
