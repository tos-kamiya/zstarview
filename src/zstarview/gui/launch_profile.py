from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from appdirs import user_config_dir

from ..cli.args import _parse_cloud_stripe, parse_args
from ..paths import APP_AUTHOR, APP_ID

GUI_LAUNCH_PROFILE_FILENAME = "gui-launch-profile.json"
_IGNORED_PROFILE_KEYS = frozenset(
    {
        "enlarge_moon",
        "cloud_voxel",
        "cloud_mode",
        "cloud_voxel_threshold",
        "cloud_missing_tint_opacity",
        "set_cache_retention_hours",
        "initial_display_mode",
        "inverted_city",
        "sky_disc_style",
        "urban_outline_feature_type",
        "urban_outline_min_height_m",
    }
)


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


def _is_supported_cloud_stripe(value: Any) -> bool:
    """Return whether a stored cloud stripe is a current cutout style.

    Profiles saved before shell clouds were removed may hold shell styles
    such as ``halftone2,30,1.7``; those values are dropped on load.
    """
    if isinstance(value, (tuple, list)):
        value = ",".join(str(item) for item in value[:2])
    if not isinstance(value, str):
        return False
    try:
        _parse_cloud_stripe(value)
    except Exception:
        return False
    return True


def _without_ignored_keys(profile: Mapping[str, Any]) -> dict[str, Any]:
    return {
        str(key): value
        for key, value in profile.items()
        if str(key) not in _IGNORED_PROFILE_KEYS
        and (str(key) != "cloud_stripe" or _is_supported_cloud_stripe(value))
    }


def _default_profile() -> dict[str, Any]:
    return _without_ignored_keys(vars(parse_args([])))


def load_gui_launch_profile() -> dict[str, Any]:
    path = _profile_file()
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}
    if not isinstance(raw, dict):
        return {}
    return _without_ignored_keys(raw)


def save_gui_launch_profile(profile: Mapping[str, Any]) -> None:
    path = _profile_file()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            _jsonable(_without_ignored_keys(profile)),
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )


def reset_gui_launch_profile() -> None:
    save_gui_launch_profile(_default_profile())


def default_gui_launch_profile() -> dict[str, Any]:
    return _default_profile()
