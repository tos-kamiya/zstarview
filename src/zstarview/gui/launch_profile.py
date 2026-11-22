from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from appdirs import user_config_dir

from ..cli.args import _parse_cloud_stripe, parse_args
from ..paths import APP_AUTHOR, APP_ID, THEME_PRESET_NAMES

GUI_LAUNCH_PROFILE_FILENAME = "gui-launch-profile.json"
# Version 2 (zstarview 2.0): -A default 5, location source, removed themes.
GUI_LAUNCH_PROFILE_VERSION = 2
_PROFILE_VERSION_KEY = "profile_version"
LOCATION_SOURCES = ("auto", "last", "city", "place")
_PRE_V2_VIEW_CENTER_ALT_DEFAULT = 90.0
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
    profile = _without_ignored_keys(vars(parse_args([])))
    profile["location_source"] = "auto"
    return profile


def location_source_from_profile(profile: Mapping[str, Any]) -> str:
    """Return the startup location source, inferring it for older profiles."""
    source = str(profile.get("location_source", "") or "").strip().lower()
    if source in LOCATION_SOURCES:
        return source
    city_value = profile.get("city")
    if isinstance(city_value, dict):
        resolver = str(city_value.get("resolver", "")).strip().lower()
        return "place" if resolver == "nominatim" else "city"
    if isinstance(city_value, str) and city_value.strip():
        return "city"
    return "auto"


def _migrate_profile(raw: dict[str, Any]) -> dict[str, Any]:
    version = raw.get(_PROFILE_VERSION_KEY, 1)
    if not isinstance(version, int):
        version = 1
    profile = dict(raw)
    if version < 2:
        # 2.0 changed the default view altitude from 90 to 5. A stored 90 is
        # indistinguishable from the old default, so use the new default.
        alt = profile.get("view_center_alt")
        if isinstance(alt, (int, float)) and float(alt) == _PRE_V2_VIEW_CENTER_ALT_DEFAULT:
            profile.pop("view_center_alt")
    if "location_source" not in profile and "city" in profile:
        # Profiles saved before 2.0 had no explicit source; keep their behavior
        # instead of switching them to the new "auto" default.
        profile["location_source"] = location_source_from_profile(profile)
    theme = profile.get("theme")
    if isinstance(theme, str) and theme not in THEME_PRESET_NAMES:
        profile.pop("theme")
    profile.pop(_PROFILE_VERSION_KEY, None)
    return profile


def load_gui_launch_profile() -> dict[str, Any]:
    path = _profile_file()
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}
    if not isinstance(raw, dict):
        return {}
    return _without_ignored_keys(_migrate_profile(raw))


def save_gui_launch_profile(profile: Mapping[str, Any]) -> None:
    path = _profile_file()
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = _without_ignored_keys(profile)
    payload[_PROFILE_VERSION_KEY] = GUI_LAUNCH_PROFILE_VERSION
    path.write_text(
        json.dumps(
            _jsonable(payload),
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
