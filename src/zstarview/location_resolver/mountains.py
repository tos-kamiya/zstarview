from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from ..paths import MOUNTAIN_VIEWPOINTS_FILE
from .viewpoints import (
    Viewpoint,
    build_viewpoint,
    list_viewpoint_all_names,
    list_viewpoint_primary_names,
    load_viewpoints,
    resolve_viewpoint,
    viewpoint_to_dict,
)

MountainViewpoint = Viewpoint

_MOUNTAIN_META_KEYS = (
    "wikidata_url",
    "wikipedia_urls",
    "location_arg",
    "slug",
    "elevation_m",
    "region_hint",
)


def _as_mountain(item: dict[str, object]) -> MountainViewpoint:
    return build_viewpoint(
        item,
        kind="mountain",
        height_key="height_m",
        meta_keys=_MOUNTAIN_META_KEYS,
    )


@lru_cache(maxsize=1)
def load_mountain_viewpoints(
    path: str | Path = MOUNTAIN_VIEWPOINTS_FILE,
) -> tuple[MountainViewpoint, ...]:
    return load_viewpoints(
        path,
        builder=_as_mountain,
        error_label="mountain viewpoints",
    )


def list_mountain_primary_names(
    mountains: tuple[MountainViewpoint, ...] | None = None,
) -> tuple[str, ...]:
    mountains = load_mountain_viewpoints() if mountains is None else mountains
    return list_viewpoint_primary_names(mountains)


def list_mountain_all_names(
    mountains: tuple[MountainViewpoint, ...] | None = None,
) -> tuple[str, ...]:
    mountains = load_mountain_viewpoints() if mountains is None else mountains
    return list_viewpoint_all_names(mountains)


def mountain_viewpoint_to_dict(mountain: MountainViewpoint) -> dict[str, object]:
    return viewpoint_to_dict(mountain)


def resolve_mountain_viewpoint(
    query: str,
    mountains: tuple[MountainViewpoint, ...] | None = None,
) -> MountainViewpoint | None:
    mountains = load_mountain_viewpoints() if mountains is None else mountains
    return resolve_viewpoint(
        query,
        mountains,
        rank_key=lambda mountain: float(mountain.meta.get("elevation_m", 0.0)),
    )
