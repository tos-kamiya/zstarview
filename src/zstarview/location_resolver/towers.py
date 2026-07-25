from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from ..paths import TOWER_VIEWPOINTS_FILE
from .viewpoints import (
    Viewpoint,
    build_viewpoint,
    list_viewpoint_all_names,
    list_viewpoint_primary_names,
    load_viewpoints,
    resolve_viewpoint,
    viewpoint_to_dict,
)

TowerViewpoint = Viewpoint

_TOWER_META_KEYS = (
    "classes",
    "class_qids",
    "wikidata_url",
    "sources",
    "viewpoint_name",
    "location_arg",
    "slug",
    "viewpoint_types",
)


def _as_tower(item: dict[str, object]) -> TowerViewpoint:
    return build_viewpoint(
        item,
        kind="tower",
        height_key="height_m",
        viewpoint_height_key="viewpoint_height_m",
        meta_keys=_TOWER_META_KEYS,
    )


@lru_cache(maxsize=1)
def load_tower_viewpoints(
    path: str | Path = TOWER_VIEWPOINTS_FILE,
) -> tuple[TowerViewpoint, ...]:
    return load_viewpoints(path, builder=_as_tower, error_label="tower viewpoints")


def list_tower_primary_names(
    towers: tuple[TowerViewpoint, ...] | None = None,
) -> tuple[str, ...]:
    towers = load_tower_viewpoints() if towers is None else towers
    return list_viewpoint_primary_names(towers)


def list_tower_all_names(
    towers: tuple[TowerViewpoint, ...] | None = None,
) -> tuple[str, ...]:
    towers = load_tower_viewpoints() if towers is None else towers
    return list_viewpoint_all_names(towers)


def tower_viewpoint_to_dict(tower: TowerViewpoint) -> dict[str, object]:
    return viewpoint_to_dict(tower)


def resolve_tower_viewpoint(
    query: str,
    towers: tuple[TowerViewpoint, ...] | None = None,
) -> TowerViewpoint | None:
    towers = load_tower_viewpoints() if towers is None else towers
    return resolve_viewpoint(query, towers, rank_key=lambda tower: tower.height_m)
