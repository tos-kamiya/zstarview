from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

from .data.derived_tile_cache import parse_derived_tile_buildings, select_derived_tile_envelopes
from .data.import_overture_buildings import DEFAULT_FETCH_RADIUS_KM
from .data.urban_outline_common import BuildingFootprint
from .data.urban_outline_from_buildings import compute_urban_outlines
from .paths import OVERTURE_DERIVED_ROOT_DIR
from .types import UrbanOutlinePolyline, ViewerData


def resolve_urban_outline_layer_for_viewer(
    viewer_data: ViewerData,
    *,
    derived_root_dir: str | Path = OVERTURE_DERIVED_ROOT_DIR,
    derived_dir: str | Path | None = None,
    derived_dirs: tuple[str | Path, ...] | None = None,
    radius_km: float = DEFAULT_FETCH_RADIUS_KM,
    min_distance_km: float = 0.0,
) -> list[UrbanOutlinePolyline] | None:
    return _build_dynamic_urban_outline_layer(
        lat_deg=float(viewer_data.lat_deg),
        lon_deg=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m),
        derived_root_dir=Path(derived_root_dir),
        derived_dir=None if derived_dir is None else Path(derived_dir),
        derived_dirs=None if derived_dirs is None else tuple(Path(path) for path in derived_dirs),
        radius_km=float(radius_km),
        min_distance_km=float(min_distance_km),
    )


@lru_cache(maxsize=64)
def _build_dynamic_urban_outline_layer(
    *,
    lat_deg: float,
    lon_deg: float,
    observer_height_m: float,
    derived_root_dir: Path,
    derived_dir: Path | None = None,
    derived_dirs: tuple[Path, ...] | None = None,
    radius_km: float = DEFAULT_FETCH_RADIUS_KM,
    min_distance_km: float = 0.0,
) -> list[UrbanOutlinePolyline] | None:
    if derived_dirs is not None:
        candidate_dirs = tuple(path for path in derived_dirs if path.exists())
    elif derived_dir is not None:
        candidate_dirs = (derived_dir,) if derived_dir.exists() else ()
    elif derived_root_dir.exists():
        candidate_dirs = _list_derived_dirs(derived_root_dir)
    else:
        candidate_dirs = ()
    if not candidate_dirs:
        return None

    buildings = []
    for derived_dir in candidate_dirs:
        try:
            envelopes = select_derived_tile_envelopes(
                derived_dir,
                observer_lat_deg=lat_deg,
                observer_lon_deg=lon_deg,
                radius_km=radius_km,
            )
        except ValueError:
            continue
        for envelope in envelopes:
            buildings.extend(parse_derived_tile_buildings(envelope.path))
    buildings = _merge_building_footprints(tuple(buildings))
    if not buildings:
        return None

    result = compute_urban_outlines(
        SimpleNamespace(
            id="coords",
            name="coords",
            latitude_deg=lat_deg,
            longitude_deg=lon_deg,
            viewpoint_height_m=observer_height_m,
            observer_height_m=observer_height_m,
        ),
        tuple(buildings),
        radius_km=radius_km,
        min_distance_km=min_distance_km,
        edge_sample_step_m=10.0,
    )
    outlines = [
        UrbanOutlinePolyline(
            points=[(point.altitude_deg, point.azimuth_deg) for point in outline.points],
            height_m=float(outline.height_m),
            source="base",
        )
        for outline in result.outlines
    ]
    return outlines or None


def _merge_building_footprints(buildings: tuple[BuildingFootprint, ...]) -> tuple[BuildingFootprint, ...]:
    parent_ids_with_parts = {
        building.parent_building_id
        for building in buildings
        if building.parent_building_id
    }
    if not parent_ids_with_parts:
        return buildings
    return tuple(
        building
        for building in buildings
        if building.parent_building_id is not None or building.building_id not in parent_ids_with_parts
    )


@lru_cache(maxsize=8)
def _list_derived_dirs(derived_root_dir: Path) -> tuple[Path, ...]:
    return tuple(
        path
        for path in sorted(derived_root_dir.glob("*/bldg"))
        if path.is_dir()
    )
