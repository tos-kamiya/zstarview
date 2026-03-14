from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

from .data.plateau_derived_tiles import parse_derived_tile_buildings, select_derived_tile_envelopes
from .data.urban_debug_layer_from_citygml import compute_debug_outlines
from .paths import PLATEAU_DERIVED_ROOT_DIR
from .types import ViewerData


def resolve_urban_debug_layer_for_viewer(
    viewer_data: ViewerData,
    *,
    derived_root_dir: str | Path = PLATEAU_DERIVED_ROOT_DIR,
) -> list[list[tuple[float, float]]] | None:
    return _build_dynamic_urban_outline_layer(
        lat_deg=float(viewer_data.lat_deg),
        lon_deg=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m),
        derived_root_dir=Path(derived_root_dir),
    )


@lru_cache(maxsize=64)
def _build_dynamic_urban_outline_layer(
    *,
    lat_deg: float,
    lon_deg: float,
    observer_height_m: float,
    derived_root_dir: Path,
) -> list[list[tuple[float, float]]] | None:
    if not derived_root_dir.exists():
        return None
    candidate_dirs = _list_derived_dirs(derived_root_dir)
    if not candidate_dirs:
        return None

    buildings = []
    for derived_dir in candidate_dirs:
        try:
            envelopes = select_derived_tile_envelopes(
                derived_dir,
                observer_lat_deg=lat_deg,
                observer_lon_deg=lon_deg,
                radius_km=3.0,
            )
        except ValueError:
            continue
        for envelope in envelopes:
            buildings.extend(parse_derived_tile_buildings(envelope.path))
    if not buildings:
        return None

    result = compute_debug_outlines(
        SimpleNamespace(
            id="coords",
            name="coords",
            latitude_deg=lat_deg,
            longitude_deg=lon_deg,
            observer_height_m=observer_height_m,
        ),
        tuple(buildings),
        radius_km=3.0,
        edge_sample_step_m=10.0,
    )
    outlines = [
        [(point.altitude_deg, point.azimuth_deg) for point in outline]
        for outline in result.outlines
    ]
    return outlines or None


@lru_cache(maxsize=8)
def _list_derived_dirs(derived_root_dir: Path) -> tuple[Path, ...]:
    return tuple(
        path
        for path in sorted(derived_root_dir.glob("*/bldg"))
        if path.is_dir()
    )
