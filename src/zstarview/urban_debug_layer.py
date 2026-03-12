from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

from .data.plateau_derived_tiles import parse_derived_tile_buildings, select_derived_tile_envelopes
from .data.urban_debug_layer_from_citygml import compute_debug_outlines
from .paths import TOKYO23_PLATEAU_DERIVED_DIR
from .types import ViewerData


def resolve_urban_debug_layer_for_viewer(
    viewer_data: ViewerData,
    *,
    tokyo23_derived_dir: str | Path = TOKYO23_PLATEAU_DERIVED_DIR,
) -> list[list[tuple[float, float]]] | None:
    return _build_dynamic_tokyo23_debug_layer(
        lat_deg=float(viewer_data.lat_deg),
        lon_deg=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m),
        derived_dir=Path(tokyo23_derived_dir),
    )


@lru_cache(maxsize=64)
def _build_dynamic_tokyo23_debug_layer(
    *,
    lat_deg: float,
    lon_deg: float,
    observer_height_m: float,
    derived_dir: Path,
) -> list[list[tuple[float, float]]] | None:
    if not derived_dir.exists():
        return None
    try:
        envelopes = select_derived_tile_envelopes(
            derived_dir,
            observer_lat_deg=lat_deg,
            observer_lon_deg=lon_deg,
            radius_km=3.0,
        )
    except ValueError:
        return None

    buildings = []
    for envelope in envelopes:
        buildings.extend(parse_derived_tile_buildings(envelope.path, min_building_height_m=40.0))

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
