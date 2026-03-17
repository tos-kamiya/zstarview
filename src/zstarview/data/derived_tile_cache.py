from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from zstarview.data.urban_outline_common import BuildingFootprint


@dataclass(frozen=True)
class TileEnvelope:
    path: Path
    min_lat_deg: float
    min_lon_deg: float
    max_lat_deg: float
    max_lon_deg: float


def approx_distance_km(
    lat0_deg: float,
    lon0_deg: float,
    lat1_deg: float,
    lon1_deg: float,
) -> float:
    dlon_km = (lon1_deg - lon0_deg) * 111.32 * math.cos(math.radians((lat0_deg + lat1_deg) * 0.5))
    dlat_km = (lat1_deg - lat0_deg) * 111.32
    return math.hypot(dlon_km, dlat_km)


def envelope_min_distance_km(
    envelope: TileEnvelope,
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
) -> float:
    nearest_lat = min(max(observer_lat_deg, envelope.min_lat_deg), envelope.max_lat_deg)
    nearest_lon = min(max(observer_lon_deg, envelope.min_lon_deg), envelope.max_lon_deg)
    return approx_distance_km(observer_lat_deg, observer_lon_deg, nearest_lat, nearest_lon)


def envelope_max_distance_km(
    envelope: TileEnvelope,
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
) -> float:
    return max(
        approx_distance_km(observer_lat_deg, observer_lon_deg, lat_deg, lon_deg)
        for lat_deg, lon_deg in (
            (envelope.min_lat_deg, envelope.min_lon_deg),
            (envelope.min_lat_deg, envelope.max_lon_deg),
            (envelope.max_lat_deg, envelope.min_lon_deg),
            (envelope.max_lat_deg, envelope.max_lon_deg),
        )
    )


def load_derived_tile_envelope(path: Path) -> TileEnvelope | None:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return None
    tile = payload.get("tile")
    if not isinstance(tile, dict):
        return None
    bbox = tile.get("bbox")
    if not isinstance(bbox, dict):
        return None
    try:
        min_lat_deg = float(bbox["min_lat"])
        min_lon_deg = float(bbox["min_lon"])
        max_lat_deg = float(bbox["max_lat"])
        max_lon_deg = float(bbox["max_lon"])
    except (KeyError, TypeError, ValueError):
        return None
    return TileEnvelope(
        path=path,
        min_lat_deg=min_lat_deg,
        min_lon_deg=min_lon_deg,
        max_lat_deg=max_lat_deg,
        max_lon_deg=max_lon_deg,
    )


def load_derived_tile_index(derived_dir: Path) -> tuple[TileEnvelope | None, tuple[TileEnvelope, ...]] | None:
    index_path = derived_dir / "tile_index.json"
    if not index_path.exists():
        return None
    payload = json.loads(index_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return None

    city_envelope = _parse_index_bbox(payload.get("bbox"), path=index_path)
    raw_tiles = payload.get("tiles")
    if not isinstance(raw_tiles, list):
        return None
    tiles: list[TileEnvelope] = []
    for row in raw_tiles:
        if not isinstance(row, dict):
            continue
        tile_path = row.get("path")
        if not isinstance(tile_path, str) or not tile_path:
            continue
        envelope = _parse_index_bbox(row.get("bbox"), path=derived_dir / tile_path)
        if envelope is None:
            continue
        tiles.append(envelope)
    return city_envelope, tuple(tiles)


def select_derived_tile_envelopes(
    derived_dir: Path,
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    radius_km: float,
) -> tuple[TileEnvelope, ...]:
    indexed = load_derived_tile_index(derived_dir)
    if indexed is not None:
        city_envelope, indexed_tiles = indexed
        if city_envelope is not None and envelope_min_distance_km(
            city_envelope,
            observer_lat_deg=observer_lat_deg,
            observer_lon_deg=observer_lon_deg,
        ) > radius_km:
            raise ValueError(f"No derived building tiles found within {radius_km:.1f} km.")
        envelopes = tuple(
            envelope
            for envelope in indexed_tiles
            if envelope_min_distance_km(
                envelope,
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
            ) <= radius_km
        )
        if not envelopes:
            raise ValueError(f"No derived building tiles found within {radius_km:.1f} km.")
        return envelopes

    envelopes: list[TileEnvelope] = []
    for path in sorted(derived_dir.glob("*.json")):
        if path.name == "tile_index.json":
            continue
        envelope = load_derived_tile_envelope(path)
        if envelope is None:
            continue
        if envelope_min_distance_km(
            envelope,
            observer_lat_deg=observer_lat_deg,
            observer_lon_deg=observer_lon_deg,
        ) <= radius_km:
            envelopes.append(envelope)
    if not envelopes:
        raise ValueError(f"No derived building tiles found within {radius_km:.1f} km.")
    return tuple(envelopes)


def parse_derived_tile_buildings(
    path: Path,
    *,
    min_building_height_m: float | None = None,
) -> tuple[BuildingFootprint, ...]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return ()
    raw_buildings = payload.get("buildings")
    if not isinstance(raw_buildings, list):
        return ()

    buildings: list[BuildingFootprint] = []
    for index, row in enumerate(raw_buildings):
        if not isinstance(row, dict):
            continue
        try:
            height_m = float(row["height_m"])
        except (KeyError, TypeError, ValueError):
            continue
        if min_building_height_m is not None and height_m < float(min_building_height_m):
            continue
        rings = _parse_rings_lonlat(row.get("rings"))
        if not rings:
            continue
        building_id = row.get("id")
        if not isinstance(building_id, str) or not building_id:
            building_id = f"bldg-{index}"
        parent_building_id = row.get("parent_building_id")
        if not isinstance(parent_building_id, str) or not parent_building_id:
            parent_building_id = None
        buildings.append(
            BuildingFootprint(
                building_id=building_id,
                height_m=height_m,
                rings_lonlat=rings,
                parent_building_id=parent_building_id,
            )
        )
    return tuple(buildings)


def _parse_rings_lonlat(raw_rings: object) -> tuple[tuple[tuple[float, float], ...], ...]:
    if not isinstance(raw_rings, list):
        return ()
    rings: list[tuple[tuple[float, float], ...]] = []
    for raw_ring in raw_rings:
        if not isinstance(raw_ring, list):
            continue
        points: list[tuple[float, float]] = []
        for raw_point in raw_ring:
            if not isinstance(raw_point, Sequence) or len(raw_point) != 2:
                continue
            try:
                lon = float(raw_point[0])
                lat = float(raw_point[1])
            except (TypeError, ValueError):
                continue
            points.append((lon, lat))
        if len(points) >= 4:
            rings.append(tuple(points))
    return tuple(rings)


def _parse_index_bbox(raw_bbox: object, *, path: Path) -> TileEnvelope | None:
    if not isinstance(raw_bbox, dict):
        return None
    try:
        min_lat_deg = float(raw_bbox["min_lat"])
        min_lon_deg = float(raw_bbox["min_lon"])
        max_lat_deg = float(raw_bbox["max_lat"])
        max_lon_deg = float(raw_bbox["max_lon"])
    except (KeyError, TypeError, ValueError):
        return None
    return TileEnvelope(
        path=path,
        min_lat_deg=min_lat_deg,
        min_lon_deg=min_lon_deg,
        max_lat_deg=max_lat_deg,
        max_lon_deg=max_lon_deg,
    )
