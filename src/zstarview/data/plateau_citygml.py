from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

from zstarview.data.urban_outline_common import (
    BuildingFootprint,
    DEFAULT_MIN_BUILDING_HEIGHT_M,
)

_ENVELOPE_RE = re.compile(r"<gml:(lowerCorner|upperCorner)>([^<]+)</gml:")
_GML_NS = {
    "bldg": "http://www.opengis.net/citygml/building/2.0",
    "gml": "http://www.opengis.net/gml",
}
DEFAULT_STOREY_HEIGHT_M = 3.5
_INVALID_CITYGML_NUMERIC_VALUES = {-9999.0}


@dataclass(frozen=True)
class TileEnvelope:
    path: Path
    min_lat_deg: float
    min_lon_deg: float
    max_lat_deg: float
    max_lon_deg: float


def load_tile_envelope(path: Path) -> TileEnvelope | None:
    with path.open("r", encoding="utf-8-sig", errors="ignore") as handle:
        head_lines: list[str] = []
        for _ in range(8):
            line = handle.readline()
            if not line:
                break
            head_lines.append(line)
    text = "".join(head_lines)
    matches = _ENVELOPE_RE.findall(text)
    if len(matches) < 2:
        return None
    lower = [float(value) for value in matches[0][1].split()]
    upper = [float(value) for value in matches[1][1].split()]
    return TileEnvelope(
        path=path,
        min_lat_deg=min(lower[0], upper[0]),
        min_lon_deg=min(lower[1], upper[1]),
        max_lat_deg=max(lower[0], upper[0]),
        max_lon_deg=max(lower[1], upper[1]),
    )


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


def parse_citygml_buildings(
    path: Path,
    *,
    min_building_height_m: float = DEFAULT_MIN_BUILDING_HEIGHT_M,
) -> tuple[BuildingFootprint, ...]:
    root = ET.parse(path).getroot()
    buildings: list[BuildingFootprint] = []
    for index, building in enumerate(root.findall(".//bldg:Building", _GML_NS)):
        height_m = resolve_building_height_m(building)
        if height_m is None:
            continue
        if not math.isfinite(height_m) or height_m < float(min_building_height_m):
            continue

        rings: list[tuple[tuple[float, float], ...]] = []
        for pos in building.findall(".//bldg:lod0RoofEdge//gml:posList", _GML_NS):
            ring = parse_pos_list_ring(pos.text)
            if ring:
                rings.append(ring)
        if not rings:
            continue
        building_id = building.attrib.get("{http://www.opengis.net/gml}id", f"bldg-{index}")
        buildings.append(
            BuildingFootprint(
                building_id=building_id,
                height_m=height_m,
                rings_lonlat=tuple(rings),
            )
        )
    return tuple(buildings)


def resolve_building_height_m(building: ET.Element) -> float | None:
    measured_height_m = parse_citygml_numeric(
        building.findtext("bldg:measuredHeight", default="", namespaces=_GML_NS)
    )
    if measured_height_m is not None:
        return measured_height_m

    storeys = parse_citygml_numeric(
        building.findtext("bldg:storeysAboveGround", default="", namespaces=_GML_NS)
    )
    if storeys is None:
        return None
    if storeys <= 0 or storeys >= 9999 or not float(storeys).is_integer():
        return None
    return storeys * DEFAULT_STOREY_HEIGHT_M


def parse_citygml_numeric(text: str | None) -> float | None:
    if text is None:
        return None
    try:
        value = float(text.strip())
    except (AttributeError, ValueError):
        return None
    if not math.isfinite(value) or value in _INVALID_CITYGML_NUMERIC_VALUES:
        return None
    return value


def parse_pos_list_ring(text: str | None) -> tuple[tuple[float, float], ...]:
    if not text:
        return ()
    values = [float(value) for value in text.split()]
    if len(values) < 12:
        return ()
    coords: list[tuple[float, float]] = []
    for index in range(0, len(values), 3):
        lat_deg = values[index]
        lon_deg = values[index + 1]
        point = (lon_deg, lat_deg)
        if not coords or coords[-1] != point:
            coords.append(point)
    if len(coords) < 4:
        return ()
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    return tuple(coords)
