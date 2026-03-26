from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
from pyproj import CRS, Transformer

from zstarview.location_resolver import TowerViewpoint


DEFAULT_MIN_BUILDING_HEIGHT_M = 40.0


@dataclass(frozen=True)
class BuildingFootprint:
    building_id: str
    height_m: float
    rings_lonlat: tuple[tuple[tuple[float, float], ...], ...]
    parent_building_id: str | None = None


def sanitize_slug(text: str) -> str:
    collapsed = re.sub(r"[^0-9A-Za-z]+", "_", text.strip())
    collapsed = collapsed.strip("_")
    return collapsed.lower() or "tower"


def make_local_transformer(tower: TowerViewpoint) -> Transformer:
    local_crs = CRS.from_proj4(
        f"+proj=aeqd +lat_0={tower.latitude_deg} +lon_0={tower.longitude_deg} "
        "+datum=WGS84 +units=m +no_defs"
    )
    return Transformer.from_crs("EPSG:4326", local_crs, always_xy=True)


def project_ring_xy(
    ring_lonlat: Sequence[tuple[float, float]],
    transformer: Transformer,
) -> np.ndarray:
    lon = np.array([point[0] for point in ring_lonlat], dtype=np.float64)
    lat = np.array([point[1] for point in ring_lonlat], dtype=np.float64)
    x, y = transformer.transform(lon, lat)
    return np.column_stack((np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)))


def bbox_min_distance_m(points_xy: np.ndarray) -> float:
    min_x = float(np.min(points_xy[:, 0]))
    max_x = float(np.max(points_xy[:, 0]))
    min_y = float(np.min(points_xy[:, 1]))
    max_y = float(np.max(points_xy[:, 1]))
    nearest_x = 0.0 if min_x <= 0.0 <= max_x else min(abs(min_x), abs(max_x))
    nearest_y = 0.0 if min_y <= 0.0 <= max_y else min(abs(min_y), abs(max_y))
    return math.hypot(nearest_x, nearest_y)


def sample_segment_points_xy(
    start_xy: np.ndarray,
    end_xy: np.ndarray,
    *,
    sample_step_m: float,
) -> np.ndarray:
    delta = end_xy - start_xy
    length = float(np.hypot(delta[0], delta[1]))
    count = max(2, int(math.ceil(length / max(sample_step_m, 0.1))) + 1)
    t = np.linspace(0.0, 1.0, num=count, dtype=np.float64)
    return start_xy[None, :] + t[:, None] * delta[None, :]


def sample_ring_points_xy(ring_xy: np.ndarray, *, sample_step_m: float) -> np.ndarray:
    segments: list[np.ndarray] = []
    for start_xy, end_xy in zip(ring_xy[:-1], ring_xy[1:]):
        segment = sample_segment_points_xy(start_xy, end_xy, sample_step_m=sample_step_m)
        if segments:
            segment = segment[1:]
        segments.append(segment)
    if not segments:
        return np.empty((0, 2), dtype=np.float64)
    return np.vstack(segments)


def iter_true_runs(mask: np.ndarray) -> Iterable[slice]:
    if mask.ndim != 1:
        raise ValueError("mask must be 1-dimensional.")
    start: int | None = None
    for index, flag in enumerate(mask.tolist()):
        if flag:
            if start is None:
                start = index
            continue
        if start is not None:
            yield slice(start, index)
            start = None
    if start is not None:
        yield slice(start, mask.size)
