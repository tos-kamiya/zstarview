"""Shared opacity and conservative observer-centric wall visibility helpers."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

NEARBY_WALL_MIN_DISTANCE_M = 5.0
DEFAULT_WALL_RADIUS_KM = 0.5
NEARBY_WALL_MAX_DISTANCE_M = DEFAULT_WALL_RADIUS_KM * 1000.0
NEARBY_WALL_FULL_OPACITY_DISTANCE_M = 10.0
NEARBY_WALL_OUTER_OPACITY = 0.1
WALL_OPACITY_LEVELS = 32
_AZIMUTH_BIN_DEG = 5.0
_AZIMUTH_BIN_COUNT = round(360.0 / _AZIMUTH_BIN_DEG)


def wall_opacity_level(
    distance_m: float, max_distance_m: float = NEARBY_WALL_MAX_DISTANCE_M
) -> int:
    if max_distance_m <= NEARBY_WALL_MIN_DISTANCE_M or distance_m > max_distance_m:
        return 0
    near = max(
        0.0,
        min(
            1.0,
            (distance_m - NEARBY_WALL_MIN_DISTANCE_M)
            / (NEARBY_WALL_FULL_OPACITY_DISTANCE_M - NEARBY_WALL_MIN_DISTANCE_M),
        ),
    )
    fade_start = NEARBY_WALL_MIN_DISTANCE_M
    if max_distance_m > NEARBY_WALL_FULL_OPACITY_DISTANCE_M:
        fade_start = NEARBY_WALL_FULL_OPACITY_DISTANCE_M
    far = max(
        0.0,
        min(
            1.0,
            (max_distance_m - distance_m) / (max_distance_m - fade_start),
        ),
    )
    if max_distance_m > NEARBY_WALL_FULL_OPACITY_DISTANCE_M:
        far = NEARBY_WALL_OUTER_OPACITY + (1.0 - NEARBY_WALL_OUTER_OPACITY) * far
    return round(WALL_OPACITY_LEVELS * min(near, far))


@dataclass(frozen=True)
class WallPrism:
    """One vertical extrusion of a plan ring, heights relative to the ground."""

    ring_lonlat: tuple[tuple[float, float], ...]
    bottom_m: float
    top_m: float


def roof_floor_prisms(
    roof_surfaces_lonlat: tuple[tuple[tuple[float, float, float], ...], ...],
    ground_elevation_m: float,
) -> tuple[WallPrism, ...]:
    """Extrude each LOD2 roof polygon down to the next lower floor of its building.

    A roof polygon's floor is its mean vertex elevation above the ground,
    rounded to the nearest metre.  Floors are ordered over the whole building
    (whether or not they overlap in plan) and each polygon is extruded from
    the next lower floor, or from the ground for the lowest floor, up to its
    own floor.  Polygons at or below the ground are skipped.
    """

    surfaces: list[tuple[int, tuple[tuple[float, float], ...]]] = []
    for surface in roof_surfaces_lonlat:
        points = surface[:-1] if len(surface) > 1 and surface[0] == surface[-1] else surface
        if len(points) < 3:
            continue
        mean_elevation = sum(point[2] for point in points) / len(points)
        floor = round(mean_elevation - ground_elevation_m)
        if floor <= 0:
            continue
        ring = tuple((point[0], point[1]) for point in points)
        surfaces.append((floor, ring + (ring[0],)))
    floors = sorted({floor for floor, _ring in surfaces})
    lower_floor = {
        floor: (floors[index - 1] if index else 0) for index, floor in enumerate(floors)
    }
    return tuple(
        WallPrism(
            ring_lonlat=ring,
            bottom_m=float(lower_floor[floor]),
            top_m=float(floor),
        )
        for floor, ring in surfaces
    )


@dataclass(frozen=True)
class WallSegment:
    start: tuple[float, float]
    end: tuple[float, float]
    bottom: float
    top: float
    distance_m: float
    height_m: float
    step_m: float
    building_id: str | None = None


_OCCLUSION_EPSILON = 1e-7
_OCCLUSION_CHUNK_SIZE = 256


def _azimuth(point: tuple[float, float]) -> float:
    return math.degrees(math.atan2(point[0], point[1])) % 360.0


def visible_wall_segments(
    segments: list[WallSegment], max_distance_m: float = NEARBY_WALL_MAX_DISTANCE_M
) -> list[WallSegment]:
    """Remove walls provably covered by a single higher-opacity wall.

    Azimuth bins only narrow the search.  Exact 3D ray/rectangle tests make
    the final decision, and ambiguous or partially covered walls are retained.
    A target is covered when both bottom-edge endpoints, cast as rays from the
    observer, cross the interior of one cover wall strictly below its top and
    above its bottom.  The tests are evaluated as NumPy arrays per azimuth bin.
    """
    count = len(segments)
    if not count:
        return []
    start_x = np.fromiter((s.start[0] for s in segments), dtype=np.float64, count=count)
    start_y = np.fromiter((s.start[1] for s in segments), dtype=np.float64, count=count)
    end_x = np.fromiter((s.end[0] for s in segments), dtype=np.float64, count=count)
    end_y = np.fromiter((s.end[1] for s in segments), dtype=np.float64, count=count)
    bottoms = np.fromiter((s.bottom for s in segments), dtype=np.float64, count=count)
    tops = np.fromiter((s.top for s in segments), dtype=np.float64, count=count)
    levels = np.fromiter(
        (wall_opacity_level(s.distance_m, max_distance_m) for s in segments),
        dtype=np.int64,
        count=count,
    )

    bins: list[list[int]] = [[] for _ in range(_AZIMUTH_BIN_COUNT)]
    target_bins = np.empty(count, dtype=np.intp)
    for i, segment in enumerate(segments):
        a, b = _azimuth(segment.start), _azimuth(segment.end)
        delta = (b - a + 180.0) % 360.0 - 180.0
        start = a if delta >= 0.0 else b
        span = abs(delta)
        first = math.floor(start / _AZIMUTH_BIN_DEG)
        last = math.floor((start + span) / _AZIMUTH_BIN_DEG)
        for index in range(first, last + 1):
            bins[index % _AZIMUTH_BIN_COUNT].append(i)
        midpoint = (
            (segment.start[0] + segment.end[0]) * 0.5,
            (segment.start[1] + segment.end[1]) * 0.5,
        )
        target_bins[i] = int(_azimuth(midpoint) // _AZIMUTH_BIN_DEG)

    eps = _OCCLUSION_EPSILON
    covered = np.zeros(count, dtype=bool)
    for bin_index, members in enumerate(bins):
        targets = np.flatnonzero(target_bins == bin_index)
        if not members or not targets.size:
            continue
        covers = np.asarray(members, dtype=np.intp)
        cover_x, cover_y = start_x[covers], start_y[covers]
        edge_x, edge_y = end_x[covers] - cover_x, end_y[covers] - cover_y
        # cross(cover.start, edge)
        numerator = cover_x * edge_y - cover_y * edge_x
        cover_levels = levels[covers]
        cover_bottoms = bottoms[covers] + eps
        cover_tops = tops[covers] - eps
        for offset in range(0, targets.size, _OCCLUSION_CHUNK_SIZE):
            chunk = targets[offset:offset + _OCCLUSION_CHUNK_SIZE]
            candidate = (cover_levels[None, :] >= levels[chunk][:, None]) & (
                covers[None, :] != chunk[:, None]
            )
            target_bottoms = bottoms[chunk][:, None]
            target_tops = tops[chunk][:, None]
            for point_x, point_y in (
                (start_x[chunk][:, None], start_y[chunk][:, None]),
                (end_x[chunk][:, None], end_y[chunk][:, None]),
            ):
                # cross(point, edge); near-parallel rays are never covered.
                denominator = point_x * edge_y - point_y * edge_x
                with np.errstate(divide="ignore", invalid="ignore"):
                    t = numerator / denominator
                    u = (cover_x * point_y - cover_y * point_x) / denominator
                # Keep boundary cases: a tiny numerical discrepancy must not
                # remove a visible edge or a wall touching/intersecting the
                # occluder.
                candidate &= (
                    (np.abs(denominator) >= 1e-9)
                    & (t > eps) & (t < 1.0 - eps)
                    & (u > eps) & (u < 1.0 - eps)
                    & (cover_bottoms < target_bottoms * t)
                    & (target_tops * t < cover_tops)
                )
            covered[chunk] = candidate.any(axis=1)
    return [segment for segment, hidden in zip(segments, covered.tolist(), strict=True)
            if not hidden]
