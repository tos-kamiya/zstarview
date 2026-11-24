"""Shared opacity and conservative observer-centric wall visibility helpers."""

from __future__ import annotations

import math
from dataclasses import dataclass

NEARBY_WALL_MIN_DISTANCE_M = 5.0
NEARBY_WALL_MAX_DISTANCE_M = 200.0
NEARBY_WALL_FULL_OPACITY_DISTANCE_M = 10.0
WALL_OPACITY_LEVELS = 32
_AZIMUTH_BIN_DEG = 5.0
_AZIMUTH_BIN_COUNT = round(360.0 / _AZIMUTH_BIN_DEG)


def wall_opacity_level(
    distance_m: float, max_distance_m: float = NEARBY_WALL_MAX_DISTANCE_M
) -> int:
    if max_distance_m <= NEARBY_WALL_MIN_DISTANCE_M:
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


def _cross(a: tuple[float, float], b: tuple[float, float]) -> float:
    return a[0] * b[1] - a[1] * b[0]


def _covered_by_wall(
    target: WallSegment, cover: WallSegment, max_distance_m: float
) -> bool:
    """Require all target corners to lie strictly behind one opaque-enough wall."""
    if wall_opacity_level(cover.distance_m, max_distance_m) < wall_opacity_level(
        target.distance_m, max_distance_m
    ):
        return False
    edge = (cover.end[0] - cover.start[0], cover.end[1] - cover.start[1])
    numerator = _cross(cover.start, edge)
    for point in (target.start, target.end):
        denominator = _cross(point, edge)
        if abs(denominator) < 1e-9:
            return False
        t = numerator / denominator
        u = _cross(cover.start, point) / denominator
        # Keep boundary cases: a tiny numerical discrepancy must not remove a
        # visible edge or a wall touching/intersecting the occluder.
        if not (1e-7 < t < 1.0 - 1e-7 and 1e-7 < u < 1.0 - 1e-7):
            return False
        if not (
            cover.bottom + 1e-7 < target.bottom * t
            and target.top * t < cover.top - 1e-7
        ):
            return False
    return True


def _azimuth(point: tuple[float, float]) -> float:
    return math.degrees(math.atan2(point[0], point[1])) % 360.0


def visible_wall_segments(
    segments: list[WallSegment], max_distance_m: float = NEARBY_WALL_MAX_DISTANCE_M
) -> list[WallSegment]:
    """Remove walls provably covered by a single higher-opacity wall.

    Azimuth bins only narrow the search.  Exact 3D ray/rectangle tests make
    the final decision, and ambiguous or partially covered walls are retained.
    """
    bins: list[list[int]] = [[] for _ in range(_AZIMUTH_BIN_COUNT)]
    for i, segment in enumerate(segments):
        a, b = _azimuth(segment.start), _azimuth(segment.end)
        delta = (b - a + 180.0) % 360.0 - 180.0
        start = a if delta >= 0.0 else b
        span = abs(delta)
        first = math.floor(start / _AZIMUTH_BIN_DEG)
        last = math.floor((start + span) / _AZIMUTH_BIN_DEG)
        for index in range(first, last + 1):
            bins[index % _AZIMUTH_BIN_COUNT].append(i)

    visible: list[WallSegment] = []
    for i, target in enumerate(segments):
        midpoint = (
            (target.start[0] + target.end[0]) * 0.5,
            (target.start[1] + target.end[1]) * 0.5,
        )
        bin_index = int(_azimuth(midpoint) // _AZIMUTH_BIN_DEG)
        if not any(
            j != i and _covered_by_wall(target, segments[j], max_distance_m)
            for j in bins[bin_index]
        ):
            visible.append(target)
    return visible
