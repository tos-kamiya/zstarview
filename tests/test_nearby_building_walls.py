from __future__ import annotations

import math
import random
from dataclasses import replace
from types import SimpleNamespace

from PySide6.QtGui import QImage, QPainter

from zstarview.astro import altaz_to_normalized_xy, is_in_fov
from zstarview.data.urban_outline_common import BuildingFootprint
from zstarview.nearby_walls import (
    WallSegment,
    visible_wall_segments,
    wall_opacity_level,
)
from zstarview.render import terrain
from zstarview.render.terrain import draw_urban_outlines
from zstarview.types import ScreenGeometry, UrbanOutlinePolyline, ViewerData
from zstarview.urban_outline_layer import (
    _nearby_wall_faces,
    _project_observer_centric_urban_outline_result,
)


def test_nearby_walls_use_footprints_and_exclude_observer_building() -> None:
    def building(name: str, west: float, east: float) -> BuildingFootprint:
        return BuildingFootprint(
            building_id=name,
            height_m=8.0,
            rings_lonlat=(
                (
                    (west, -0.00003),
                    (east, -0.00003),
                    (east, 0.00003),
                    (west, 0.00003),
                    (west, -0.00003),
                ),
            ),
        )

    nearby = building("near", 0.0001, 0.00015)
    too_close = building("too_close", 0.00002, 0.00003)
    inside = building("inside", -0.00002, 0.00002)
    near_40 = building("near_40", 0.00027, 0.00030)
    near_100 = building("near_100", 0.00070, 0.00075)
    near_200 = building("near_200", 0.0013, 0.00135)
    distant = building("far", 0.005, 0.0051)
    kwargs = dict(lat_deg=0.0, lon_deg=0.0, observer_elevation_m=1.7)
    walls = _nearby_wall_faces(
        (nearby, near_40, near_100, near_200, too_close, inside, distant), **kwargs
    )
    assert walls
    assert walls == _nearby_wall_faces((nearby, near_40, near_100, near_200), **kwargs)
    assert all(len(face.points) == 4 and face.distance_km < 0.2 for face in walls)
    outer_walls = _nearby_wall_faces((near_200,), **kwargs)
    assert outer_walls and all(face.distance_km > 0.1 for face in outer_walls)
    assert wall_opacity_level(5.0) == 0
    assert wall_opacity_level(7.5) == 16
    assert wall_opacity_level(10.0) == 32
    assert wall_opacity_level(255.0) == 18
    assert wall_opacity_level(500.0) == 3
    assert wall_opacity_level(500.01) == 0

    beyond_default = building("beyond_default", 0.0054, 0.00545)
    assert not _nearby_wall_faces((beyond_default,), **kwargs)
    extended = _nearby_wall_faces((beyond_default,), max_distance_m=750.0, **kwargs)
    assert extended and all(
        face.wall_opacity_level is not None and face.wall_opacity_level > 0
        for face in extended
    )
    assert not _nearby_wall_faces((nearby,), max_distance_m=0.0, **kwargs)
    assert wall_opacity_level(150.0, 500.0) > wall_opacity_level(150.0, 200.0)
    projected = _project_observer_centric_urban_outline_result(
        SimpleNamespace(outlines=(), walls=extended), view_center=(0.0, 90.0),
        edge_fov_deg=90.0,
    )
    assert projected and projected[0].wall_opacity_level == extended[0].wall_opacity_level
    assert projected[0].building_id == beyond_default.building_id


def test_wall_subdivision_uses_larger_steps_farther_away() -> None:
    ring = (
        (0.00017, -0.00003),
        (0.00023, -0.00003),
        (0.00023, 0.00003),
        (0.00017, 0.00003),
        (0.00017, -0.00003),
    )
    near = BuildingFootprint(building_id="wall", height_m=8.0, rings_lonlat=(ring,))
    far_ring = tuple((lon + 0.0002, lat) for lon, lat in ring)
    far = replace(near, rings_lonlat=(far_ring,))
    kwargs = dict(lat_deg=0.0, lon_deg=0.0, observer_elevation_m=1.7)

    near_faces = _nearby_wall_faces((near,), **kwargs)
    far_faces = _nearby_wall_faces((far,), **kwargs)

    # Three surviving walls, each split into 2 x 2 patches near and 1 x 1 far.
    assert len(near_faces) == 12
    assert len(far_faces) == 3


def test_part_wall_respects_min_height_and_containing_rings_are_excluded() -> None:
    ring = ((0.0002, -0.0001), (0.0003, -0.0001), (0.0003, 0.0001),
            (0.0002, 0.0001), (0.0002, -0.0001))
    part = BuildingFootprint(building_id="part", parent_building_id="parent", height_m=8.0,
                             min_height_m=4.0, rings_lonlat=(ring,))
    kwargs = dict(lat_deg=0.0, lon_deg=0.0, observer_elevation_m=0.0)
    raised = _nearby_wall_faces((part,), **kwargs)
    grounded = _nearby_wall_faces((replace(part, min_height_m=0.0),), **kwargs)
    assert raised and grounded
    assert min(p.altitude_deg for face in raised for p in face.points) > 0.0
    assert min(p.altitude_deg for face in grounded for p in face.points) == 0.0

    enclosing = ((-0.0004, -0.0004), (0.0004, -0.0004), (0.0004, 0.0004),
                 (-0.0004, 0.0004), (-0.0004, -0.0004))
    # Flattened courtyard or multipolygon rings containing the observer exclude the building.
    assert not _nearby_wall_faces((replace(part, rings_lonlat=(enclosing, ring)),), **kwargs)


def test_overlapping_wall_patches_are_composited_once() -> None:
    wall = UrbanOutlinePolyline(
        points=[(-5.0, 85.0), (-5.0, 95.0), (5.0, 95.0), (5.0, 85.0)],
        height_m=10.0,
        distance_km=0.01,
        source="nearby_wall",
        building_id="near",
    )
    viewer = ViewerData(
        location=(0.0, 0.0),
        timezone_name="UTC",
        city_name="test",
        view_center=(0.0, 90.0),
    )

    def render(
        walls: list[UrbanOutlinePolyline], *, scalar: bool = False, inverted: bool = False
    ) -> QImage:
        image = QImage(200, 200, QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(0)
        painter = QPainter(image)
        projection_args = {}
        if scalar:
            projection_args = {
                "is_in_fov_func": lambda *args, **kwargs: is_in_fov(*args, **kwargs),
                "altaz_to_normalized_xy_func": lambda *args, **kwargs: (
                    altaz_to_normalized_xy(*args, **kwargs)
                ),
            }
        draw_urban_outlines(
            painter,
            ScreenGeometry(center=(100, 100), radius=100),
            viewer,
            walls,
            opacity=0.2,
            inverted_city=inverted,
            split_by_gaps_func=lambda points: [points],
            **projection_args,
        )
        painter.end()
        return image

    once = render([wall])
    twice = render([wall, wall])
    far = replace(wall, distance_km=0.08)
    faded = render([far])
    wall_color = once.pixelColor(100, 100)
    assert wall_color.alpha() == round(255 * 0.2 * 0.5)
    # Premultiplied storage rounds the RGB channels at partial opacity.
    assert all(
        abs(actual - expected) <= 4
        for actual, expected in zip(wall_color.getRgb()[:3], (207, 229, 255))
    )
    assert 0 < faded.pixelColor(100, 100).alpha() < once.pixelColor(100, 100).alpha()
    assert once == twice
    assert render([wall, far]) == once
    assert render([far, wall]) == once
    assert render([wall, far], scalar=True) == once
    roof = replace(wall, source="base", points=[*wall.points, wall.points[0]])
    roof_only = render([roof])
    assert roof_only.pixelColor(100, 100).alpha() > 0
    for scalar in (False, True):
        outlined_wall = render([wall, roof], scalar=scalar)
        assert outlined_wall.pixelColor(100, 100) == wall_color
        assert outlined_wall != once
    assert render([wall, replace(roof, building_id="other")]).pixelColor(
        100, 100
    ).alpha() > wall_color.alpha()
    out_of_view_wall = replace(wall, points=[(alt, az + 180) for alt, az in wall.points])
    assert render([out_of_view_wall, roof]) == roof_only
    extended_level = wall_opacity_level(150.0, 500.0)
    extended_wall = replace(wall, distance_km=0.15, wall_opacity_level=extended_level)
    assert render([extended_wall]).pixelColor(100, 100).alpha() > 0

    below_horizon_wall = replace(wall, points=[(alt - 10.0, az) for alt, az in wall.points])
    below_horizon_roof = replace(
        roof, points=[*below_horizon_wall.points, below_horizon_wall.points[0]]
    )
    for scalar in (False, True):
        empty = render([], scalar=scalar, inverted=True)
        inverted_roof = render([below_horizon_roof], scalar=scalar, inverted=True)
        assert inverted_roof != empty
        assert render([below_horizon_wall], scalar=scalar, inverted=True) == empty
        assert render(
            [below_horizon_wall, below_horizon_roof], scalar=scalar, inverted=True
        ) == inverted_roof

    image = QImage(200, 200, QImage.Format.Format_ARGB32_Premultiplied)
    image.setDevicePixelRatio(2.0)
    image.fill(0)
    painter = QPainter(image)
    painter.translate(10.0, 0.0)
    draw_urban_outlines(
        painter, ScreenGeometry(center=(50, 50), radius=50), viewer, [wall], opacity=0.2
    )
    painter.end()
    assert image.pixelColor(120, 100).alpha() > 0
    assert image.pixelColor(80, 100).alpha() == 0


def test_wall_occlusion_requires_full_coverage_and_sufficient_opacity() -> None:
    cover = WallSegment((-3, 10), (3, 10), -5, 5, 10, 10, 2)
    target = WallSegment((-1, 30), (1, 30), -5, 5, 30, 10, 2)
    assert visible_wall_segments([target, cover]) == [cover]
    # Partial coverage, lower opacity, or a wall behind the target cannot cull it.
    partial = replace(cover, start=(0, 10))
    assert target in visible_wall_segments([target, partial])
    faint = replace(cover, distance_m=5.1)
    assert wall_opacity_level(faint.distance_m) < wall_opacity_level(target.distance_m)
    assert target in visible_wall_segments([target, faint])
    behind = replace(cover, start=(-3, 40), end=(3, 40))
    assert target in visible_wall_segments([target, behind])


def test_wall_occlusion_handles_azimuth_wrap() -> None:
    cover = WallSegment((-3, 10), (3, 10), -5, 5, 10, 10, 2)
    target = WallSegment((-1, 30), (1, 30), -5, 5, 30, 10, 2)
    assert visible_wall_segments([cover, target]) == [cover]


def test_occlusion_uses_configured_wall_radius_for_opacity() -> None:
    cover = WallSegment((-3, 10), (3, 10), -5, 5, 7, 10, 2)
    target = WallSegment((-1, 150), (1, 150), -5, 5, 150, 10, 2)
    assert visible_wall_segments([target, cover], 200) == [cover]
    assert target in visible_wall_segments([target, cover], 500)


def _reference_covered_by_wall(
    target: WallSegment, cover: WallSegment, max_distance_m: float
) -> bool:
    """Scalar form of the exact occlusion test used before vectorization."""

    def cross(a: tuple[float, float], b: tuple[float, float]) -> float:
        return a[0] * b[1] - a[1] * b[0]

    if wall_opacity_level(cover.distance_m, max_distance_m) < wall_opacity_level(
        target.distance_m, max_distance_m
    ):
        return False
    edge = (cover.end[0] - cover.start[0], cover.end[1] - cover.start[1])
    numerator = cross(cover.start, edge)
    for point in (target.start, target.end):
        denominator = cross(point, edge)
        if abs(denominator) < 1e-9:
            return False
        t = numerator / denominator
        u = cross(cover.start, point) / denominator
        if not (1e-7 < t < 1.0 - 1e-7 and 1e-7 < u < 1.0 - 1e-7):
            return False
        if not (cover.bottom + 1e-7 < target.bottom * t and target.top * t < cover.top - 1e-7):
            return False
    return True


def _reference_visible_wall_segments(
    segments: list[WallSegment], max_distance_m: float
) -> list[WallSegment]:
    bin_deg = 5.0
    bin_count = round(360.0 / bin_deg)

    def azimuth(point: tuple[float, float]) -> float:
        return math.degrees(math.atan2(point[0], point[1])) % 360.0

    bins: list[list[int]] = [[] for _ in range(bin_count)]
    for i, segment in enumerate(segments):
        a, b = azimuth(segment.start), azimuth(segment.end)
        delta = (b - a + 180.0) % 360.0 - 180.0
        start = a if delta >= 0.0 else b
        for index in range(
            math.floor(start / bin_deg), math.floor((start + abs(delta)) / bin_deg) + 1
        ):
            bins[index % bin_count].append(i)
    visible = []
    for i, target in enumerate(segments):
        midpoint = (
            (target.start[0] + target.end[0]) * 0.5,
            (target.start[1] + target.end[1]) * 0.5,
        )
        if not any(
            j != i and _reference_covered_by_wall(target, segments[j], max_distance_m)
            for j in bins[int(azimuth(midpoint) // bin_deg)]
        ):
            visible.append(target)
    return visible


def test_vectorized_wall_occlusion_matches_scalar_reference() -> None:
    rng = random.Random(20261127)
    segments: list[WallSegment] = []
    for _ in range(80):
        # Building-like rectangles around the observer, plus a few segments that
        # share endpoints or lie on rays through the observer.
        cx = rng.uniform(-150.0, 150.0)
        cy = rng.uniform(-150.0, 150.0)
        half_w, half_h = rng.uniform(1.0, 15.0), rng.uniform(1.0, 15.0)
        corners = [
            (cx - half_w, cy - half_h), (cx + half_w, cy - half_h),
            (cx + half_w, cy + half_h), (cx - half_w, cy + half_h),
        ]
        distance = max(5.0, math.hypot(cx, cy) - math.hypot(half_w, half_h))
        bottom = rng.choice((-1.7, -5.0, rng.uniform(-10.0, 0.0)))
        top = bottom + rng.uniform(3.0, 60.0)
        for start, end in zip(corners, corners[1:] + corners[:1]):
            segments.append(WallSegment(start, end, bottom, top, distance, top, 2.0))
    segments.append(WallSegment((0.0, 20.0), (0.0, 40.0), -2.0, 10.0, 20.0, 10.0, 2.0))
    segments.append(WallSegment((-5.0, 10.0), (5.0, 10.0), -2.0, 10.0, 10.0, 10.0, 2.0))
    segments.append(WallSegment((-5.0, 10.0), (5.0, 10.0), -2.0, 10.0, 10.0, 10.0, 2.0))
    for max_distance_m in (200.0, 500.0):
        expected = _reference_visible_wall_segments(segments, max_distance_m)
        assert 0 < len(expected) < len(segments)
        assert visible_wall_segments(segments, max_distance_m) == expected


def test_wall_occlusion_handles_empty_input() -> None:
    assert visible_wall_segments([]) == []


def test_wall_panorama_crosses_north_and_is_reused_across_view_changes(monkeypatch) -> None:
    wall = UrbanOutlinePolyline(
        points=[(-5.0, 355.0), (-5.0, 5.0), (5.0, 5.0), (5.0, 355.0)],
        height_m=10.0,
        distance_km=0.01,
        source="nearby_wall",
    )
    atlas_sizes: list[tuple[int, ...]] = []
    original = terrain._build_wall_panorama

    def track(*args, **kwargs):
        panorama = original(*args, **kwargs)
        atlas_sizes.append(panorama.levels.shape)
        return panorama

    monkeypatch.setattr(terrain, "_build_wall_panorama", track)
    walls = [wall]
    geometry = ScreenGeometry((100, 100), 100)
    viewer = ViewerData(
        location=(0.0, 0.0),
        timezone_name="UTC",
        city_name="test",
        view_center=(0.0, 359.0),
    )

    def paint(outlines: list[UrbanOutlinePolyline], view: ViewerData, size: int) -> int:
        image = QImage(size, size, QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(0)
        painter = QPainter(image)
        draw_urban_outlines(painter, geometry, view, outlines)
        painter.end()
        return image.pixelColor(100, 100).alpha()

    assert paint(walls, viewer, 200) > 0
    assert paint(walls, replace(viewer, view_center=(0.0, 1.0)), 240) > 0
    assert paint(walls, replace(viewer, view_center=(0.0, 180.0)), 200) == 0
    assert atlas_sizes == [(44, 1440)]
    assert paint(list(walls), viewer, 200) > 0
    assert atlas_sizes == [(44, 1440), (44, 1440)]


def test_wall_panorama_screen_crop_preserves_visible_pixels(monkeypatch) -> None:
    wall = UrbanOutlinePolyline(
        points=[(-2.0, 355.0), (-2.0, 5.0), (2.0, 5.0), (2.0, 355.0)],
        height_m=10.0,
        distance_km=0.01,
        source="nearby_wall",
    )
    outlines = [wall]
    geometry = ScreenGeometry((200, 200), 200)
    viewer = ViewerData(
        location=(0.0, 0.0), timezone_name="UTC", city_name="test",
        view_center=(0.0, 0.0),
    )

    def paint(altitude: float, size: int) -> QImage:
        terrain._project_wall_panorama.cache_clear()
        image = QImage(size, size, QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(0)
        painter = QPainter(image)
        draw_urban_outlines(
            painter, geometry, replace(viewer, view_center=(altitude, 0.0)), outlines
        )
        painter.end()
        return image

    for size in (400, 401):
        for altitude in (0.0, 45.0, -45.0, 89.0):
            cropped = paint(altitude, size)
            with monkeypatch.context() as patch:
                patch.setattr(
                    terrain, "_wall_panorama_screen_bounds",
                    lambda size, *args: (0, 0, *size),
                )
                full = paint(altitude, size)
            assert cropped == full

    bounds = terrain._wall_panorama_screen_bounds(
        (100, 100), (4.0, 4.0), geometry, (0.0, 0.0),
        90.0, 100.0, (-2.5, 2.5),
    )
    assert bounds is not None
    assert bounds[3] - bounds[1] < 50


def _square(lon0: float, lat0: float, half: float, z: float) -> tuple[tuple[float, float, float], ...]:
    corners = (
        (lon0 - half, lat0 - half),
        (lon0 + half, lat0 - half),
        (lon0 + half, lat0 + half),
        (lon0 - half, lat0 + half),
    )
    ring = tuple((lon, lat, z) for lon, lat in corners)
    return ring + (ring[0],)


def test_roof_floor_prisms_extrude_to_next_lower_floor() -> None:
    from zstarview.nearby_walls import roof_floor_prisms

    ground = 18.4
    base = _square(0.0, 0.0, 0.0003, ground + 6.2)
    # Not above `base` in plan: floors are ordered over the whole building.
    deck = _square(0.001, 0.0, 0.0001, ground + 129.8)
    sloped = (
        (0.0, 0.0, ground + 228.0),
        (0.00005, 0.0, ground + 230.0),
        (0.00005, 0.00005, ground + 230.0),
        (0.0, 0.00005, ground + 228.0),
        (0.0, 0.0, ground + 228.0),
    )
    same_floor = _square(0.002, 0.0, 0.0001, ground + 6.4)
    below_ground = _square(0.003, 0.0, 0.0001, ground - 1.0)

    prisms = roof_floor_prisms((deck, base, sloped, same_floor, below_ground), ground)

    assert [(p.bottom_m, p.top_m) for p in prisms] == [
        (6.0, 130.0),
        (0.0, 6.0),
        (130.0, 229.0),
        (0.0, 6.0),
    ]
    assert prisms[1].ring_lonlat[0] == prisms[1].ring_lonlat[-1]
    assert prisms[1].ring_lonlat[0] == base[0][:2]
    assert roof_floor_prisms((), ground) == ()


def _tower(roofs, ring_half: float = 0.0003) -> BuildingFootprint:
    footprint = tuple(p[:2] for p in _square(0.002, 0.0, ring_half, 0.0))
    return BuildingFootprint(
        building_id="tower",
        height_m=230.0,
        rings_lonlat=(footprint,),
        geometry_lod=2,
        roof_surfaces_lonlat=tuple(roofs),
        ground_elevation_m=10.0,
    )


def test_plateau_roof_surfaces_replace_full_height_footprint_walls() -> None:
    kwargs = dict(lat_deg=0.0, lon_deg=0.0, observer_elevation_m=11.7, max_distance_m=500.0)
    roofs = (
        _square(0.002, 0.0, 0.0003, 10.0 + 20.0),  # wide low floor
        _square(0.002, 0.0, 0.00005, 10.0 + 230.0),  # thin spire
    )
    tower = _tower(roofs)
    walls = _nearby_wall_faces((tower,), **kwargs)
    approximated = _nearby_wall_faces(
        (replace(tower, geometry_lod=1, roof_surfaces_lonlat=()),), **kwargs
    )
    assert walls and approximated

    def max_altitude_beyond(faces, azimuth_offset_deg: float) -> float:
        return max(
            (
                p.altitude_deg
                for face in faces
                for p in face.points
                if abs(p.azimuth_deg - 90.0) > azimuth_offset_deg
            ),
            default=-90.0,
        )

    # The wide floor stops at 20 m; only the thin spire reaches 230 m.
    spire_half_width_deg = 2.0
    assert max_altitude_beyond(walls, spire_half_width_deg) < 10.0
    assert max_altitude_beyond(approximated, spire_half_width_deg) > 45.0
    assert max(p.altitude_deg for face in walls for p in face.points) > 45.0


def test_roof_floor_prisms_skip_rings_around_observer_and_fall_back() -> None:
    kwargs = dict(lat_deg=0.0, lon_deg=0.0, observer_elevation_m=11.7, max_distance_m=500.0)
    # A roof polygon covering the observer must not wall in the observer.
    covering = _square(0.0, 0.0, 0.001, 10.0 + 20.0)
    spire = _square(0.002, 0.0, 0.00005, 10.0 + 230.0)
    walls = _nearby_wall_faces((_tower((covering, spire)),), **kwargs)
    assert walls
    assert min(p.altitude_deg for face in walls for p in face.points) > 1.0

    # Roofs at or below the ground leave no prisms; the footprint is used.
    buried = _tower((_square(0.002, 0.0, 0.0003, 5.0),))
    fallback = _nearby_wall_faces((buried,), **kwargs)
    expected = _nearby_wall_faces(
        (replace(buried, geometry_lod=1, roof_surfaces_lonlat=()),), **kwargs
    )
    assert fallback == expected
