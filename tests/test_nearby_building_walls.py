from __future__ import annotations

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
    distant = building("far", 0.002, 0.0021)
    kwargs = dict(lat_deg=0.0, lon_deg=0.0, observer_elevation_m=1.7)
    walls = _nearby_wall_faces(
        (nearby, near_40, near_100, near_200, too_close, inside, distant), **kwargs
    )
    assert walls
    assert walls == _nearby_wall_faces((nearby, near_40, near_100, near_200), **kwargs)
    assert all(len(face.points) == 4 and face.distance_km < 0.2 for face in walls)
    outer_walls = _nearby_wall_faces((near_200,), **kwargs)
    assert outer_walls and all(face.distance_km > 0.1 for face in outer_walls)
    assert wall_opacity_level(100.0) > wall_opacity_level(150.0) > 0
    assert wall_opacity_level(200.0) == 0

    beyond_default = building("beyond_default", 0.0027, 0.00275)
    assert not _nearby_wall_faces((beyond_default,), **kwargs)
    extended = _nearby_wall_faces((beyond_default,), max_distance_m=500.0, **kwargs)
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
    )
    viewer = ViewerData(
        location=(0.0, 0.0),
        timezone_name="UTC",
        city_name="test",
        view_center=(0.0, 90.0),
    )

    def render(walls: list[UrbanOutlinePolyline], *, scalar: bool = False) -> QImage:
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
            **projection_args,
        )
        painter.end()
        return image

    once = render([wall])
    twice = render([wall, wall])
    far = replace(wall, distance_km=0.08)
    faded = render([far])
    assert once.pixelColor(100, 100).alpha() > 0
    assert 0 < faded.pixelColor(100, 100).alpha() < once.pixelColor(100, 100).alpha()
    assert once == twice
    assert render([wall, far]) == once
    assert render([far, wall]) == once
    assert render([wall, far], scalar=True) == once
    extended_level = wall_opacity_level(150.0, 500.0)
    extended_wall = replace(wall, distance_km=0.15, wall_opacity_level=extended_level)
    assert render([extended_wall]).pixelColor(100, 100).alpha() > 0

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


def test_wall_projection_cache_tracks_scene_and_view(monkeypatch) -> None:
    wall = UrbanOutlinePolyline(
        points=[(-5.0, 85.0), (-5.0, 95.0), (5.0, 95.0), (5.0, 85.0)],
        height_m=10.0,
        distance_km=0.01,
        source="nearby_wall",
    )
    calls: list[tuple[float, float]] = []
    original = terrain._nearby_wall_paths_vectorized

    def track(*args, **kwargs):
        calls.append(kwargs["view_center"])
        return original(*args, **kwargs)

    monkeypatch.setattr(terrain, "_nearby_wall_paths_vectorized", track)
    walls = [wall]
    geometry = ScreenGeometry((100, 100), 100)
    viewer = ViewerData(
        location=(0.0, 0.0),
        timezone_name="UTC",
        city_name="test",
        view_center=(0.0, 90.0),
    )

    def paint(outlines: list[UrbanOutlinePolyline], view: ViewerData) -> None:
        image = QImage(200, 200, QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(0)
        painter = QPainter(image)
        draw_urban_outlines(painter, geometry, view, outlines)
        painter.end()

    paint(walls, viewer)
    paint(walls, viewer)
    assert len(calls) == 1
    paint(walls, replace(viewer, view_center=(0.0, 91.0)))
    paint(list(walls), viewer)
    assert len(calls) == 3


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
