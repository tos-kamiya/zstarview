from __future__ import annotations

from dataclasses import replace

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
from zstarview.urban_outline_layer import _nearby_wall_faces


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
                "altaz_to_normalized_xy_func": lambda *args, **kwargs: altaz_to_normalized_xy(*args, **kwargs),
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
    viewer = ViewerData(location=(0.0, 0.0), timezone_name="UTC", city_name="test",
                        view_center=(0.0, 90.0))

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
