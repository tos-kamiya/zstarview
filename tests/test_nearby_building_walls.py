from __future__ import annotations

from dataclasses import replace

from PySide6.QtGui import QImage, QPainter

from zstarview.data.urban_outline_common import BuildingFootprint
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
    distant = building("far", 0.001, 0.0011)
    kwargs = dict(lat_deg=0.0, lon_deg=0.0, observer_elevation_m=1.7)
    walls = _nearby_wall_faces((nearby, near_40, near_100, too_close, inside, distant), **kwargs)
    assert walls
    assert walls == _nearby_wall_faces((nearby, near_40, near_100), **kwargs)
    assert all(len(face.points) == 4 and face.distance_km < 0.1 for face in walls)
    assert any(face.distance_km > 0.04 for face in walls)


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

    def render(walls: list[UrbanOutlinePolyline]) -> QImage:
        image = QImage(200, 200, QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(0)
        painter = QPainter(image)
        draw_urban_outlines(
            painter,
            ScreenGeometry(center=(100, 100), radius=100),
            viewer,
            walls,
            opacity=0.2,
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
