from PySide6.QtCore import QPointF

from zstarview.render import aircraft as render_aircraft
from zstarview.types import ScreenGeometry


def test_aircraft_ribbon_polygons_expand_a_linear_trail() -> None:
    geometry = ScreenGeometry(center=(100, 100), radius=100)
    screen_points = [
        QPointF(20.0, 40.0),
        QPointF(60.0, 40.0),
        QPointF(100.0, 40.0),
    ]

    polygons = render_aircraft._aircraft_ribbon_polygons(
        screen_points,
        ribbon_width_px=6.0,
        geometry=geometry,
    )

    assert len(polygons) == 1
    bounds = polygons[0].boundingRect()
    assert bounds.width() > 0.0
    assert bounds.height() > 0.0


def test_aircraft_ribbon_polygons_split_on_large_gaps() -> None:
    geometry = ScreenGeometry(center=(100, 100), radius=100)
    screen_points = [
        QPointF(20.0, 40.0),
        QPointF(60.0, 40.0),
        QPointF(320.0, 40.0),
        QPointF(360.0, 40.0),
    ]

    polygons = render_aircraft._aircraft_ribbon_polygons(
        screen_points,
        ribbon_width_px=6.0,
        geometry=geometry,
    )

    assert len(polygons) == 2
