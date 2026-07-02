import astropy.time
from PySide6.QtCore import QPointF
from PySide6.QtGui import QFont

from zstarview.aircraft.types import AircraftOverlayPoint
from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.render import aircraft as render_aircraft
from zstarview.types import ScreenGeometry, ViewerData


def test_aircraft_draw_overlay_fills_and_outlines_ribbon(monkeypatch) -> None:
    class _Painter:
        def __init__(self) -> None:
            self.draw_polygon_calls = 0

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def font(self) -> QFont:
            return QFont()

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolygon(self, *_args, **_kwargs) -> None:
            self.draw_polygon_calls += 1

        def drawPolyline(self, *_args, **_kwargs) -> None:
            raise AssertionError("drawPolyline should not be used when ribbon polygons exist")

    painter = _Painter()
    aircraft_points = [
        AircraftOverlayPoint(
            icao24="abc123",
            callsign="TEST123",
            alt_deg=10.0,
            az_deg=151.0,
            trail_alt_az_points=((10.0, 151.0), (10.2, 151.2), (10.4, 151.5)),
            distance_km=5.0,
            age_seconds=10.0,
            alpha_scale=1.0,
        )
    ]
    monkeypatch.setattr(
        render_aircraft,
        "project_aircraft_snapshots",
        lambda *_args, **_kwargs: aircraft_points,
    )

    render_aircraft.draw_aircraft_overlay(
        painter=painter,
        geometry=ScreenGeometry(center=(100, 100), radius=100),
        aircraft_snapshots={"dummy": True},
        viewer_data=ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(10.0, 151.0),
            edge_fov_deg=110.0,
            content_fov_deg=180.0,
            observer_height_m=1.7,
        ),
        time_obj=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
        opacity=1.0,
        label_candidates=[],
        theme=THEME_STYLES_BY_PRESET["day"],
    )

    assert painter.draw_polygon_calls >= 2


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
