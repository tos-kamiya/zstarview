from __future__ import annotations

import math
from datetime import datetime, timezone

import astropy.time
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, get_sun
from PySide6.QtCore import QPointF
from PySide6.QtGui import QColor, QImage, QPainter

from zstarview.astro import (
    altaz_to_normalized_xy,
    calculate_moon_north_up_screen_rotation,
    calculate_moon_render_data,
    calculate_solar_north_up_screen_rotation,
)
from zstarview.moon_hover import MoonHoverImage
from zstarview.render.solar_system import draw_nasa_moon_image


def _signed_delta_deg(new_angle: float, old_angle: float) -> float:
    return (float(new_angle) - float(old_angle) + 180.0) % 360.0 - 180.0


def test_moon_screen_rotation_follows_zenith_view_azimuth() -> None:
    _sun_dir, base_rotation = calculate_moon_render_data(
        (20.0, 100.0),
        (45.0, 180.0),
        (90.0, 0.0),
    )
    _sun_dir, rotated_view_rotation = calculate_moon_render_data(
        (20.0, 100.0),
        (45.0, 180.0),
        (90.0, 45.0),
    )

    assert math.isclose(
        _signed_delta_deg(rotated_view_rotation, base_rotation),
        45.0,
        abs_tol=1.0e-6,
    )


def test_moon_screen_rotation_stays_stable_for_same_view() -> None:
    _sun_dir, rotation_a = calculate_moon_render_data(
        (15.0, 80.0),
        (40.0, 170.0),
        (60.0, 120.0),
    )
    _sun_dir, rotation_b = calculate_moon_render_data(
        (15.0, 80.0),
        (40.0, 170.0),
        (60.0, 120.0),
    )

    assert rotation_a == rotation_b


def test_north_up_rotation_differs_from_zenith_up_when_needed() -> None:
    _sun_dir, local_rotation = calculate_moon_render_data(
        None,
        (45.0, 90.0),
        (45.0, 180.0),
    )
    north_up_rotation = calculate_moon_north_up_screen_rotation(
        (45.0, 90.0),
        (45.0, 180.0),
        observer_latitude_deg=35.0,
    )
    assert north_up_rotation != local_rotation


def test_north_up_rotation_maps_image_top_to_zenith_at_north_pole() -> None:
    moon_altaz = (45.0, 90.0)
    view_center = (45.0, 180.0)
    north_up_rotation = calculate_moon_north_up_screen_rotation(
        moon_altaz,
        view_center,
        observer_latitude_deg=90.0,
    )
    moon_x, moon_y = altaz_to_normalized_xy(*moon_altaz, view_center)
    north_x, north_y = altaz_to_normalized_xy(45.01, 90.0, view_center)
    screen_north_angle = math.degrees(math.atan2(north_y - moon_y, north_x - moon_x))
    rendered_image_top_angle = north_up_rotation - 90.0

    assert math.isclose(
        _signed_delta_deg(rendered_image_top_angle, screen_north_angle),
        0.0,
        abs_tol=0.05,
    )


def test_north_up_rotation_depends_on_observer_latitude() -> None:
    equator_rotation = calculate_moon_north_up_screen_rotation(
        (30.0, 120.0),
        (45.0, 180.0),
        observer_latitude_deg=0.0,
    )
    mid_latitude_rotation = calculate_moon_north_up_screen_rotation(
        (30.0, 120.0),
        (45.0, 180.0),
        observer_latitude_deg=35.0,
    )

    assert not math.isclose(equator_rotation, mid_latitude_rotation, abs_tol=1.0e-6)


def test_solar_north_rotation_changes_with_season() -> None:
    location = EarthLocation(lat=35.0 * u.deg, lon=139.0 * u.deg)
    march_time = astropy.time.Time("2024-03-20T12:00:00", scale="utc")
    september_time = astropy.time.Time("2024-09-22T12:00:00", scale="utc")
    march_sun = get_sun(march_time).transform_to(
        AltAz(obstime=march_time, location=location)
    )
    september_sun = get_sun(september_time).transform_to(
        AltAz(obstime=september_time, location=location)
    )
    march_altaz = (float(march_sun.alt.deg), float(march_sun.az.deg))
    september_altaz = (float(september_sun.alt.deg), float(september_sun.az.deg))
    march_rotation = calculate_solar_north_up_screen_rotation(
        march_altaz,
        march_altaz,
        time_obj=march_time,
        observer_latitude_deg=35.0,
        observer_longitude_deg=139.0,
    )
    september_rotation = calculate_solar_north_up_screen_rotation(
        september_altaz,
        september_altaz,
        time_obj=september_time,
        observer_latitude_deg=35.0,
        observer_longitude_deg=139.0,
    )
    assert abs(_signed_delta_deg(september_rotation, march_rotation)) > 40.0


def test_solar_north_rotation_follows_view_projection() -> None:
    kwargs = {
        "time_obj": astropy.time.Time("2024-03-20T12:00:00", scale="utc"),
        "observer_latitude_deg": 35.0,
        "observer_longitude_deg": 139.0,
    }
    base = calculate_solar_north_up_screen_rotation(
        (45.0, 180.0),
        (90.0, 0.0),
        **kwargs,
    )
    rotated = calculate_solar_north_up_screen_rotation(
        (45.0, 180.0),
        (90.0, 45.0),
        **kwargs,
    )
    assert math.isclose(_signed_delta_deg(rotated, base), 45.0, abs_tol=0.05)


def test_nasa_moon_image_masks_black_canvas() -> None:
    source = QImage(730, 730, QImage.Format.Format_RGB32)
    source.fill(QColor(255, 255, 255))
    image_data = MoonHoverImage(
        image=source,
        time_utc=datetime(2026, 8, 28, 4, tzinfo=timezone.utc),
        diameter_arcsec=1835.7,
    )
    target = QImage(120, 120, QImage.Format.Format_ARGB32)
    target.fill(QColor(0, 0, 0, 0))
    painter = QPainter(target)
    draw_nasa_moon_image(painter, QPointF(60.0, 60.0), 40.0, image_data, 0.0)
    painter.end()

    assert target.pixelColor(0, 0).alpha() == 0
    assert target.pixelColor(60, 60).alpha() > 0
