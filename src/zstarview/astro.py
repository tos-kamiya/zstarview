import math
import time as _time
from pathlib import Path
from typing import Any, Dict, List, Tuple, cast

from appdirs import user_cache_dir
import astropy
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, GeocentricTrueEcliptic, SkyCoord
from skyfield.api import Loader, Topos
import skyfield.api

from .paths import (
    APP_AUTHOR,
    APP_ID,
    ANGLE_BELOW_HORIZON,
    FIELD_OF_VIEW_DEG,
    PLANET_SYMBOLS,
)
from .types import PlanetBody, StarData


# Skyfield ephemeris cache loader (separate from UI)
_cache_path = Path(user_cache_dir(appname=APP_ID, appauthor=APP_AUTHOR))
_cache_path.mkdir(parents=True, exist_ok=True)
_starfield_load = Loader(str(_cache_path))


def altaz_to_normalized_xy(alt: float, az: float, view_center: Tuple[float, float]) -> Tuple[float, float]:
    """Convert alt/az to normalized screen coordinates relative to a view center.

    Returns (nx, ny) where 1.0 equals 90 degrees of angular distance.
    """
    center_alt, center_az = view_center
    alt1 = math.radians(center_alt)
    az1 = math.radians(center_az)
    alt2 = math.radians(alt)
    az2 = math.radians(az)

    cos_theta = math.sin(alt1) * math.sin(alt2) + math.cos(alt1) * math.cos(alt2) * math.cos(az2 - az1)
    theta = math.acos(max(-1.0, min(1.0, cos_theta)))

    r = theta / (math.pi / 2)

    dx = math.cos(alt2) * math.sin(az2 - az1)
    dy = math.cos(alt1) * math.sin(alt2) - math.sin(alt1) * math.cos(alt2) * math.cos(az2 - az1)
    length = math.hypot(dx, dy)
    if length != 0:
        dx /= length
        dy /= length
    nx = r * dx
    ny = -r * dy
    return (nx, ny)


def is_in_fov(alt: float, az: float, view_center: Tuple[float, float]) -> bool:
    """Check if a target at (alt, az) is within the field of view relative to view_center."""
    center_alt, center_az = view_center
    alt1, az1 = math.radians(center_alt), math.radians(center_az)
    alt2, az2 = math.radians(alt), math.radians(az)
    cos_theta = math.sin(alt1) * math.sin(alt2) + math.cos(alt1) * math.cos(alt2) * math.cos(az2 - az1)
    theta = math.acos(min(1.0, max(-1.0, cos_theta)))
    return math.degrees(theta) <= FIELD_OF_VIEW_DEG / 2


def load_star_catalog(filename: str) -> List[Dict[str, Any]]:
    # Re-exported from catalog; avoid import cycles if needed elsewhere
    from .catalog import load_star_catalog as _load

    return _load(filename)


def calculate_visible_stars(
    star_catalog: List[Dict[str, Any]],
    lat: float,
    lon: float,
    time_obj: astropy.time.Time,
    view_center: Tuple[float, float],
) -> Tuple[List[StarData], EarthLocation]:
    """Compute visible stars and return them with the observer location."""
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    visible_stars: List[StarData] = []
    for i, star in enumerate(star_catalog):
        coord = cast(SkyCoord, star["coord"])
        altaz = coord.transform_to(AltAz(obstime=time_obj, location=location))
        if altaz.alt.deg > -ANGLE_BELOW_HORIZON and is_in_fov(altaz.alt.deg, altaz.az.deg, view_center):
            visible_stars.append(StarData(name=star["name"], alt=altaz.alt.deg, az=altaz.az.deg, vmag=star["vmag"], bv=star["bv"]))
        if (i + 1) % 500 == 0:
            _time.sleep(0)
    return (visible_stars, location)


def calculate_moon_phase_angle(observer: Topos, t: skyfield.timelib.Time, planets: Any) -> float:
    """Calculate the phase angle of the Moon (Sun-Moon separation at observer)."""
    moon = planets["moon"]
    sun = planets["sun"]
    e_to_moon = observer.at(t).observe(moon).apparent()
    e_to_sun = observer.at(t).observe(sun).apparent()
    return e_to_moon.separation_from(e_to_sun).degrees


def calculate_visible_planets(
    lat: float,
    lon: float,
    astropy_time: astropy.time.Time,
    view_center: Tuple[float, float],
) -> List[PlanetBody]:
    """Calculate visible bodies (Sun, Moon, planets)."""
    ts = skyfield.api.load.timescale()
    t = ts.from_astropy(astropy_time)
    planets = _starfield_load("de421.bsp")
    observer = planets["earth"] + Topos(latitude_degrees=lat, longitude_degrees=lon)

    visible_bodies: List[PlanetBody] = []
    for name, symbol in PLANET_SYMBOLS.items():
        planet = planets[name]
        astrometric = observer.at(t).observe(planet).apparent()
        alt, az, _ = astrometric.altaz()
        if alt.degrees > -ANGLE_BELOW_HORIZON and is_in_fov(alt.degrees, az.degrees, view_center):
            phase_angle = None
            if name == "moon":
                phase_angle = calculate_moon_phase_angle(observer, t, planets)
            visible_bodies.append(
                PlanetBody(
                    name=name,
                    alt=alt.degrees,
                    az=az.degrees,
                    symbol=symbol,
                    phase_angle=phase_angle,
                )
            )
    return visible_bodies


def calculate_horizon_points(location: EarthLocation, time: astropy.time.Time, view_center: Tuple[float, float]) -> List[Tuple[float, float]]:
    """Generate points along the horizon for drawing."""
    points: List[Tuple[float, float]] = []
    alt = 0.0
    for az in range(0, 360 + 5, 5):
        if not is_in_fov(alt, az, view_center):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        points.append((nx, ny))
    return points


def calculate_celestial_equator_points(location: EarthLocation, time: astropy.time.Time, view_center: Tuple[float, float]) -> List[Tuple[float, float]]:
    """Generate points along the celestial equator for drawing."""
    a = 5
    points: List[Tuple[float, float]] = []
    for ra_deg in range(0, 360 + 5, 5):
        coord = SkyCoord(ra=ra_deg * u.deg, dec=0 * u.deg, frame="icrs")
        altaz = coord.transform_to(AltAz(obstime=time, location=location))
        if altaz.alt.deg > -a and is_in_fov(altaz.alt.deg, altaz.az.deg, view_center):
            nx, ny = altaz_to_normalized_xy(altaz.alt.deg, altaz.az.deg, view_center)
            points.append((nx, ny))
    return points


def calculate_ecliptic_points(location: EarthLocation, time: astropy.time.Time, view_center: Tuple[float, float]) -> List[Tuple[float, float]]:
    """Generate points along the ecliptic for drawing."""
    points: List[Tuple[float, float]] = []
    for lon_deg in range(0, 360 + 5, 5):
        ecl = SkyCoord(lon=lon_deg * u.deg, lat=0 * u.deg, frame=GeocentricTrueEcliptic(obstime=time))
        icrs = ecl.transform_to("icrs")
        altaz = icrs.transform_to(AltAz(obstime=time, location=location))
        if altaz.alt.deg > -ANGLE_BELOW_HORIZON and is_in_fov(altaz.alt.deg, altaz.az.deg, view_center):
            nx, ny = altaz_to_normalized_xy(altaz.alt.deg, altaz.az.deg, view_center)
            points.append((nx, ny))
    return points


def calculate_sun_angle_on_moon(moon_altaz: Tuple[float, float], sun_altaz: Tuple[float, float]) -> float:
    """Angle of the sun relative to the moon on screen (radians)."""
    m_alt, m_az = moon_altaz
    s_alt, s_az = sun_altaz
    d_az = math.radians(s_az - m_az)
    d_alt = math.radians(s_alt - m_alt)
    dx = -math.sin(d_az) * math.cos(math.radians(s_alt))
    dy = math.sin(d_alt)
    angle = math.atan2(dx, dy)
    return angle
