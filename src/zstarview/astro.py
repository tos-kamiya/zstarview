import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from appdirs import user_cache_dir
import astropy
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, GeocentricTrueEcliptic, SkyCoord
from skyfield.api import Loader, Topos
import skyfield.api
import numpy as np
import polars as pl

from .paths import (
    APP_AUTHOR,
    APP_ID,
    ANGLE_BELOW_HORIZON,
    FIELD_OF_VIEW_DEG,
    PLANET_SYMBOLS,
    PLANET_IDS,
)
from .types import LunarEclipseInfo, PlanetBody, SolarEclipseInfo


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


def is_in_fov_vectorized(alt: np.ndarray, az: np.ndarray, view_center: Tuple[float, float]) -> np.ndarray:
    """Vectorized check if targets at (alt, az) are within the field of view."""
    center_alt, center_az = view_center
    alt1, az1 = np.radians(center_alt), np.radians(center_az)
    alt2, az2 = np.radians(alt), np.radians(az)

    cos_theta = np.sin(alt1) * np.sin(alt2) + np.cos(alt1) * np.cos(alt2) * np.cos(az2 - az1)

    # Clip to avoid domain errors with arccos
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)

    return np.degrees(theta) <= FIELD_OF_VIEW_DEG / 2


def load_star_catalog(filename: str) -> List[Dict[str, Any]]:
    """Load the star catalog from file.

    This function re-exports the implementation from `catalog` to
    avoid circular imports when used elsewhere.
    """
    from .catalog import load_star_catalog as _load

    return _load(filename)


def calculate_visible_stars(
    star_df: pl.DataFrame,
    lat: float,
    lon: float,
    time_obj: astropy.time.Time,
    view_center: Tuple[float, float],
) -> Tuple[Dict[str, np.ndarray], EarthLocation]:
    """Compute visible stars and return them with the observer location."""
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    altaz_frame = AltAz(obstime=time_obj, location=location)

    # Get data as numpy arrays, ensure numeric types
    ra_h = star_df["RAh"].cast(pl.Float64, strict=False).to_numpy()
    dec = star_df["Dec"].cast(pl.Float64, strict=False).to_numpy()
    vmag = star_df["Vmag"].cast(pl.Float64, strict=False).to_numpy()
    # B-V can be NaN if empty
    bv = star_df["B-V"].cast(pl.Float64, strict=False).fill_null(np.nan).to_numpy()
    name = star_df["Name"].to_numpy()

    # Create a single SkyCoord object for all stars
    coords = SkyCoord(ra=(ra_h * 15.0) * u.deg, dec=dec * u.deg, frame="icrs")

    # Transform all coordinates at once
    altaz_coords = coords.transform_to(altaz_frame)
    alt = altaz_coords.alt.deg
    az = altaz_coords.az.deg

    # Vectorized visibility check
    in_view_mask = (alt > -ANGLE_BELOW_HORIZON) & is_in_fov_vectorized(alt, az, view_center)

    # Filter the results using the boolean mask
    visible_stars = {
        "name": name[in_view_mask],
        "alt": alt[in_view_mask],
        "az": az[in_view_mask],
        "vmag": vmag[in_view_mask],
        "bv": bv[in_view_mask],
    }

    return (visible_stars, location)


def calculate_moon_phase_angle(observer: Topos, t: skyfield.timelib.Time, planets: Any) -> float:
    """Calculate the phase angle of the Moon (Sun-Moon separation at observer)."""
    moon = planets["moon"]
    sun = planets["sun"]
    e_to_moon = observer.at(t).observe(moon).apparent()
    e_to_sun = observer.at(t).observe(sun).apparent()
    return e_to_moon.separation_from(e_to_sun).degrees


def calculate_lunar_eclipse_data(t: astropy.time.Time, observer) -> LunarEclipseInfo:
    """Calculate lunar eclipse parameters for a given time and observer.

    Returns an `EclipseInfo` object containing:
    - whether an eclipse occurs,
    - the eclipse type (penumbral, partial, total),
    - shadow center alt/az,
    - umbra and penumbra angular radii,
    - apparent moon radius.
    """
    # Load planetary ephemerides
    planets = _starfield_load("de440s.bsp")
    earth = planets["earth"]
    sun = planets["sun"]
    moon = planets["moon"]

    # Sun and Moon positions from Earth's center in GCRS
    sun_pos = earth.at(t).observe(sun)
    moon_pos = earth.at(t).observe(moon)

    # Exclude if Sun-Moon separation is not close to 180°
    separation_deg = sun_pos.separation_from(moon_pos).degrees
    if abs(separation_deg - 180.0) > 3.0:
        return LunarEclipseInfo(is_eclipse=False)

    # Astronomical constants
    R_earth_km = 6371.0
    R_sun_km = 696340.0
    R_moon_km = 1737.4

    # Distances from Earth's center
    D_sun_km = sun_pos.distance().km
    D_moon_km = moon_pos.distance().km

    # Angular radii
    earth_ang_rad = math.asin(R_earth_km / D_moon_km)
    sun_ang_rad = math.asin(R_sun_km / D_sun_km)
    moon_ang_rad = math.asin(R_moon_km / D_moon_km)

    # Umbra and penumbra radii
    umbra_radius_rad = max(0.0, earth_ang_rad - sun_ang_rad)
    penumbra_radius_rad = earth_ang_rad + sun_ang_rad

    # Angular distance between shadow center and Moon
    d_rad = math.radians(180.0 - separation_deg)

    # Classify eclipse type
    if d_rad > (penumbra_radius_rad + moon_ang_rad):
        eclipse_type = "none"
    elif d_rad > (umbra_radius_rad + moon_ang_rad):
        eclipse_type = "penumbral"
    elif d_rad > abs(umbra_radius_rad - moon_ang_rad):
        eclipse_type = "partial"
    else:
        eclipse_type = "total"

    # Shadow center as seen from the observer (opposite Sun direction)
    sun_apparent = observer.at(t).observe(sun).apparent()
    s_alt, s_az, _ = sun_apparent.altaz()
    shadow_center_alt = -s_alt.degrees
    shadow_center_az = (s_az.degrees + 180.0) % 360.0

    return LunarEclipseInfo(
        is_eclipse=(eclipse_type != "none"),
        eclipse_type=eclipse_type,
        shadow_center_alt=shadow_center_alt,
        shadow_center_az=shadow_center_az,
        umbra_radius_deg=math.degrees(umbra_radius_rad),
        penumbra_radius_deg=math.degrees(penumbra_radius_rad),
        moon_radius_deg=math.degrees(moon_ang_rad),
    )


SUN_RADIUS_KM = 695700.0
MOON_RADIUS_KM = 1737.4


def _angular_radius_rad(radius_km: float, distance_km: float) -> float:
    """Apparent angular radius [rad] = asin(R / D) (more precise than small-angle approximation)."""
    x = radius_km / max(1e-9, distance_km)
    return math.asin(max(-1.0, min(1.0, x)))


def _circle_overlap_area_fraction(R: float, r: float, d: float) -> float:
    """
    Calculates the overlapping area of two circles (radii R, r; center distance d)
    divided by the area of the solar disk (πR^2).
    R, r, and d are in the same units (here, angular radius in radians).
    """
    # Separated
    if d >= R + r:
        return 0.0
    # Contained (the smaller circle is completely inside the larger one)
    if d <= abs(R - r):
        return (min(R, r) ** 2) / (R**2)

    R2, r2 = R * R, r * r
    # Angles (segment angles)
    alpha = math.acos((d * d + R2 - r2) / (2.0 * d * R))
    beta = math.acos((d * d + r2 - R2) / (2.0 * d * r))
    # Lens area (Brahmagupta's formula-like square root term)
    lens = R2 * alpha + r2 * beta - 0.5 * math.sqrt(max(0.0, (-d + R + r) * (d + R - r) * (d - R + r) * (d + R + r)))
    return lens / (math.pi * R2)


def calculate_solar_eclipse_data(t: astropy.time.Time, observer) -> SolarEclipseInfo:
    planets = _starfield_load("de440s.bsp")
    sun = planets["sun"]
    moon = planets["moon"]

    # Observer's topocentric apparent position
    obs = observer.at(t)
    sun_app = obs.observe(sun).apparent()
    moon_app = obs.observe(moon).apparent()

    # Angular distance between centers
    sep = sun_app.separation_from(moon_app)
    sep_deg = sep.degrees
    sep_rad = sep.radians

    # Apparent angular radius (calculated from distance each time)
    D_sun_km = sun_app.distance().km
    D_moon_km = moon_app.distance().km
    sun_rad = _angular_radius_rad(SUN_RADIUS_KM, D_sun_km)  # [rad]
    moon_rad = _angular_radius_rad(MOON_RADIUS_KM, D_moon_km)  # [rad]

    # Obscuration (area fraction of the Sun obscured by the Moon)
    obscuration = _circle_overlap_area_fraction(sun_rad, moon_rad, sep_rad)

    # Determine eclipse type
    if sep_rad >= sun_rad + moon_rad or obscuration <= 0.0:
        eclipse_type = "none"
    elif sep_rad <= abs(sun_rad - moon_rad):
        # The smaller circle is completely contained: total if Moon > Sun, annular if Sun > Moon
        eclipse_type = "total" if moon_rad >= sun_rad else "annular"
    else:
        eclipse_type = "partial"

    return SolarEclipseInfo(
        is_eclipse=(eclipse_type != "none"),
        eclipse_type=eclipse_type,
        sep_deg=float(sep_deg),
        obscuration=obscuration,
    )


def eclipse_factor_from_info(info: Optional[SolarEclipseInfo]) -> float:
    """Return dimming factor in [0,1]. 1=no eclipse (no dim)."""
    if not info or not info.is_eclipse or info.obscuration <= 1e-3:
        return 1.0
    obsc = max(0.0, min(1.0, info.obscuration))

    # Perceptual model: gets dark suddenly after ~90%
    f0 = 0.92
    is_total = info.eclipse_type == "total"
    k = 10.0 if is_total else 7.5

    s = 1.0 / (1.0 + math.exp(k * (obsc - f0)))  # 1→0
    # Totality is darker
    min_l = 0.02 if is_total else 0.15
    return min_l + (1.0 - min_l) * s


def calculate_planets(
    lat: float,
    lon: float,
    astropy_time: astropy.time.Time,
    view_center: Tuple[float, float],
) -> List[PlanetBody]:
    """Calculate all planetary bodies (Sun, Moon, planets)."""
    ts = skyfield.api.load.timescale()
    t = ts.from_astropy(astropy_time)
    planets = _starfield_load("de440s.bsp")
    observer = planets["earth"] + Topos(latitude_degrees=lat, longitude_degrees=lon)

    bodies: List[PlanetBody] = []
    for name, symbol in PLANET_SYMBOLS.items():
        planet_id = PLANET_IDS[name]
        planet = planets[planet_id]
        astrometric = observer.at(t).observe(planet).apparent()
        alt, az, _ = astrometric.altaz()
        is_visible = alt.degrees > -ANGLE_BELOW_HORIZON and is_in_fov(alt.degrees, az.degrees, view_center)

        pa = lei = None
        if name == "moon":
            pa = calculate_moon_phase_angle(observer, t, planets)
            lei = calculate_lunar_eclipse_data(t, observer)
        sei = None
        if name == "sun":
            sei = calculate_solar_eclipse_data(t, observer)

        bodies.append(
            PlanetBody(
                name=name,
                alt=alt.degrees,
                az=az.degrees,
                symbol=symbol,
                is_visible=is_visible,
                phase_angle=pa,
                lunar_eclipse_info=lei,
                solar_eclipse_info=sei,
            )
        )
    return bodies


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


def calculate_moon_render_data(
    sun_altaz: Optional[Tuple[float, float]],
    moon_altaz: Optional[Tuple[float, float]],
    view_center: Tuple[float, float],
) -> Tuple[np.ndarray, float]:
    """
    Calculates all necessary data for rendering the moon.

    Returns a tuple of:
    - sun_dir_in_moon_frame: 3D vector of the sun's direction in the moon's local frame.
    - screen_rotation_deg: The angle to rotate the moon image on the screen.
    """

    # Helper to convert Alt/Az to a 3D vector in the observer's reference frame.
    # Y-up, X-East, Z-North (Right-handed system)
    def altaz_to_cartesian(alt: float, az: float) -> np.ndarray:
        alt_rad = math.radians(alt)
        az_rad = math.radians(az)  # Azimuth is measured from North (0) towards East (90)
        x = math.cos(alt_rad) * math.sin(az_rad)
        y = math.sin(alt_rad)
        z = math.cos(alt_rad) * math.cos(az_rad)
        return np.array([x, y, z])

    sun_dir_in_observer_frame = None
    if sun_altaz:
        sun_dir_in_observer_frame = altaz_to_cartesian(sun_altaz[0], sun_altaz[1])

    sun_dir_in_moon_frame = np.array([0.0, 0.0, 0.0])
    if sun_dir_in_observer_frame is not None and moon_altaz is not None:
        m_alt, m_az = moon_altaz

        # If the moon is very close to zenith or nadir (alt ~ +/-90 deg),
        # the calculation for the local coordinate frame can become unstable.
        # To avoid this singularity, we perturb the altitude slightly.
        if abs(m_alt) > 89.9999:
            m_alt = math.copysign(89.9999, m_alt)

        moon_vec = altaz_to_cartesian(m_alt, m_az)
        z_axis = -moon_vec / np.linalg.norm(moon_vec)
        zenith_vec = np.array([0.0, 1.0, 0.0])
        y_axis = zenith_vec - np.dot(zenith_vec, z_axis) * z_axis
        y_axis /= np.linalg.norm(y_axis)
        x_axis = np.cross(y_axis, z_axis)

        sun_dir_in_moon_frame = np.array(
            [
                np.dot(sun_dir_in_observer_frame, x_axis),
                np.dot(sun_dir_in_observer_frame, y_axis),
                np.dot(sun_dir_in_observer_frame, z_axis),
            ]
        )

    # Calculate the screen rotation for the moon image
    screen_rotation_deg = 0
    if moon_altaz:
        # Get screen coordinates for the moon's center and a point just "above" it
        m_alt, m_az = moon_altaz
        nx_center, ny_center = altaz_to_normalized_xy(m_alt, m_az, view_center)
        nx_up, ny_up = altaz_to_normalized_xy(m_alt + 0.1, m_az, view_center)

        # Calculate the angle of the "up" direction on the screen
        dx = nx_up - nx_center
        dy = ny_up - ny_center
        if abs(dx) > 1e-6 or abs(dy) > 1e-6:
            screen_up_angle_rad = math.atan2(dy, dx)
            # The image's "up" is vertical (points to -Y), so rotate to match the screen's "up"
            # After user feedback, the correct rotation is to flip the sign.
            screen_rotation_deg = -(math.degrees(screen_up_angle_rad) - 90)

    return sun_dir_in_moon_frame, screen_rotation_deg
