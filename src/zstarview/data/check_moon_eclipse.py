"""Quick verifier for lunar-eclipse parameters around a given time window.

This script prints eclipse status and basic geometry for Tokyo at 10-minute
intervals. It serves as a lightweight check of the core math used in the app.
"""

import math
from datetime import datetime, timedelta, timezone
from typing import Optional

from skyfield.api import Topos, load
from skyfield.timelib import Time


class EclipseInfo:
    """Minimal container for computed lunar-eclipse geometry.

    Attributes:
        is_eclipse: Whether an eclipse is happening (including penumbral).
        eclipse_type: "penumbral", "partial", "total", or "none".
        shadow_center_alt: Apparent altitude (deg) of umbra center (opposite Sun).
        shadow_center_az: Apparent azimuth (deg) of umbra center (opposite Sun).
        umbra_radius_deg: Angular radius (deg) of Earth's umbra at Moon distance.
        penumbra_radius_deg: Angular radius (deg) of Earth's penumbra at Moon distance.
        moon_radius_deg: Apparent lunar angular radius (deg).
    """

    def __init__(
        self,
        is_eclipse: bool,
        eclipse_type: Optional[str] = None,
        shadow_center_alt: Optional[float] = None,
        shadow_center_az: Optional[float] = None,
        umbra_radius_deg: Optional[float] = None,
        penumbra_radius_deg: Optional[float] = None,
        moon_radius_deg: Optional[float] = None,
    ):
        self.is_eclipse = is_eclipse
        self.eclipse_type = eclipse_type
        self.shadow_center_alt = shadow_center_alt
        self.shadow_center_az = shadow_center_az
        self.umbra_radius_deg = umbra_radius_deg
        self.penumbra_radius_deg = penumbra_radius_deg
        self.moon_radius_deg = moon_radius_deg


# Load planetary ephemerides (SPICE kernel)
eph = load("de442s.bsp")
earth = eph["earth"]
sun = eph["sun"]
moon = eph["moon"]


def calculate_lunar_eclipse_data(t: Time, observer: Topos) -> EclipseInfo:
    """Compute eclipse state and basic geometry at time `t` for `observer`.

    The method:
    1) Checks Sun-Moon separation (about 180 deg near eclipse).
    2) Computes umbra/penumbra radii from simple similar-triangles geometry.
    3) Classifies eclipse by Moon-shadow center separation.
    4) Derives the shadow-center alt/az as the anti-solar direction.

    Returns:
        EclipseInfo: Minimal set of eclipse parameters for display or logging.
    """
    # Apparent positions as seen from Earth's center (for separation, distances)
    sun_pos = earth.at(t).observe(sun)
    moon_pos = earth.at(t).observe(moon)

    # Quick cull if not near opposition
    separation_deg = sun_pos.separation_from(moon_pos).degrees
    if abs(separation_deg - 180.0) > 3.0:
        return EclipseInfo(is_eclipse=False, eclipse_type="none")

    # Radii (km)
    R_earth_km = 6371.0
    R_sun_km = 696340.0
    R_moon_km = 1737.4

    # Distances (km) from Earth's center
    D_sun_km = sun_pos.distance().km
    D_moon_km = moon_pos.distance().km

    # Angular radii at the Moon's distance (for Earth) and at Earth's distance (for Sun)
    earth_angular_radius_from_moon = math.asin(R_earth_km / D_moon_km)
    sun_angular_radius_from_earth = math.asin(R_sun_km / D_sun_km)
    moon_radius_rad = math.asin(R_moon_km / D_moon_km)

    # Umbra/Penumbra radii at the Moon
    umbra_radius_rad = max(0.0, earth_angular_radius_from_moon - sun_angular_radius_from_earth)
    penumbra_radius_rad = earth_angular_radius_from_moon + sun_angular_radius_from_earth

    # Angular distance from shadow center (about 180 deg - Sun-Moon separation)
    d_rad = math.radians(180.0 - separation_deg)

    # Classify state
    if d_rad > (penumbra_radius_rad + moon_radius_rad):
        eclipse_type = "none"
    elif d_rad > (umbra_radius_rad + moon_radius_rad):
        eclipse_type = "penumbral"
    elif d_rad > abs(umbra_radius_rad - moon_radius_rad):
        eclipse_type = "partial"
    else:
        eclipse_type = "total"

    # Shadow center from the local observer: opposite of the apparent Sun
    sun_astrometric = observer.at(t).observe(sun).apparent()
    s_alt, s_az, _ = sun_astrometric.altaz()
    shadow_center_alt = -s_alt.degrees
    shadow_center_az = (s_az.degrees + 180.0) % 360.0

    return EclipseInfo(
        is_eclipse=(eclipse_type != "none"),
        eclipse_type=eclipse_type,
        shadow_center_alt=shadow_center_alt,
        shadow_center_az=shadow_center_az,
        umbra_radius_deg=math.degrees(umbra_radius_rad),
        penumbra_radius_deg=math.degrees(penumbra_radius_rad),
        moon_radius_deg=math.degrees(moon_radius_rad),
    )


# Observer: Tokyo
tokyo = earth + Topos(latitude_degrees=35.6895, longitude_degrees=139.6917)
ts = load.timescale()

# Scan window in JST (converted to UTC internally)
start_jst = datetime(2025, 9, 8, 1, 0, tzinfo=timezone(timedelta(hours=9)))
end_jst = datetime(2025, 9, 8, 5, 0, tzinfo=timezone(timedelta(hours=9)))

dt = start_jst
while dt <= end_jst:
    t = ts.from_datetime(dt.astimezone(timezone.utc))
    eclipse = calculate_lunar_eclipse_data(t, tokyo)

    time_str = dt.strftime("%Y-%m-%d %H:%M")
    if eclipse.is_eclipse:
        print(f"{time_str} JST: [eclipse] {eclipse.eclipse_type.title()} Eclipse in progress")
        print(f"  - Umbra radius:    {eclipse.umbra_radius_deg:.2f} deg")
        print(f"  - Penumbra radius: {eclipse.penumbra_radius_deg:.2f} deg")
        print(f"  - Moon radius:     {eclipse.moon_radius_deg:.2f} deg")
        print(
            "  - Shadow center:   alt="
            f"{eclipse.shadow_center_alt:.1f} deg, az={eclipse.shadow_center_az:.1f} deg"
        )
    else:
        print(f"{time_str} JST: (No eclipse)")

    dt += timedelta(minutes=10)
