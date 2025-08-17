import math
from datetime import datetime, timedelta, timezone

from skyfield.api import load, Topos
from skyfield.timelib import Time


# EclipseInfo構造体（簡略版）
class EclipseInfo:
    def __init__(
        self,
        is_eclipse: bool,
        eclipse_type: str = None,
        shadow_center_alt=None,
        shadow_center_az=None,
        umbra_radius_deg=None,
        penumbra_radius_deg=None,
        moon_radius_deg=None,
    ):
        self.is_eclipse = is_eclipse
        self.eclipse_type = eclipse_type
        self.shadow_center_alt = shadow_center_alt
        self.shadow_center_az = shadow_center_az
        self.umbra_radius_deg = umbra_radius_deg
        self.penumbra_radius_deg = penumbra_radius_deg
        self.moon_radius_deg = moon_radius_deg


# カーネル読み込み
eph = load("de440s.bsp")
earth = eph["earth"]
sun = eph["sun"]
moon = eph["moon"]


def calculate_lunar_eclipse_data(t: Time, observer: Topos) -> EclipseInfo:
    sun_pos = earth.at(t).observe(sun)
    moon_pos = earth.at(t).observe(moon)

    separation_deg = sun_pos.separation_from(moon_pos).degrees
    if abs(separation_deg - 180.0) > 3.0:
        return EclipseInfo(is_eclipse=False)

    R_earth_km = 6371.0
    R_sun_km = 696340.0
    R_moon_km = 1737.4

    D_sun_km = sun_pos.distance().km
    D_moon_km = moon_pos.distance().km

    earth_angular_radius_from_moon = math.asin(R_earth_km / D_moon_km)
    sun_angular_radius_from_earth = math.asin(R_sun_km / D_sun_km)

    umbra_radius_rad = max(0.0, earth_angular_radius_from_moon - sun_angular_radius_from_earth)
    penumbra_radius_rad = earth_angular_radius_from_moon + sun_angular_radius_from_earth
    d_rad = math.radians(180.0 - separation_deg)
    moon_radius_rad = math.asin(R_moon_km / D_moon_km)

    # 食の状態を判定
    if d_rad > (penumbra_radius_rad + moon_radius_rad):
        eclipse_type = "none"
    elif d_rad > (umbra_radius_rad + moon_radius_rad):
        eclipse_type = "penumbral"
    elif d_rad > abs(umbra_radius_rad - moon_radius_rad):
        eclipse_type = "partial"
    else:
        eclipse_type = "total"

    sun_astrometric = observer.at(t).observe(sun).apparent()
    s_alt, s_az, _ = sun_astrometric.altaz()
    shadow_center_alt = -s_alt.degrees
    shadow_center_az = (s_az.degrees + 180) % 360

    return EclipseInfo(
        is_eclipse=True,
        eclipse_type=eclipse_type,
        shadow_center_alt=shadow_center_alt,
        shadow_center_az=shadow_center_az,
        umbra_radius_deg=math.degrees(umbra_radius_rad),
        penumbra_radius_deg=math.degrees(penumbra_radius_rad),
        moon_radius_deg=math.degrees(moon_radius_rad),
    )


# 観測地：東京
tokyo = earth + Topos(latitude_degrees=35.6895, longitude_degrees=139.6917)
ts = load.timescale()

# JST → UTC に変換
start_jst = datetime(2025, 9, 8, 1, 0, tzinfo=timezone(timedelta(hours=9)))
end_jst = datetime(2025, 9, 8, 5, 0, tzinfo=timezone(timedelta(hours=9)))

dt = start_jst
while dt <= end_jst:
    t = ts.from_datetime(dt.astimezone(timezone.utc))
    eclipse = calculate_lunar_eclipse_data(t, tokyo)

    time_str = dt.strftime("%Y-%m-%d %H:%M")
    if eclipse.is_eclipse:
        print(f"{time_str} JST: 🌒 {eclipse.eclipse_type.title()} Eclipse in progress")
        print(f"  - Umbra radius:    {eclipse.umbra_radius_deg:.2f}°")
        print(f"  - Penumbra radius: {eclipse.penumbra_radius_deg:.2f}°")
        print(f"  - Moon radius:     {eclipse.moon_radius_deg:.2f}°")
        print(f"  - Shadow center:   alt={eclipse.shadow_center_alt:.1f}°, az={eclipse.shadow_center_az:.1f}°")
    else:
        print(f"{time_str} JST: (No eclipse)")

    dt += timedelta(minutes=10)
