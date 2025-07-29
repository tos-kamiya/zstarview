import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import math
import os.path
import sys
import threading
import time

import pygame
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, GeocentricTrueEcliptic
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord
from skyfield.api import Topos
import skyfield.api

_dir = os.path.dirname(os.path.abspath(__file__))

EMOJI_FONT_PATH = os.path.join(_dir, "data", "Noto_Sans_Symbols", "NotoSansSymbols-VariableFont_wght.ttf")
CITY_COORD_FILE = os.path.join(_dir, "data", "cities1000.txt")
STARS_CSV_FILE = os.path.join(_dir, "data", "stars.csv")

TEXT_COLOR = (120, 120, 120)
TEXT_FONT_SIZE = 24

HORIZON_LINE_COLOR = (40, 40, 40)
CELESTIAL_EQUATOR_COLOR = (40, 40, 40)
ECLIPTIC_COLOR = (80, 60, 0)


def render_emoji(emoji, font_path: str, size: int = TEXT_FONT_SIZE) -> pygame.Surface:
    """
    Render emoji.

    Args:
        emoji (Any): Description.
        font_path (str): Description.
        size (int): Description.

    Returns:
        pygame.Surface: Description.
    """
    font = pygame.font.Font(font_path, size)
    surface = font.render(emoji, True, TEXT_COLOR)
    return surface


emoji_surfaces = {}


def load_city_coords(filename: str) -> dict[str, tuple[float, float]]:
    """
    Loads data: load city coords.

    Args:
        filename (str): Description.

    Returns:
        dict[str, tuple[float, float]]: Description.
    """
    city_table = {}
    with open(filename, encoding="utf-8") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 15:
                continue
            name = cols[1]
            lat = float(cols[4])
            lon = float(cols[5])
            country = cols[8]
            key = f"{country.lower()}/{name.lower()}"
            city_table[key] = (lat, lon)
    return city_table


def bv_to_color(bv: float) -> tuple[int, int, int]:
    """
    Bv to color.

    Args:
        bv (float): Description.

    Returns:
        tuple[int, int, int]: Description.
    """
    if bv < 0.0:
        return (170, 191, 255)
    elif bv < 0.3:
        return (202, 215, 255)
    elif bv < 0.6:
        return (248, 247, 255)
    elif bv < 1.0:
        return (255, 210, 161)
    else:
        return (255, 204, 111)


def mag_to_radius(vmag: float) -> float:
    """
    Mag to radius.

    Args:
        vmag (float): Description.

    Returns:
        float: Description.
    """
    base_radius = 7.0
    return max(0.1, base_radius * 10 ** (-0.2 * vmag))


def altaz_to_normalized_xy(alt: float, az: float) -> tuple[float, float]:
    """
    Altaz to normalized xy.

    Args:
        alt (float): Description.
        az (float): Description.

    Returns:
        tuple[float, float]: Description.
    """
    alt = float(alt)
    az = math.radians(float(az))
    r_norm = (90 - alt) / 90.0
    nx = -r_norm * math.sin(az)
    ny = -r_norm * math.cos(az)
    return (nx, ny)


def load_star_catalog(filename: str) -> list[dict[str, object]]:
    """
    Loads data: load star catalog.

    Args:
        filename (str): Description.

    Returns:
        list[dict[str, object]]: Description.
    """
    star_catalog = []
    with open(filename, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                ra = float(row["RA"]) * 15
                name = row["Name"]
                dec = float(row["Dec"])
                vmag = float(row["Vmag"])
                bv = float(row["B-V"])
                coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
                star_catalog.append({"name": name, "coord": coord, "vmag": vmag, "bv": bv})
            except Exception:
                continue
    return star_catalog


def update_star_positions(
    star_catalog: list[dict[str, object]], lat: float, lon: float, time_obj: Time
) -> tuple[list[dict[str, object]], Time, EarthLocation]:
    """
    Updates update star positions.

    Args:
        star_catalog (list[dict[str, object]]): Description.
        lat (float): Description.
        lon (float): Description.

    Returns:
        tuple[list[dict[str, object]], Time, EarthLocation]: Description.
    """
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    visible_stars = []
    for star in star_catalog:
        altaz = star["coord"].transform_to(AltAz(obstime=time_obj, location=location))
        if altaz.alt.deg > 0:
            visible_stars.append(
                {"alt": altaz.alt.deg, "az": altaz.az.deg, "vmag": star["vmag"], "bv": star["bv"], "name": star["name"]}
            )
    return (visible_stars, location)


def estimate_magnitude(body_name, observer, t: Time, planets: list[dict[str, object]]) -> float | None:
    """
    Estimate magnitude.

    Args:
        body_name (Any): Description.
        observer (Any): Description.
        t (Time): Description.
        planets (list[dict[str, object]]): Description.

    Returns:
        float | None: Description.
    """
    earth = planets["earth"]
    sun = planets["sun"]
    planet = planets[body_name]
    obs_to_planet = observer.at(t).observe(planet).apparent()
    obs_to_sun = observer.at(t).observe(sun).apparent()
    phase_angle = obs_to_planet.separation_from(planet.at(t).observe(sun).apparent())
    r = planet.at(t).observe(sun).apparent().distance().au
    delta = obs_to_planet.distance().au
    i = phase_angle.degrees
    if body_name == "mercury":
        mag = -0.42 + 0.038 * i - 0.000273 * i**2 + 2e-06 * i**3
    elif body_name == "venus":
        mag = -4.4 + 0.0009 * i + 0.000239 * i**2 - 6.5e-07 * i**3
    elif body_name == "mars":
        mag = -1.52 + 0.016 * i
    elif body_name == "jupiter barycenter":
        mag = -9.4 + 0.005 * i
    elif body_name == "saturn barycenter":
        B = math.radians(26.73)
        mag = -8.88 + 0.044 * i - 2.6 * math.sin(B) + 1.25 * math.sin(B) ** 2
    else:
        return None
    mag += 5 * math.log10(r * delta)
    return mag


PLANET_SYMBOLS = {
    "sun": "☀",
    "moon": "🌛",
    "mercury": "☿",
    "venus": "♀",
    "mars": "♂",
    "jupiter barycenter": "♃",
    "saturn barycenter": "♄",
}


def get_visible_planets(lat: float, lon: float, astropy_time: Time) -> list[dict[str, object]]:
    """
    Gets visible planets using a given astropy Time object.

    Args:
        lat (float): Latitude in degrees.
        lon (float): Longitude in degrees.
        astropy_time (Time): Astropy Time object.

    Returns:
        list[dict[str, object]]: List of visible planets with their positions.
    """
    ts = skyfield.api.load.timescale()
    t = ts.from_astropy(astropy_time)  # ← astropy -> skyfield の変換
    planets = skyfield.api.load("de421.bsp")
    observer = planets["earth"] + Topos(latitude_degrees=lat, longitude_degrees=lon)

    visible_bodies = []
    for name, symbol in PLANET_SYMBOLS.items():
        planet = planets[name]
        astrometric = observer.at(t).observe(planet).apparent()
        alt, az, _ = astrometric.altaz()
        if alt.degrees > 0:
            mag = None
            if name != "sun":
                mag = estimate_magnitude(name, observer, t, planets)
            body = {"alt": alt.degrees, "az": az.degrees, "symbol": symbol, "mag": mag, "name": name}
            if name == "moon":
                body["phase_angle"] = moon_phase_angle(observer, t, planets)
            visible_bodies.append(body)
    return visible_bodies



def moon_phase_angle(observer, t: Time, planets: list[dict[str, object]]) -> float:
    """
    Moon phase angle.

    Args:
        observer (Any): Description.
        t (Time): Description.
        planets (list[dict[str, object]]): Description.

    Returns:
        float: Description.
    """
    moon = planets["moon"]
    sun = planets["sun"]
    e_to_moon = observer.at(t).observe(moon).apparent()
    e_to_sun = observer.at(t).observe(sun).apparent()
    phase_angle = e_to_moon.separation_from(e_to_sun)
    return phase_angle.degrees


def draw_moon(screen: pygame.Surface, x: int, y: int, radius: float, phase_angle_deg: float):
    """
    Draws graphical elements: draw moon.

    Args:
        screen (pygame.Surface): Description.
        x (int): Description.
        y (int): Description.
        radius (float): Description.
        phase_angle_deg (float): Description.
    """
    phase = phase_angle_deg % 360
    brightness = (1 + math.cos(math.radians(phase))) / 2
    pygame.draw.circle(screen, (230, 230, 230), (int(x), int(y)), int(radius))
    mask = pygame.Surface((radius * 2, radius * 2), pygame.SRCALPHA)
    mask.fill((0, 0, 0, 0))
    if phase < 180:
        pygame.draw.ellipse(mask, (0, 0, 0, 255), (0, 0, radius * 2, radius * 2))
        pygame.draw.ellipse(mask, (0, 0, 0, 0), (radius * (1 - brightness), 0, radius * 2, radius * 2))
    else:
        pygame.draw.ellipse(mask, (0, 0, 0, 255), (0, 0, radius * 2, radius * 2))
        pygame.draw.ellipse(mask, (0, 0, 0, 0), (radius * (brightness - 1), 0, radius * 2, radius * 2))
    screen.blit(mask, (int(x - radius), int(y - radius)), special_flags=pygame.BLEND_RGBA_SUB)


def draw_cross_gauge(screen: pygame.Surface, color: tuple[int, int, int], x: int, y: int):
    """
    Draws graphical elements: draw cross gauge.

    Args:
        screen (pygame.Surface): Description.
        color (tuple[int, int, int]): Description.
        x (int): Description.
        y (int): Description.
    """
    cross_outer_len = 15
    cross_inner_len = 4
    pygame.draw.line(screen, color, (x - cross_outer_len, y), (x - cross_inner_len, y), 1)
    pygame.draw.line(screen, color, (x + cross_inner_len, y), (x + cross_outer_len, y), 1)
    pygame.draw.line(screen, color, (x, y - cross_outer_len), (x, y - cross_inner_len), 1)
    pygame.draw.line(screen, color, (x, y + cross_inner_len), (x, y + cross_outer_len), 1)


def rotate_points_at_max_azimuth_gap(
    points: list[tuple[float, float]], azimuths: list[float]
) -> list[tuple[float, float]]:
    """
    Rotates points to avoid a large azimuth jump in the middle of the line.

    Args:
        points: List of (x, y) tuples.
        azimuths: Corresponding azimuth values in degrees.

    Returns:
        Rotated list of (x, y) points.
    """
    if len(points) < 2:
        return points

    max_gap = -1
    max_index = 0
    for i in range(1, len(azimuths)):
        delta = abs(azimuths[i] - azimuths[i - 1])
        delta = min(delta, 360 - delta)  # wrap-around考慮
        if delta > max_gap:
            max_gap = delta
            max_index = i

    rotated_points = points[max_index:] + points[:max_index]
    return rotated_points


def calculate_celestial_equator_points_norm(location: EarthLocation, time: Time) -> list[tuple[float, float]]:
    """
    Calculates calculate celestial equator points norm.

    Args:
        location (EarthLocation): Description.
        time (Time): Description.

    Returns:
        list[tuple[float, float]]: Description.
    """
    points = []
    azimuths = []
    for ra_deg in range(0, 360, 5):
        coord = SkyCoord(ra=ra_deg * u.deg, dec=0 * u.deg, frame="icrs")
        altaz = coord.transform_to(AltAz(obstime=time, location=location))
        if altaz.alt.deg > 0:
            nx, ny = altaz_to_normalized_xy(altaz.alt.deg, altaz.az.deg)
            points.append((nx, ny))
            azimuths.append(altaz.az.deg)
    return rotate_points_at_max_azimuth_gap(points, azimuths)


def calculate_ecliptic_points_norm(location: EarthLocation, time: Time) -> list[tuple[float, float]]:
    """
    Calculates calculate ecliptic points norm.

    Args:
        location (EarthLocation): Description.
        time (Time): Description.

    Returns:
        list[tuple[float, float]]: Description.
    """
    points = []
    azimuths = []
    for lon_deg in range(0, 360, 5):
        ecl = SkyCoord(lon=lon_deg * u.deg, lat=0 * u.deg, frame=GeocentricTrueEcliptic(obstime=time))
        icrs = ecl.transform_to("icrs")
        altaz = icrs.transform_to(AltAz(obstime=time, location=location))
        if altaz.alt.deg > 0:
            nx, ny = altaz_to_normalized_xy(altaz.alt.deg, altaz.az.deg)
            points.append((nx, ny))
            azimuths.append(altaz.az.deg)
    return rotate_points_at_max_azimuth_gap(points, azimuths)


@dataclass
class SkyData:
    location: EarthLocation
    time: Time
    planets: list[dict[str, object]]
    stars: list[dict[str, object]]
    celestial_equator_points: list[tuple[float, float]]
    ecliptic_points: list[tuple[float, float]]


def draw_sky(
    screen: pygame.Surface,
    sky_data: SkyData,
):
    """
    Draws graphical elements: draw sky.

    Args:
        screen (pygame.Surface): Description.
        sky_data (SkyData): Data for stars, planets, etc.
    """
    stars = sky_data.stars
    planets = sky_data.planets

    font = pygame.font.SysFont(None, TEXT_FONT_SIZE)
    screen.fill((0, 0, 0))
    center, radius = get_screen_geometry(screen)
    pygame.draw.circle(screen, HORIZON_LINE_COLOR, center, radius, 1)

    def to_screen_xy(nx, ny):
        return (center[0] + nx * radius, center[1] + ny * radius)

    time_text = font.render(sky_data.time.strftime("%Y-%m-%d %H:%M:%S UTC"), True, (180, 180, 180))
    screen.blit(time_text, (10, 10))
    directions = {"N": 0, "NE": 45, "E": 90, "SE": 135, "S": 180, "SW": 225, "W": 270, "NW": 315}
    for label, angle in directions.items():
        nx, ny = altaz_to_normalized_xy(0, angle)
        x, y = to_screen_xy(nx, ny)
        text = font.render(label, True, (150, 150, 150))
        text_rect = text.get_rect(center=(x, y))
        screen.blit(text, text_rect)
    if len(sky_data.celestial_equator_points) >= 2:
        points = [to_screen_xy(nx, ny) for nx, ny in sky_data.celestial_equator_points]
        pygame.draw.lines(screen, CELESTIAL_EQUATOR_COLOR, False, points, 1)
    if len(sky_data.ecliptic_points) >= 2:
        points = [to_screen_xy(nx, ny) for nx, ny in sky_data.ecliptic_points]
        pygame.draw.lines(screen, ECLIPTIC_COLOR, False, points, 1)

    mouse_x, mouse_y = pygame.mouse.get_pos()
    highlight = None
    min_dist = 9999
    planet_screen_coords = {
        body["name"]: to_screen_xy(*altaz_to_normalized_xy(body["alt"], body["az"])) for body in planets
    }
    star_screen_coords = {
        star["name"]: to_screen_xy(*altaz_to_normalized_xy(star["alt"], star["az"])) for star in stars
    }

    def check_proximity(x, y, obj):
        nonlocal highlight, min_dist
        dist = math.hypot(mouse_x - x, mouse_y - y)
        if dist < 30 and dist < min_dist:
            highlight = (obj, x, y)
            min_dist = dist

    for body in planets:
        x, y = planet_screen_coords[body["name"]]
        check_proximity(x, y, body)
    for star in stars:
        x, y = star_screen_coords[star["name"]]
        check_proximity(x, y, star)
    for body in planets:
        x, y = planet_screen_coords[body["name"]]
        name = body.get("name") or ""
        if name == "sun":
            draw_cross_gauge(screen, TEXT_COLOR, x, y)
        elif name == "moon":
            moon_radius = 0.5 / 2 * (radius / 90)
            draw_moon(screen, x, y, moon_radius, body["phase_angle"])
            draw_cross_gauge(screen, TEXT_COLOR, x, y)
        else:
            emoji_surface = emoji_surfaces.get(name)
            if emoji_surface:
                rect = emoji_surface.get_rect(center=(x, y))
                screen.blit(emoji_surface, rect)
    for star in stars:
        x, y = star_screen_coords[star["name"]]
        color = bv_to_color(star["bv"])
        rad = mag_to_radius(star["vmag"]) * radius / 2000
        size = int(math.ceil(rad * 2))
        if size < 3:
            alpha = int(255 * rad)
            dot = pygame.Surface((3, 3), pygame.SRCALPHA)
            dot.fill((*color, alpha))
            screen.blit(dot, (int(x) - 1, int(y) - 1))
        else:
            surface = pygame.Surface((size, size), pygame.SRCALPHA)
            center_px = size // 2
            for i in range(size):
                for j in range(size):
                    dx = i - center_px
                    dy = j - center_px
                    dist_sq = dx * dx + dy * dy
                    if dist_sq <= rad * rad:
                        alpha = int(255 * math.exp(-dist_sq / (0.5 * rad * rad)))
                        surface.set_at((i, j), (*color, alpha))
            screen.blit(surface, (int(x - center_px), int(y - center_px)), special_flags=pygame.BLEND_RGBA_ADD)
    if highlight:
        obj, x, y = highlight
        pygame.draw.circle(screen, TEXT_COLOR, (int(x), int(y)), 10, 2)
        name = obj.get("name") or ""
        label = font.render(name, True, TEXT_COLOR)
        screen.blit(label, (x + 15, y - 15))
    pygame.display.flip()


def get_screen_geometry(screen: pygame.Surface) -> tuple[tuple[int, int], int]:
    """
    Gets get screen geometry.

    Args:
        screen (pygame.Surface): Description.

    Returns:
        tuple[tuple[int, int], int]: Description.
    """
    width, height = screen.get_size()
    center = (width // 2, height // 2)
    s1, s2 = (width, height)
    if s1 < s2:
        s1, s2 = (s2, s1)
    radius = max(50, int(s1 * 0.7 + s2 * (1.0 - 0.7)) // 2 - 10)
    return (center, radius)


def show_splash_until_thread_finishes(screen, font, message: str, thread: threading.Thread):
    """
    Show splash screen with message while a background thread is running.
    Respond to OS events to avoid hang detection.

    Args:
        screen: Pygame screen surface.
        font: Pygame font object.
        message: Message to display.
        thread: Thread to monitor for completion.
    """
    dot_count = 0
    last_update = 0

    while thread.is_alive():
        now = time.time()
        if now - last_update > 0.5:
            dot_count = (dot_count + 1) % 4
            dots = "." * dot_count
            text = message + dots
            screen.fill((0, 0, 0))
            label = font.render(text, True, (255, 255, 255))
            rect = label.get_rect(center=(screen.get_width() // 2, screen.get_height() // 2))
            screen.blit(label, rect)
            pygame.display.flip()
            last_update = now

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
        time.sleep(0.01)


def main():
    """
    Main entry point for the star sky visualizer.
    """
    city_table = load_city_coords(CITY_COORD_FILE)
    city = "Tokyo"
    if len(sys.argv) >= 2:
        city = sys.argv[1]
    city = city.lower()
    if city not in city_table:
        candidate_cities = []
        for c in city_table.keys():
            if c.endswith("/" + city):
                candidate_cities.append(c)
        if len(candidate_cities) == 0:
            print(f"Unknown city: {city}")
            return
        elif len(candidate_cities) >= 2:
            print(f"Specify explicit country name: {candidate_cities}")
            return
        else:
            city = candidate_cities[0]

    lat, lon = city_table[city]
    star_catalog = load_star_catalog(STARS_CSV_FILE)

    sky_data = None

    def background_data_loader():
        nonlocal sky_data
        try:
            now = datetime.now(timezone.utc)
            time_obj = Time(now)
            stars, loc = update_star_positions(star_catalog, lat, lon, time_obj)
            planets = get_visible_planets(lat, lon, time_obj)
            celestial_equator_points = calculate_celestial_equator_points_norm(loc, time_obj)
            ecliptic_points = calculate_ecliptic_points_norm(loc, time_obj)
            sky_data = SkyData(loc, time_obj, planets, stars, celestial_equator_points, ecliptic_points)
        except Exception as e:
            import traceback

            print("Error in background_data_loader():", file=sys.stderr)
            traceback.print_exc()

    loading_thread = threading.Thread(target=background_data_loader)
    loading_thread.start()

    pygame.init()
    for p, s in PLANET_SYMBOLS.items():
        emoji_surfaces[p] = render_emoji(s, EMOJI_FONT_PATH)

    splash_size = (600, 300)
    splash = pygame.display.set_mode(splash_size)
    pygame.display.set_caption(f"Zenith Star View - {city}")
    font = pygame.font.SysFont(None, 32)
    show_splash_until_thread_finishes(splash, font, "Loading star data...", loading_thread)

    assert sky_data is not None

    screen_size = (800, 800)
    screen = pygame.display.set_mode(screen_size, pygame.RESIZABLE)
    pygame.display.set_caption(f"Zenith Star View - {city}")

    sky_data_lock = threading.Lock()

    with sky_data_lock:
        draw_sky(screen, sky_data)

    new_sky_data = None

    def background_updater():
        nonlocal new_sky_data
        while running:
            time.sleep(300)
            now = datetime.now(timezone.utc)
            time_obj = Time(now)
            stars, loc = update_star_positions(star_catalog, lat, lon, time_obj)
            planets = get_visible_planets(lat, lon, time_obj)
            celestial_equator_points = calculate_celestial_equator_points_norm(loc, time_obj)
            ecliptic_points = calculate_ecliptic_points_norm(loc, time_obj)
            with sky_data_lock:
                new_sky_data = SkyData(loc, time_obj, planets, stars, celestial_equator_points, ecliptic_points)

    fullscreen = False
    running = True
    thread = threading.Thread(target=background_updater, daemon=True)
    thread.start()
    last_mouse_pos = (-1, -1)
    while running:
        redraw_needed = False
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.VIDEORESIZE:
                if not fullscreen:
                    screen_size = event.size
                redraw_needed = True
            elif event.type == pygame.MOUSEMOTION:
                if pygame.mouse.get_pos() != last_mouse_pos:
                    last_mouse_pos = pygame.mouse.get_pos()
                    redraw_needed = True
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_F11:
                    fullscreen = not fullscreen
                    if fullscreen:
                        screen = pygame.display.set_mode(screen_size, pygame.RESIZABLE | pygame.NOFRAME)
                    else:
                        screen = pygame.display.set_mode(screen_size, pygame.RESIZABLE)
                    redraw_needed = True
                elif event.key == pygame.K_ESCAPE:
                    if fullscreen:
                        fullscreen = False
                        screen = pygame.display.set_mode(screen_size, pygame.RESIZABLE)
                        redraw_needed = True
        with sky_data_lock:
            if new_sky_data is not None:
                sky_data, new_sky_data = new_sky_data, None
                redraw_needed = True
            if redraw_needed:
                draw_sky(screen, sky_data)
        time.sleep(0.05)
    pygame.quit()


if __name__ == "__main__":
    main()
