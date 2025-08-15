import math
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from PIL import Image
from zoneinfo import ZoneInfo

from PySide6.QtCore import QPoint, QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF, QRadialGradient

from ..paths import (
    CELESTIAL_EQUATOR_COLOR,
    DIRECTIONS,
    ECLIPTIC_COLOR,
    FIELD_OF_VIEW_DEG,
    HORIZON_LINE_COLOR,
    TEXT_COLOR,
)
from ..types import PlanetBody, ScreenGeometry, SkyData, ViewerData
from ..astro import altaz_to_normalized_xy, calculate_sun_angle_on_moon, is_in_fov
from ..utils.image import generate_moon_phase_image
from ..utils.qt import pil2qpixmap


def bv_to_rgb_vectorized(bv: np.ndarray) -> np.ndarray:
    """Vectorized conversion of B-V color index to RGB tuples."""
    # First, initialize all stars to the default color (Orange-ish)
    # The output will be an array of shape (number_of_stars, 3)
    rgb = np.full((bv.shape[0], 3), [255, 204, 111], dtype=int)

    # Overwrite rows with the correct color based on conditions
    rgb[bv < 0.0] = [170, 191, 255]  # Blueish
    rgb[(bv >= 0.0) & (bv < 0.3)] = [202, 215, 255]  # White-Blue
    rgb[(bv >= 0.3) & (bv < 0.6)] = [248, 247, 255]  # White
    rgb[(bv >= 0.6) & (bv < 1.0)] = [255, 210, 161]  # Yellowish
    return rgb


def altaz_to_normalized_xy_vectorized(
    alt: np.ndarray, az: np.ndarray, view_center: Tuple[float, float]
) -> Tuple[np.ndarray, np.ndarray]:
    """Vectorized conversion of alt/az to normalized screen coordinates."""
    center_alt, center_az = view_center
    alt1, az1 = np.radians(center_alt), np.radians(center_az)
    alt2, az2 = np.radians(alt), np.radians(az)

    cos_theta = np.sin(alt1) * np.sin(alt2) + np.cos(alt1) * np.cos(alt2) * np.cos(az2 - az1)
    theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))

    r = theta / (math.pi / 2)

    dx = np.cos(alt2) * np.sin(az2 - az1)
    dy = np.cos(alt1) * np.sin(alt2) - np.sin(alt1) * np.cos(alt2) * np.cos(az2 - az1)
    length = np.hypot(dx, dy)
    # Avoid division by zero for objects at the pole
    length[length == 0] = 1.0
    nx = r * dx / length
    ny = -r * dy / length
    return (nx, ny)


def normalized_to_screen_xy_vectorized(nx: np.ndarray, ny: np.ndarray, geometry: ScreenGeometry) -> Tuple[np.ndarray, np.ndarray]:
    """Vectorized conversion of normalized coordinates to screen coordinates."""
    return (geometry.center[0] + nx * geometry.radius, geometry.center[1] + ny * geometry.radius)


def find_highlighted_object(
    sky_data: Optional[SkyData],
    viewer_data: ViewerData,
    mouse_pos: QPoint,
    geometry: ScreenGeometry,
) -> Optional[Tuple[Dict[str, Union[str, float]], QPointF]]:
    """Find the nearest celestial object to the mouse cursor (vectorized for stars)."""
    min_dist_sq = 30**2  # squared pixels
    highlighted_object = None

    if not sky_data:
        return None

    # Handle stars first (vectorized)
    stars = sky_data.stars
    if stars["alt"].size > 0:
        nx, ny = altaz_to_normalized_xy_vectorized(stars["alt"], stars["az"], viewer_data.view_center)
        x, y = normalized_to_screen_xy_vectorized(nx, ny, geometry)
        dist_sq = (mouse_pos.x() - x) ** 2 + (mouse_pos.y() - y) ** 2
        closest_star_idx = np.argmin(dist_sq)
        if dist_sq[closest_star_idx] < min_dist_sq:
            min_dist_sq = dist_sq[closest_star_idx]
            # Reconstruct a dictionary for the single highlighted star
            highlighted_star = {key: val[closest_star_idx] for key, val in stars.items()}
            highlighted_object = (highlighted_star, QPointF(x[closest_star_idx], y[closest_star_idx]))

    # Handle planets (scalar)
    for body in sky_data.planets:
        nx, ny = altaz_to_normalized_xy(body.alt, body.az, viewer_data.view_center)
        pos = normalized_to_screen_xy(nx, ny, geometry)
        dist_sq = (mouse_pos.x() - pos.x()) ** 2 + (mouse_pos.y() - pos.y()) ** 2
        if dist_sq < min_dist_sq:
            min_dist_sq = dist_sq
            highlighted_object = (body, pos)  # body is already an object/dataclass

    return highlighted_object


def draw_radial_background(painter: QPainter, rect: QRectF, geometry: ScreenGeometry):
    assert geometry.radius >= 10
    fov_middle = 90 + (FIELD_OF_VIEW_DEG / 2 - 90) / 2
    r90 = float(geometry.radius)
    r_fov = float(geometry.radius * (fov_middle / 90))
    r_max = float(r_fov * 1.4)
    step_px = 0.5

    def pos(r):
        return max(0.0, min(1.0, r / r_max))

    def col(r, s):
        return QColor(0, 0, 0, max(0, 255 - (s + int(150 * (r - r90) / r_max))))

    c = geometry.center
    g = QRadialGradient(QPoint(c[0], c[1]), r_max)
    g.setColorAt(pos(0), col(r90, 0))
    g.setColorAt(pos(r90), col(r90, 0))
    g.setColorAt(pos(r90 + step_px), col(r90, 10))
    g.setColorAt(pos(r_fov), col(r_fov, 10))
    g.setColorAt(pos(r_fov + step_px), col(r_fov, 20))
    g.setColorAt(1.0, col(r_max, 20))

    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(g)
    painter.drawRect(rect)
    painter.restore()


def split_by_gaps(points: List[Tuple[float, float]]) -> List[List[Tuple[float, float]]]:
    """Split a polyline by gaps to avoid drawing large jumps."""

    def dist(p1: Tuple[float, float], p2: Tuple[float, float]) -> float:
        return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

    fragments: List[List[Tuple[float, float]]] = [[]]
    for p in points:
        if not fragments[-1] or dist(p, fragments[-1][-1]) < 0.3:
            fragments[-1].append(p)
        else:
            fragments.append([p])
    return fragments


def draw_sky_reference_lines(painter: QPainter, geometry: ScreenGeometry, sky_data: SkyData):
    """Draw equator, ecliptic, and horizon polylines."""
    point_list_pen_styles = [
        (sky_data.celestial_equator_points, (CELESTIAL_EQUATOR_COLOR, 2, Qt.PenStyle.DashLine)),
        (sky_data.ecliptic_points, (ECLIPTIC_COLOR, 2, Qt.PenStyle.DotLine)),
        (sky_data.horizon_points, (HORIZON_LINE_COLOR, 2, Qt.PenStyle.SolidLine)),
    ]
    for points, pen_style in point_list_pen_styles:
        for frag in split_by_gaps(points):
            if len(frag) >= 2:
                pts = [normalized_to_screen_xy(nx, ny, geometry) for nx, ny in frag]
                poly = QPolygonF(pts)
                painter.setPen(QPen(*pen_style))
                painter.drawPolyline(poly)


def draw_stars(painter: QPainter, geometry: ScreenGeometry, sky_data: SkyData, viewer_data: ViewerData, star_base_radius: float):
    """Draw stars using vectorized calculations and batched drawing for small stars."""
    stars = sky_data.stars
    if stars["alt"].size == 0:
        return

    # 1. Vectorized calculations for all stars
    nx, ny = altaz_to_normalized_xy_vectorized(stars["alt"], stars["az"], viewer_data.view_center)
    x, y = normalized_to_screen_xy_vectorized(nx, ny, geometry)
    sizes = np.maximum(0.1, star_base_radius * 10 ** (-0.2 * stars["vmag"])) * geometry.radius / 500
    bv_clamped = np.nan_to_num(stars["bv"], nan=0.45)
    rgb_colors = bv_to_rgb_vectorized(bv_clamped)

    # 2. Split stars into small and large
    large_star_threshold = 4.0
    small_star_mask = sizes < large_star_threshold
    large_star_mask = ~small_star_mask

    painter.setPen(Qt.PenStyle.NoPen)

    # 3. Draw large stars individually for high quality
    if np.any(large_star_mask):
        large_x, large_y, large_sizes, large_rgb = (
            x[large_star_mask],
            y[large_star_mask],
            sizes[large_star_mask],
            rgb_colors[large_star_mask],
        )
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
        for i in range(len(large_x)):
            pos = QPointF(large_x[i], large_y[i])
            color = QColor(int(large_rgb[i][0]), int(large_rgb[i][1]), int(large_rgb[i][2]))
            radius = math.sqrt(large_sizes[i])

            gradient = QRadialGradient(pos, radius)
            gradient.setColorAt(0, color)
            color_transparent = QColor(color)
            color_transparent.setAlpha(0)
            gradient.setColorAt(1, color_transparent)
            painter.setBrush(gradient)
            painter.drawEllipse(pos, radius, radius)

    # 4. Draw small stars in batches
    if np.any(small_star_mask):
        small_x, small_y, small_sizes, small_rgb = (
            x[small_star_mask],
            y[small_star_mask],
            sizes[small_star_mask],
            rgb_colors[small_star_mask],
        )
        unique_colors = np.unique(small_rgb, axis=0)
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
        for r, g, b in unique_colors:
            color = QColor(int(r), int(g), int(b))
            painter.setBrush(color)
            mask = np.all(small_rgb == (r, g, b), axis=1)
            rects = [
                QRectF(cx - s / 2, cy - s / 2, s, s)
                for cx, cy, s in zip(small_x[mask], small_y[mask], small_sizes[mask])
            ]
            painter.drawRects(rects)

    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)


def draw_gauge_cross(painter: QPainter, color: QColor, center: QPointF):
    cross_outer_len, cross_inner_len = 15, 4
    x, y = center.x(), center.y()
    painter.setPen(QPen(color, 1))
    painter.drawLine(QPointF(x - cross_outer_len, y), QPointF(x - cross_inner_len, y))
    painter.drawLine(QPointF(x + cross_inner_len, y), QPointF(x + cross_outer_len, y))
    painter.drawLine(QPointF(x, y - cross_outer_len), QPointF(x, y - cross_inner_len))
    painter.drawLine(QPointF(x, y + cross_inner_len), QPointF(x, y + cross_outer_len))


def draw_zenith_marker(painter: QPainter, geometry: ScreenGeometry, view_center: Tuple[float, float]):
    alt_zenith = 90.0
    az_ref = view_center[1]

    if not is_in_fov(alt_zenith, az_ref, view_center):
        return

    nx, ny = altaz_to_normalized_xy(alt_zenith, az_ref, view_center)
    pos = normalized_to_screen_xy(nx, ny, geometry)

    s = 7
    x, y = pos.x(), pos.y()

    painter.setPen(QPen(TEXT_COLOR, 1))
    painter.drawLine(QPointF(x - s, y - s), QPointF(x + s, y + s))
    painter.drawLine(QPointF(x - s, y + s), QPointF(x + s, y - s))


def draw_moon(
    painter: QPainter,
    center: QPointF,
    radius: float,
    phase_angle_deg: float,
    sun_altaz: Optional[Tuple[float, float]] = None,
    moon_altaz: Optional[Tuple[float, float]] = None,
    opacity: float = 1.0,
):
    """Draw the moon with phase and optional rotation to align with sun angle."""
    img_size = int(radius * 2)
    if img_size < 5:
        img_size = 5

    view_dir = np.array([0, 0, 1])
    phase_angle_rad = math.radians(phase_angle_deg)
    sun_dir = np.array([np.sin(phase_angle_rad), 0, -np.cos(phase_angle_rad)])
    sun_dir /= np.linalg.norm(sun_dir)
    moon_img_pil = generate_moon_phase_image(img_size, sun_dir, view_dir)

    rotate_deg = 0
    if sun_altaz is not None and moon_altaz is not None:
        angle = calculate_sun_angle_on_moon(moon_altaz, sun_altaz)
        rotate_deg = -math.degrees(angle) - 90

    moon_img_pil = moon_img_pil.rotate(rotate_deg, resample=Image.Resampling.BICUBIC, expand=False)

    pixmap = pil2qpixmap(moon_img_pil)
    target_rect = QRectF(center.x() - img_size / 2, center.y() - img_size / 2, img_size, img_size)
    painter.save()
    painter.setOpacity(opacity)
    painter.drawPixmap(target_rect, pixmap, QRectF(0, 0, img_size, img_size))
    painter.restore()


def draw_planets(
    painter: QPainter,
    geometry: ScreenGeometry,
    sky_data: SkyData,
    viewer_data: ViewerData,
    enlarge_moon: bool,
    emoji_font: QFont,
):
    sun_altaz: Optional[Tuple[float, float]] = None
    moon_altaz: Optional[Tuple[float, float]] = None
    for body in sky_data.planets:
        if body.name == "sun":
            sun_altaz = (body.alt, body.az)
        if body.name == "moon":
            moon_altaz = (body.alt, body.az)

    for body in sky_data.planets:
        pos = normalized_to_screen_xy(*altaz_to_normalized_xy(body.alt, body.az, viewer_data.view_center), geometry)
        if body.name == "sun":
            draw_gauge_cross(painter, TEXT_COLOR, pos)
        elif body.name == "moon":
            moon_radius = 0.5 * (1 if not enlarge_moon else 3) / 2 * (geometry.radius / 90.0)
            draw_moon(
                painter,
                pos,
                moon_radius,
                body.phase_angle if body.phase_angle is not None else 0.0,
                sun_altaz=sun_altaz,
                moon_altaz=moon_altaz,
                opacity=1.0 if not enlarge_moon else 0.7,
            )
            draw_gauge_cross(painter, TEXT_COLOR, pos)
        else:
            painter.setFont(emoji_font)
            painter.setPen(TEXT_COLOR)
            painter.drawText(pos, body.symbol)


def draw_direction_labels(painter: QPainter, geometry: ScreenGeometry, view_center: Tuple[float, float], text_font: QFont):
    painter.setPen(TEXT_COLOR)
    painter.setFont(text_font)
    alt = 0
    for label, az in DIRECTIONS.items():
        if not is_in_fov(alt, az, view_center):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        pos = normalized_to_screen_xy(nx, ny, geometry)
        painter.drawText(pos, label)


def draw_overlay_info(
    painter: QPainter,
    sky_data: SkyData,
    viewer_data: ViewerData,
    highlighted_object: Optional[Tuple[Dict[str, Union[str, float]], QPointF]],
    text_font: QFont,
):
    # ---- Local time ----
    utc_time = sky_data.time
    tz_name = viewer_data.timezone_name
    try:
        local_tz = ZoneInfo(tz_name)
        local_dt = utc_time.to_datetime(timezone=local_tz)
        time_text = local_dt.strftime("%Y-%m-%d %H:%M:%S %Z")
    except Exception:
        time_text = utc_time.to_datetime().strftime("%Y-%m-%d %H:%M:%S UTC")

    painter.setPen(QColor(180, 180, 180))
    painter.setFont(text_font)
    painter.drawText(QPoint(10, 20), time_text)

    # ---- City name ----
    city_name_text = viewer_data.city_name.title()
    painter.drawText(QPoint(10, 40), city_name_text)

    # ---- View direction（Alt/Az）----
    alt_deg, az_deg = viewer_data.view_center

    def az_to_compass(az: float) -> str:
        names = [
            "N", "NNE", "NE", "ENE",
            "E", "ESE", "SE", "SSE",
            "S", "SSW", "SW", "WSW",
            "W", "WNW", "NW", "NNW",
        ]
        idx = int(((az % 360) + 11.25) // 22.5) % 16
        return names[idx]

    compass = az_to_compass(az_deg)
    deg = "\N{DEGREE SIGN}"
    view_text = f"View: Alt {alt_deg:.0f}{deg}  Az {az_deg:.0f}{deg} ({compass})"
    painter.drawText(QPoint(10, 60), view_text)

    # ---- Star/planet highlight ----
    if highlighted_object:
        obj, pos = highlighted_object
        painter.setPen(QPen(TEXT_COLOR, 2))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawEllipse(pos, 10, 10)

        # PlanetBody(dataclass) or star(dict)
        name = getattr(obj, "name", "") if hasattr(obj, "name") else obj.get("name", "")
        painter.setPen(TEXT_COLOR)
        painter.drawText(QPointF(pos.x() + 15, pos.y() - 15), str(name))


def get_screen_geometry(width: int, height: int, alt: float) -> ScreenGeometry:
    """Calculate the center and radius for drawing based on window size and view altitude."""
    margin_x = 10
    margin_y = 10
    radius = (width - margin_x * 2) // 2
    ud = 90
    dd = alt
    center = int(radius + margin_x), int((height - margin_y * 2) * ud / (ud + dd) + margin_y)
    return ScreenGeometry(center, radius)


def normalized_to_screen_xy(nx: float, ny: float, geometry: ScreenGeometry) -> QPointF:
    """Convert normalized coordinates to screen coordinates."""
    return QPointF(geometry.center[0] + nx * geometry.radius, geometry.center[1] + ny * geometry.radius)
