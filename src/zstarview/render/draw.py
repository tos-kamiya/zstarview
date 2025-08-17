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
from ..types import ScreenGeometry, SkyData, ViewerData
from ..astro import altaz_to_normalized_xy, is_in_fov, calculate_moon_render_data
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


def normalized_to_screen_xy_vectorized(
    nx: np.ndarray, ny: np.ndarray, geometry: ScreenGeometry
) -> Tuple[np.ndarray, np.ndarray]:
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
        if not body.is_visible:
            continue
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


def draw_stars(
    painter: QPainter,
    geometry: ScreenGeometry,
    sky_data: SkyData,
    viewer_data: ViewerData,
    star_base_radius: float,
):
    """Draw stars using vectorized calculations.
    Brightness is represented primarily by area (from visual magnitude),
    with additive blending. Small stars are batched as rectangles; large
    stars are drawn as soft disks with a radial gradient but with the
    same area to keep the magnitude→area mapping consistent.
    """
    stars = sky_data.stars

    # 1) mag → relative luminance → pixel area
    #    Base: L_raw = 10^(-0.4 * (vmag - v_ref))
    #    Then apply a tone curve (beta < 1) to tame very bright stars.
    nx, ny = altaz_to_normalized_xy_vectorized(stars["alt"], stars["az"], viewer_data.view_center)
    x, y = normalized_to_screen_xy_vectorized(nx, ny, geometry)

    vmag = stars["vmag"]
    bv_clamped = np.nan_to_num(stars["bv"], nan=0.45)
    rgb_colors = bv_to_rgb_vectorized(bv_clamped)  # assumes 0–255

    v_ref = 1.0  # reference mag
    L_raw = 10.0 ** (-0.4 * (vmag - v_ref))
    beta = 0.8  # 0.7–1.0: lower → less blow-up for very bright stars
    L = np.power(L_raw, beta)

    # Area scale: keep visuals relatively stable vs. widget radius.
    # base_linear_px is a "baseline linear size (px)" which we square to get area.
    base_linear_px = (geometry.radius / 500.0) * float(star_base_radius)
    base_area_px = max(0.5, base_linear_px * base_linear_px)  # avoid < 1 px²
    area_px = base_area_px * L

    # Clamp area to stabilize density and avoid extremes.
    min_area_px = 0.5  # ~1–2 px² improves visibility
    max_area_px = (geometry.radius * 0.03) ** 2  # size-dependent upper bound
    area_px = np.clip(area_px, min_area_px, max_area_px)

    # 2) Split by *linear* size while preserving the area semantics:
    #    small stars: square with side s = sqrt(area)
    #    large stars: soft disk with radius r = sqrt(area/pi)
    small_linear_px = np.sqrt(area_px)  # square side length
    large_radius_px = np.sqrt(area_px / np.pi)  # disk radius

    # Threshold is on *linear* size.
    large_star_threshold_px = 4.0  # ≥ ~4 px → draw as soft disk
    small_star_mask = small_linear_px < large_star_threshold_px
    large_star_mask = ~small_star_mask

    painter.setPen(Qt.PenStyle.NoPen)

    # 3) Large stars: soft disk with the same area, via radial gradient.
    if np.any(large_star_mask):
        lx = x[large_star_mask]
        ly = y[large_star_mask]
        lr = large_radius_px[large_star_mask]
        lrgb = rgb_colors[large_star_mask]

        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
        # Peak alpha is constant (total added light is governed by area).
        alpha_peak = 220  # tune in ~190–255
        for i in range(len(lx)):
            pos = QPointF(float(lx[i]), float(ly[i]))
            r = float(lr[i])
            base = QColor(int(lrgb[i][0]), int(lrgb[i][1]), int(lrgb[i][2]))
            c0 = QColor(base)
            c0.setAlpha(alpha_peak)
            c1 = QColor(base)
            c1.setAlpha(0)
            g = QRadialGradient(pos, r)
            g.setColorAt(0.0, c0)
            g.setColorAt(1.0, c1)
            painter.setBrush(g)
            painter.drawEllipse(pos, r, r)

    # 4) Small stars: fixed-but-gently-scaled alpha + area (squares), batched.
    if np.any(small_star_mask):
        sx = x[small_star_mask]; sy = y[small_star_mask]
        sA = area_px[small_star_mask]
        sL = np.sqrt(sA)  # side length from area
        srgb = rgb_colors[small_star_mask]

        # Base alpha kept modest so faint stars don't pop as dots.
        # Then lift α slightly with area (subtle, gamma<1).
        alpha_base  = 140                    # try 130–160
        alpha_gain  = 50                     # how much to lift at the upper end (30–60)
        alpha_scale = np.power(np.clip(sA / 4.0, 0.0, 1.0), 0.75)  # area vs. 2x2px reference
        alpha_f     = np.clip(alpha_base + alpha_gain * alpha_scale, 60, 210)

        # Quantize alpha to reduce unique RGBA buckets → faster & less banding
        # e.g., to ~8–12 steps:
        step = 12
        alpha_u8 = (np.round(alpha_f / step) * step).astype(np.uint8)

        # Batch by RGBA
        rgba = np.column_stack([srgb, alpha_u8])
        unique_rgba = np.unique(rgba, axis=0)

        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
        for r, g, b, a in unique_rgba:
            m = np.all(rgba == (r, g, b, a), axis=1)  # mask within small-star subset
            color = QColor(int(r), int(g), int(b), int(a))
            painter.setBrush(color)
            rects = [
                QRectF(float(cx) - float(s)/2.0, float(cy) - float(s)/2.0, float(s), float(s))
                for cx, cy, s in zip(sx[m], sy[m], sL[m])
            ]
            painter.drawRects(rects)

    # Reset composition mode.
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
    sun_dir_in_moon_frame: np.ndarray,
    screen_rotation_deg: float,
    opacity: float = 1.0,
):
    """Draw the moon with phase based on the sun's direction vector and apply screen rotation."""
    img_size = int(radius * 2)
    if img_size < 5:
        img_size = 5

    # The moon image is generated from a viewpoint looking along the Z-axis.
    view_dir = np.array([0, 0, 1])
    moon_img_pil = generate_moon_phase_image(img_size, sun_dir_in_moon_frame, view_dir)

    # Rotate the image to account for the projection's effect on orientation.
    if abs(screen_rotation_deg) > 0.1:
        moon_img_pil = moon_img_pil.rotate(screen_rotation_deg, resample=Image.Resampling.BICUBIC, expand=False)

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
    moon_zoom = 5 if enlarge_moon else 1
    sun_altaz: Optional[Tuple[float, float]] = None
    moon_altaz: Optional[Tuple[float, float]] = None
    moon_body = None  # Store the whole moon object

    for body in sky_data.planets:
        if body.name == "sun":
            sun_altaz = (body.alt, body.az)
        if body.name == "moon":
            moon_altaz = (body.alt, body.az)
            moon_body = body

    for body in sky_data.planets:
        if not body.is_visible:
            continue

        pos = normalized_to_screen_xy(
            *altaz_to_normalized_xy(body.alt, body.az, viewer_data.view_center), geometry
        )
        if body.name == "sun":
            draw_gauge_cross(painter, TEXT_COLOR, pos)
        elif body.name == "moon" and moon_body:  # Check moon_body exists
            sun_dir_in_moon_frame, screen_rotation_deg = calculate_moon_render_data(
                sun_altaz, moon_altaz, viewer_data.view_center
            )

            # moon's angular radius is ~0.25 deg. Screen radius is 90 deg.
            moon_radius_px = (0.25 / 90.0) * geometry.radius * moon_zoom

            # Draw the moon with its phase
            draw_moon(
                painter,
                pos,
                moon_radius_px,
                sun_dir_in_moon_frame=sun_dir_in_moon_frame,
                screen_rotation_deg=screen_rotation_deg,
                opacity=1.0 if not enlarge_moon else 0.7,
            )

            # --- Draw Lunar Eclipse Shadow ---
            eclipse_info = moon_body.eclipse_info

            if eclipse_info and eclipse_info.is_eclipse:
                # Get screen coordinates of the shadow's center
                shadow_pos = normalized_to_screen_xy(
                    *altaz_to_normalized_xy(
                        eclipse_info.shadow_center_alt,
                        eclipse_info.shadow_center_az,
                        viewer_data.view_center,
                    ),
                    geometry,
                )

                # Calculate shadow radii in pixels, scaling relative to the moon's angular radius
                if eclipse_info.moon_radius_deg > 1e-6:
                    px_per_deg = moon_radius_px / eclipse_info.moon_radius_deg
                    penumbra_radius_px = (
                        eclipse_info.penumbra_radius_deg * px_per_deg
                    )
                    umbra_radius_px = eclipse_info.umbra_radius_deg * px_per_deg

                    painter.save()
                    painter.setPen(Qt.PenStyle.NoPen)

                    # Draw Penumbra (larger, lighter shadow)
                    penumbra_color = QColor(0, 0, 0, 100)
                    painter.setBrush(penumbra_color)
                    painter.drawEllipse(
                        shadow_pos, penumbra_radius_px, penumbra_radius_px
                    )

                    # Draw Umbra (smaller, darker shadow)
                    umbra_color = QColor(40, 10, 10, 220)  # Dark reddish
                    painter.setBrush(umbra_color)
                    painter.drawEllipse(shadow_pos, umbra_radius_px, umbra_radius_px)

                    painter.restore()

            draw_gauge_cross(painter, TEXT_COLOR, pos)
        else:
            painter.setFont(emoji_font)
            painter.setPen(TEXT_COLOR)
            painter.drawText(pos, body.symbol)


def draw_direction_labels(
    painter: QPainter, geometry: ScreenGeometry, view_center: Tuple[float, float], text_font: QFont
):
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
    vmag_limit: float,
    enlarge_moon: bool,
    highlighted_object: Optional[Tuple[Dict[str, Union[str, float]], QPointF]],
    text_font: QFont,
):
    line_height = 20
    line_x = 10
    line_y = 0

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
    line_y += line_height
    painter.drawText(QPoint(line_x, line_y), time_text)

    # ---- City, view direction（Alt/Az）----
    city_name_text = viewer_data.city_name.title()
    line_y += line_height
    painter.drawText(QPoint(line_x, line_y), city_name_text)

    alt_deg, az_deg = viewer_data.view_center

    def az_to_compass(az: float) -> str:
        names = [
            "N",
            "NNE",
            "NE",
            "ENE",
            "E",
            "ESE",
            "SE",
            "SSE",
            "S",
            "SSW",
            "SW",
            "WSW",
            "W",
            "WNW",
            "NW",
            "NNW",
        ]
        idx = int(((az % 360) + 11.25) // 22.5) % 16
        return names[idx]

    compass = az_to_compass(az_deg)
    deg = "\N{DEGREE SIGN}"
    view_text = f"Alt {alt_deg:.0f}{deg}  Az {az_deg:.0f}{deg} ({compass})"
    line_y += line_height
    painter.drawText(QPoint(line_x, line_y), view_text)

    line_y += line_height
    painter.drawText(QPoint(line_x, line_y), f"Vmag limit {vmag_limit:.1f}")

    if enlarge_moon:
        line_y += line_height
        painter.drawText(QPoint(line_x, line_y), "Moon size: 5x")

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