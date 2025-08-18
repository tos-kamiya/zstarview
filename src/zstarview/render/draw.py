import math
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from PIL import Image
from zoneinfo import ZoneInfo

from PySide6.QtCore import QPoint, QPointF, QRect, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPainterPath, QPen, QPolygonF, QRadialGradient, QImage

from ..paths import (
    CELESTIAL_EQUATOR_COLOR,
    DIRECTIONS,
    ECLIPTIC_COLOR,
    FIELD_OF_VIEW_DEG,
    HORIZON_LINE_COLOR,
    TEXT_COLOR,
)
from ..types import ScreenGeometry, SkyData, ViewerData, CelestialObject
from ..astro import altaz_to_normalized_xy, is_in_fov, calculate_moon_render_data
from ..utils.image import generate_moon_phase_image
from ..utils.qt import pil2qpixmap

DEBUG_ECLIPSE = True


def bv_to_rgb_vectorized(bv: np.ndarray) -> np.ndarray:
    """
    Vectorized conversion of B-V color index to RGB tuples.

    Args:
        bv: A NumPy array of B-V color indices.

    Returns:
        A NumPy array of RGB color tuples, where each tuple is of type int.
    """
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
    """
    Vectorized conversion of altitude/azimuth to normalized screen coordinates.

    This function projects spherical coordinates (altitude and azimuth) onto a 2D
    plane using a stereographic projection. The projection is centered on the
    `view_center` coordinates.

    Args:
        alt: A NumPy array of altitude values in degrees.
        az: A NumPy array of azimuth values in degrees.
        view_center: A tuple containing the (altitude, azimuth) of the view center in degrees.

    Returns:
        A tuple of two NumPy arrays (nx, ny), representing the normalized x and y
        coordinates on the screen. The coordinates are in the range [-1, 1].
    """
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
    """
    Vectorized conversion of normalized coordinates to screen coordinates.

    This function maps normalized coordinates (in the range [-1, 1]) to the
    actual pixel coordinates on the screen, based on the provided screen geometry.

    Args:
        nx: A NumPy array of normalized x coordinates.
        ny: A NumPy array of normalized y coordinates.
        geometry: A ScreenGeometry object containing the screen's center and radius.

    Returns:
        A tuple of two NumPy arrays (x, y), representing the screen coordinates.
    """
    return (geometry.center[0] + nx * geometry.radius, geometry.center[1] + ny * geometry.radius)


def find_highlighted_object(
    sky_data: Optional[SkyData],
    viewer_data: ViewerData,
    mouse_pos: QPoint,
    geometry: ScreenGeometry,
) -> Optional[Tuple[CelestialObject, QPointF]]:
    """
    Find the nearest celestial object to the mouse cursor.

    This function searches for the closest star or planet to the given mouse
    position. It performs a vectorized search for stars and a scalar search for
    planets to optimize performance.

    Args:
        sky_data: A SkyData object containing information about celestial objects.
        viewer_data: A ViewerData object containing the viewer's location and time.
        mouse_pos: The QPoint representing the current mouse cursor position.
        geometry: A ScreenGeometry object for coordinate transformations.

    Returns:
        A tuple containing the highlighted celestial object (as a dictionary or
        a CelestialObject) and its screen position (QPointF), or None if no
        object is close enough to the cursor.
    """
    min_dist_sq = 30**2  # squared pixels
    highlighted_object: Optional[Tuple[CelestialObject, QPointF]] = None

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
            highlighted_star: CelestialObject = {key: val[closest_star_idx] for key, val in stars.items()}
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


def draw_radial_background(painter: QPainter, rect: QRectF, geometry: ScreenGeometry) -> None:
    """
    Draws a radial gradient background to represent the sky.

    This function creates a subtle radial gradient to simulate the sky's
    appearance, with a darker tone towards the zenith and a slightly
    lighter tone towards the horizon.

    Args:
        painter: The QPainter object to use for drawing.
        rect: The QRectF defining the area to fill with the background.
        geometry: A ScreenGeometry object for calculating gradient parameters.
    """
    assert geometry.radius >= 10
    fov_middle = 90 + (FIELD_OF_VIEW_DEG / 2 - 90) / 2
    r90 = float(geometry.radius)
    r_fov = float(geometry.radius * (fov_middle / 90))
    r_max = float(r_fov * 1.4)
    step_px = 0.5

    def pos(r: float) -> float:
        return max(0.0, min(1.0, r / r_max))

    def col(r: float, s: float) -> QColor:
        return QColor(0, 0, 0, max(0, 255 - (s + int(150 * (r - r90) / r_max))))

    c = geometry.center
    g = QRadialGradient(QPointF(c[0], c[1]), r_max)
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
    """
    Split a polyline by large gaps to avoid drawing long, straight lines
    across the screen when a celestial path wraps around.

    Args:
        points: A list of (x, y) tuples representing the polyline.

    Returns:
        A list of polyline fragments, where each fragment is a list of points.
    """

    def dist(p1: Tuple[float, float], p2: Tuple[float, float]) -> float:
        return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

    fragments: List[List[Tuple[float, float]]] = [[]]
    for p in points:
        if not fragments[-1] or dist(p, fragments[-1][-1]) < 0.3:
            fragments[-1].append(p)
        else:
            fragments.append([p])
    return fragments


def draw_sky_reference_lines(painter: QPainter, geometry: ScreenGeometry, sky_data: SkyData) -> None:
    """
    Draw celestial reference lines like the equator, ecliptic, and horizon.

    Args:
        painter: The QPainter to use for drawing.
        geometry: The screen geometry for coordinate conversion.
        sky_data: The data containing the points for the reference lines.
    """
    point_list_pen_styles: List[Tuple[List[Tuple[float, float]], Tuple[QColor, int, Qt.PenStyle]]] = [
        (sky_data.celestial_equator_points, (CELESTIAL_EQUATOR_COLOR, 2, Qt.PenStyle.DashLine)),
        (sky_data.ecliptic_points, (ECLIPTIC_COLOR, 2, Qt.PenStyle.DotLine)),
        (sky_data.horizon_points, (HORIZON_LINE_COLOR, 2, Qt.PenStyle.SolidLine)),
    ]
    for points, pen_style in point_list_pen_styles:
        for frag in split_by_gaps(points):
            if len(frag) >= 2:
                pts = [normalized_to_screen_xy(nx, ny, geometry) for nx, ny in frag]
                poly = QPolygonF(pts)
                painter.setPen(QPen(pen_style[0], pen_style[1], pen_style[2]))
                painter.drawPolyline(poly)


def draw_stars(
    painter: QPainter,
    geometry: ScreenGeometry,
    sky_data: SkyData,
    viewer_data: ViewerData,
    star_base_radius: float,
) -> None:
    """
    Draw stars using vectorized calculations for efficiency.

    Brightness is represented primarily by area (calculated from visual magnitude),
    with additive blending for a more realistic effect. Small stars are batched
    and drawn as rectangles, while larger stars are rendered as soft disks with
    a radial gradient, maintaining a consistent magnitude-to-area mapping.

    Args:
        painter: The QPainter object for drawing.
        geometry: The screen geometry for coordinate conversions.
        sky_data: The data containing star information (positions, magnitudes, etc.).
        viewer_data: The viewer's data, including the view center.
        star_base_radius: A base radius for scaling star sizes.
    """
    stars = sky_data.stars

    # 1) mag -> relative luminance -> pixel area
    #    Base: L_raw = 10^(-0.4 * (vmag - v_ref))
    #    Then apply a tone curve (beta < 1) to tame very bright stars.
    nx, ny = altaz_to_normalized_xy_vectorized(stars["alt"], stars["az"], viewer_data.view_center)
    x, y = normalized_to_screen_xy_vectorized(nx, ny, geometry)

    vmag = stars["vmag"]
    bv_clamped = np.nan_to_num(stars["bv"], nan=0.45)
    rgb_colors = bv_to_rgb_vectorized(bv_clamped)  # assumes 0-255

    v_ref = 1.0  # reference mag
    L_raw = 10.0 ** (-0.4 * (vmag - v_ref))
    beta = 0.8  # 0.7-1.0: lower -> less blow-up for very bright stars
    L = np.power(L_raw, beta)

    # Area scale: keep visuals relatively stable vs. widget radius.
    # base_linear_px is a "baseline linear size (px)" which we square to get area.
    base_linear_px = (geometry.radius / 500.0) * float(star_base_radius)
    base_area_px = max(0.5, base_linear_px * base_linear_px)  # avoid < 1 px^2
    area_px = base_area_px * L

    # Clamp area to stabilize density and avoid extremes.
    min_area_px = 0.5  # ~1-2 px^2 improves visibility
    max_area_px = (geometry.radius * 0.03) ** 2  # size-dependent upper bound
    area_px = np.clip(area_px, min_area_px, max_area_px)

    # 2) Split by *linear* size while preserving the area semantics:
    #    small stars: square with side s = sqrt(area)
    #    large stars: soft disk with radius r = sqrt(area/pi)
    small_linear_px = np.sqrt(area_px)  # square side length
    large_radius_px = np.sqrt(area_px / np.pi)  # disk radius

    # Threshold is on *linear* size.
    large_star_threshold_px = 4.0  # >= ~4 px -> draw as soft disk
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
        alpha_peak = 220  # tune in ~190-255
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
        sx = x[small_star_mask]
        sy = y[small_star_mask]
        sA = area_px[small_star_mask]
        sL = np.sqrt(sA)  # side length from area
        srgb = rgb_colors[small_star_mask]

        # Base alpha kept modest so faint stars don't pop as dots.
        # Then lift alpha slightly with area (subtle, gamma<1).
        alpha_base = 140  # try 130-160
        alpha_gain = 50  # how much to lift at the upper end (30-60)
        alpha_scale = np.power(np.clip(sA / 4.0, 0.0, 1.0), 0.75)  # area vs. 2x2px reference
        alpha_f = np.clip(alpha_base + alpha_gain * alpha_scale, 60, 210)

        # Quantize alpha to reduce unique RGBA buckets -> faster & less banding
        # e.g., to ~8-12 steps:
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
                QRectF(float(cx) - float(s) / 2.0, float(cy) - float(s) / 2.0, float(s), float(s))
                for cx, cy, s in zip(sx[m], sy[m], sL[m])
            ]
            painter.drawRects(rects)

    # Reset composition mode.
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)


def draw_gauge_cross(painter: QPainter, color: QColor, center: QPointF) -> None:
    """
    Draws a cross-shaped gauge marker.

    This is used to indicate the position of certain celestial objects, like
    the Sun or the Moon.

    Args:
        painter: The QPainter to use for drawing.
        color: The color of the cross.
        center: The center point (QPointF) of the cross.
    """
    cross_outer_len, cross_inner_len = 15, 4
    x, y = center.x(), center.y()
    painter.setPen(QPen(color, 1))
    painter.drawLine(QPointF(x - cross_outer_len, y), QPointF(x - cross_inner_len, y))
    painter.drawLine(QPointF(x + cross_inner_len, y), QPointF(x + cross_outer_len, y))
    painter.drawLine(QPointF(x, y - cross_outer_len), QPointF(x, y - cross_inner_len))
    painter.drawLine(QPointF(x, y + cross_inner_len), QPointF(x, y + cross_outer_len))


def draw_zenith_marker(painter: QPainter, geometry: ScreenGeometry, view_center: Tuple[float, float]) -> None:
    """
    Draws a marker at the zenith (the point directly overhead).

    Args:
        painter: The QPainter to use for drawing.
        geometry: The screen geometry for coordinate conversion.
        view_center: The current view center (altitude, azimuth).
    """
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
    base_color: Optional[QColor] = None,
) -> None:
    """
    Draw the moon with its correct phase and orientation.

    The moon's phase is determined by the sun's direction vector relative to
    the moon. The image is also rotated to match the screen orientation.

    Args:
        painter: The QPainter for drawing.
        center: The center position of the moon on the screen.
        radius: The radius of the moon in pixels.
        sun_dir_in_moon_frame: The direction vector of the sun in the moon's reference frame.
        screen_rotation_deg: The rotation angle of the screen in degrees.
        opacity: The opacity of the moon image.
        base_color: An optional QColor to tint the moon, used for eclipses.
    """
    img_size = max(5, int(math.ceil(radius * 2.0)))
    view_dir = np.array([0, 0, 1], dtype=float)
    if base_color is not None:
        tint_rgba = (base_color.red(), base_color.green(), base_color.blue(), base_color.alpha())
    else:
        tint_rgba = None
    moon_img_pil = generate_moon_phase_image(img_size, sun_dir_in_moon_frame, view_dir, tint_color=tint_rgba)

    if abs(screen_rotation_deg) > 0.1:
        moon_img_pil = moon_img_pil.rotate(screen_rotation_deg, resample=Image.Resampling.BICUBIC, expand=False)

    pixmap = pil2qpixmap(moon_img_pil)
    target_rect = QRectF(center.x() - img_size / 2, center.y() - img_size / 2, img_size, img_size)

    painter.save()
    painter.setOpacity(opacity)
    painter.drawPixmap(target_rect, pixmap, QRectF(0, 0, img_size, img_size))

    if base_color is not None:
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceAtop)
        painter.setBrush(base_color)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawEllipse(center, radius, radius)

    painter.restore()


def draw_planets(
    painter: QPainter,
    geometry: ScreenGeometry,
    sky_data: SkyData,
    viewer_data: ViewerData,
    enlarge_moon: bool,
    emoji_font: QFont,
) -> None:
    """
    Draw the planets, including the Sun and Moon.

    This function iterates through the planets in the sky data, calculates their
    screen positions, and draws them. The Sun is drawn as a gauge cross, the
    Moon is drawn with its phase, and other planets are represented by emoji symbols.

    Args:
        painter: The QPainter for drawing.
        geometry: The screen geometry for coordinate conversion.
        sky_data: The data containing planet information.
        viewer_data: The viewer's data for position calculations.
        enlarge_moon: A boolean indicating whether to draw the moon larger.
        emoji_font: The QFont to use for drawing planet symbols.
    """
    moon_zoom = 5 if enlarge_moon else 1
    moon_body = None
    sun_altaz: Optional[Tuple[float, float]] = None
    moon_altaz: Optional[Tuple[float, float]] = None

    for body in sky_data.planets:
        if body.name == "sun":
            sun_altaz = (body.alt, body.az)
        elif body.name == "moon":
            moon_altaz = (body.alt, body.az)
            moon_body = body

    for body in sky_data.planets:
        if not body.is_visible:
            continue

        pos = normalized_to_screen_xy(
            *altaz_to_normalized_xy(body.alt, body.az, viewer_data.view_center),
            geometry,
        )

        if body.name == "sun":
            draw_gauge_cross(painter, TEXT_COLOR, pos)

        elif body.name == "moon" and moon_body and sun_altaz and moon_altaz:
            sun_dir_in_moon_frame, screen_rotation_deg = calculate_moon_render_data(
                sun_altaz, moon_altaz, viewer_data.view_center
            )
            moon_radius_px = (0.25 / 90.0) * geometry.radius * moon_zoom

            eclipse = body.eclipse_info
            base_color: Optional[QColor] = None
            if eclipse and eclipse.is_eclipse:
                if eclipse.eclipse_type == "partial":
                    base_color = QColor(30, 0, 0, 60)
                elif eclipse.eclipse_type == "total":
                    base_color = QColor(40, 10, 10, 180)

            draw_moon(
                painter,
                pos,
                moon_radius_px,
                sun_dir_in_moon_frame=sun_dir_in_moon_frame,
                screen_rotation_deg=screen_rotation_deg,
                opacity=1.0 if not enlarge_moon else 0.7,
                base_color=base_color,
            )
            draw_gauge_cross(painter, TEXT_COLOR, pos)

        else:
            painter.setFont(emoji_font)
            painter.setPen(TEXT_COLOR)
            painter.drawText(pos, body.symbol)


def draw_direction_labels(
    painter: QPainter, geometry: ScreenGeometry, view_center: Tuple[float, float], text_font: QFont
) -> None:
    """
    Draw compass direction labels (N, S, E, W) on the horizon.

    Args:
        painter: The QPainter for drawing.
        geometry: The screen geometry for coordinate conversion.
        view_center: The current view center to determine which labels are visible.
        text_font: The QFont to use for the labels.
    """
    painter.setPen(TEXT_COLOR)
    painter.setFont(text_font)
    alt = 0.0
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
    highlighted_object: Optional[Tuple[CelestialObject, QPointF]],
    text_font: QFont,
) -> None:
    """
    Draws overlay text information on the screen.

    This includes the current time, location, view direction, magnitude limit,
    and information about any highlighted celestial object.

    Args:
        painter: The QPainter for drawing.
        sky_data: The data containing the current time.
        viewer_data: The data containing viewer location and view direction.
        vmag_limit: The current visual magnitude limit for stars.
        enlarge_moon: A boolean indicating if the moon is enlarged.
        highlighted_object: The currently highlighted object, if any.
        text_font: The QFont to use for the text.
    """
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

    # ---- City, view direction (Alt/Az) ----
    city_name_text = viewer_data.city_name.title()
    line_y += line_height
    painter.drawText(QPoint(line_x, line_y), city_name_text)

    alt_deg, az_deg = viewer_data.view_center

    def az_to_compass(az: float) -> str:
        """Converts azimuth in degrees to a compass direction string."""
        names = [
            "N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
            "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW",
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
    """
    Calculate the center and radius for drawing based on window size and view altitude.

    Args:
        width: The width of the drawing area.
        height: The height of the drawing area.
        alt: The altitude of the view center.

    Returns:
        A ScreenGeometry object with the calculated center and radius.
    """
    margin_x = 10
    margin_y = 10
    radius = (width - margin_x * 2) // 2
    ud = 90.0
    dd = alt
    center = (
        int(radius + margin_x),
        int((height - margin_y * 2) * ud / (ud + dd) + margin_y)
    )
    return ScreenGeometry(center, radius)


def normalized_to_screen_xy(nx: float, ny: float, geometry: ScreenGeometry) -> QPointF:
    """
    Convert normalized coordinates to screen coordinates.

    Args:
        nx: Normalized x-coordinate (-1 to 1).
        ny: Normalized y-coordinate (-1 to 1).
        geometry: The screen geometry.

    Returns:
        The corresponding screen coordinates as a QPointF.
    """
    return QPointF(geometry.center[0] + nx * geometry.radius, geometry.center[1] + ny * geometry.radius)


def _deg2rad(d: float) -> float:
    """Converts degrees to radians."""
    return d * math.pi / 180.0

def _clamp01(x: float) -> float:
    """Clamps a value to the range [0, 1]."""
    return max(0.0, min(1.0, x))

def _angle_between(alt1_deg: float, az1_deg: float, alt2_deg: float, az2_deg: float) -> float:
    """
    Calculates the angular distance between two points on the celestial sphere.

    Args:
        alt1_deg: Altitude of the first point in degrees.
        az1_deg: Azimuth of the first point in degrees.
        alt2_deg: Altitude of the second point in degrees.
        az2_deg: Azimuth of the second point in degrees.

    Returns:
        The angular distance in radians.
    """
    a1, z1 = _deg2rad(alt1_deg), _deg2rad(az1_deg)
    a2, z2 = _deg2rad(alt2_deg), _deg2rad(az2_deg)
    cos_g = (math.sin(a1)*math.sin(a2) + math.cos(a1)*math.cos(a2)*math.cos(z2 - z1))
    cos_g = max(-1.0, min(1.0, cos_g))
    return math.acos(cos_g)

# --- ▼ Phenomenon Imitation Model ▼ ---

def get_sun_color(sun_alt_deg: float) -> Tuple[float, float, float]:
    """
    Step 1: Determine the color of sunlight based on the sun's altitude.

    Args:
        sun_alt_deg: The sun's altitude in degrees.

    Returns:
        A tuple of (r, g, b) float values representing the sun's color.
    """
    # Define colors
    zenith_color = (0.3, 0.5, 1.0)  # Color at zenith (blue)
    horizon_color = (1.0, 0.8, 0.4) # Color at horizon (orange)
    night_color = (0.01, 0.02, 0.05) # Night color (dark blue)

    # Normalize sun altitude from -5 degrees (sunset) to 90 degrees (zenith) to a 0-1 range
    t = _clamp01((sun_alt_deg + 5.0) / 95.0)

    # Daytime color (horizon to zenith)
    day_color = _lerp_color(horizon_color, zenith_color, t)

    # Mix day and night colors (to represent the rapid darkening near the horizon)
    fade = _clamp01(sun_alt_deg / 10.0) # Fade between 0 and 10 degrees

    return _lerp_color(night_color, day_color, fade)


def get_sky_color(view_altaz: Tuple[float, float], sun_altaz: Tuple[float, float]) -> Tuple[float, float, float]:
    """
    Calculates the sky color for a given viewing direction and sun position.

    Args:
        view_altaz: A tuple of (altitude, azimuth) for the viewing direction.
        sun_altaz: A tuple of (altitude, azimuth) for the sun's position.

    Returns:
        A tuple of (r, g, b) float values for the sky color.
    """
    # After sunset, it gradually darkens (completely dark at -10°, towards day at 0°)
    sun_alt_deg, sun_az_deg = sun_altaz
    view_alt_deg, view_az_deg = view_altaz

    if sun_alt_deg <= -10.0:
        return (0.0, 0.0, 0.0)

    # Basic sun color (assumed 0..1: preferably linear RGB)
    sun_color = get_sun_color(sun_alt_deg)

    # 1) Brightness based on angle to the sun (stable even at zenith)
    gamma = _angle_between(view_alt_deg, view_az_deg, sun_alt_deg, sun_az_deg)  # [0..pi]
    cosg = math.cos(gamma)
    brightness = (1.0 + cosg) * 0.5          # 0..1
    brightness = brightness ** 2.0            # Emphasize the sun-facing direction

    # 2) Tone based on altitude (darker at zenith, whitish at horizon)
    t = _clamp01(view_alt_deg / 90.0)         # 0(horizon) -> 1(zenith)
    zenith_darkness = 0.5 + 0.5 * t           # 0.5..1.0 (darker towards zenith)
    horizon_whiteness = (1.0 - t) * 0.3       # 0.3..0.0 (whiter towards horizon)

    # 3) Twilight correction (interpolate -10..0° to 0..1)
    if sun_alt_deg < 0.0:
        twilight = _smoothstep(-10.0, 0.0, sun_alt_deg)  # 0..1
    else:
        twilight = 1.0

    # Composite (simple: mix of white + sun contribution)
    r = sun_color[0] * brightness * zenith_darkness * twilight + horizon_whiteness * twilight
    g = sun_color[1] * brightness * zenith_darkness * twilight + horizon_whiteness * twilight
    b = sun_color[2] * brightness * zenith_darkness * twilight + horizon_whiteness * twilight

    # 4) Clip (final safety)
    return (_clamp01(r), _clamp01(g), _clamp01(b))

def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    """
    Performs a smooth Hermite interpolation between 0 and 1 when edge0 < x < edge1.
    """
    # Hermite smoothstep
    t = _clamp01((x - edge0) / (edge1 - edge0))
    return t * t * (3.0 - 2.0 * t)


def _lerp_color(a: Tuple[float,float,float],
                b: Tuple[float,float,float], t: float) -> Tuple[float,float,float]:
    """
    Linearly interpolates between two colors.
    """
    return (a[0] + (b[0]-a[0]) * t,
            a[1] + (b[1]-a[1]) * t,
            a[2] + (b[2]-a[2]) * t)

def _fwd_offset_altaz(alt_c: float, az_c: float,
                      theta_deg: float, psi_deg: float) -> Tuple[float, float]:
    """
    Calculates the (alt, az) of a point that is at an angular distance 'theta'
    and direction 'psi' from a center point (alt_c, az_c).

    Args:
        alt_c: Center altitude in degrees.
        az_c: Center azimuth in degrees.
        theta_deg: Angular distance from the center in degrees.
        psi_deg: Direction from the center in degrees (0° North, 90° East, clockwise).

    Returns:
        A tuple of (altitude, azimuth) in degrees for the new point.
    """
    phi1  = math.radians(alt_c)
    lambda1  = math.radians(az_c)
    theta   = math.radians(theta_deg)
    psi   = math.radians(psi_deg)

    sin_phi1, cos_phi1 = math.sin(phi1), math.cos(phi1)
    sin_theta, cos_theta   = math.sin(theta), math.cos(theta)

    sin_phi2 = sin_phi1 * cos_theta + cos_phi1 * sin_theta * math.cos(psi)
    sin_phi2 = max(-1.0, min(1.0, sin_phi2))
    phi2 = math.asin(sin_phi2)

    y = math.sin(psi) * sin_theta * cos_phi1
    x = cos_theta - sin_phi1 * sin_phi2
    lambda2 = lambda1 + math.atan2(y, x)

    alt = math.degrees(phi2)
    az  = (math.degrees(lambda2) + 360.0) % 360.0
    return alt, az

def _alpha_from_alt(alt: float, alpha: float, fade_hi: float = 0.0, fade_lo: float = -2.0) -> float:
    """
    Returns an alpha value based on altitude, for fading near the horizon.
       alt >= fade_hi  -> alpha (equivalent to 1)
       fade_lo < alt < fade_hi -> linear fade
       alt <= fade_lo -> 0 (not drawn)

    Args:
        alt: The altitude in degrees.
        alpha: The base alpha value.
        fade_hi: The altitude above which alpha is at its maximum.
        fade_lo: The altitude below which alpha is zero.

    Returns:
        The calculated alpha value.
    """
    if alt <= fade_lo:
        return 0.0
    if alt >= fade_hi:
        return alpha
    t = (alt - fade_lo) / (fade_hi - fade_lo)   # [fade_hi, fade_lo] -> [0,1]
    return alpha * t

def draw_sky_color_disc(
    painter: QPainter,
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],   # (alt, az)
    sun_altaz: Tuple[float, float],     # (alt, az)
    *,
    exposure: float = 1.0,
    saturation: float = 1.2,
    alpha: float = 1.0,
    # --- Sampling density (knobs for quality vs. speed) ---
    ring_step_px: int = 6,            # Target pixel interval in the radial direction (basis for determining Δθ)
    sample_pitch_px: float = 6.0,     # Target sample interval on each ring (px)
    min_ang_samples: int = 8,         # Minimum number of samples for each ring
    dot_size: int = 12,               # Side length of the tile (square) in px: roughly 2*ring_step_px is a guideline
    # --- Parameters for stabilizing Δθ estimation ---
    deriv_probe_deg: float = 0.25,    # Small angle (degrees) for finite difference of dr/dθ
    min_theta_step_deg: float = 0.2,  # Lower limit for Δθ (degrees)
    max_theta_step_deg: float = 6.0,  # Upper limit for Δθ (degrees)
) -> None:
    """
    Draws the sky color disc with dynamic sampling.

    The radial step (Δθ) is dynamically determined so that the on-screen
    pixel distance is around `ring_step_px`. The number of samples in the
    circumferential direction is determined by the "measured radius" of the
    ring and the sample pitch.

    Args:
        painter: The QPainter to draw with.
        geometry: The screen geometry.
        view_center: The (alt, az) of the view center.
        sun_altaz: The (alt, az) of the sun.
        exposure: Exposure adjustment for the final color.
        saturation: Saturation adjustment for the final color.
        alpha: Overall alpha transparency.
        ring_step_px: Target pixel distance for radial steps.
        sample_pitch_px: Target pixel distance for circumferential samples.
        min_ang_samples: Minimum number of circumferential samples.
        dot_size: The size of the colored dots to draw.
        deriv_probe_deg: Probe angle for derivative estimation.
        min_theta_step_deg: Minimum angular step.
        max_theta_step_deg: Maximum angular step.
    """
    assert altaz_to_normalized_xy and normalized_to_screen_xy and get_sky_color

    R = int(geometry.radius)
    if R < 2:
        return

    # Extract center (supports QPoint/QPointF/Tuple)
    c = geometry.center
    cx = c.x() if hasattr(c, "x") else c[0]
    cy = c.y() if hasattr(c, "y") else c[1]

    img_w = img_h = R * 2
    img = QImage(img_w, img_h, QImage.Format.Format_ARGB32)
    img.fill(QColor(0, 0, 0, 0))  # Transparent

    ip = QPainter(img)
    ip.setRenderHint(QPainter.RenderHint.Antialiasing, False)
    ip.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, False)

    # Clip to a circle
    path = QPainterPath()
    path.addEllipse(0, 0, 2*R, 2*R)  # Circle at (0,0) with radius R in image coords
    ip.setClipPath(path)

    # Clamp to avoid zenith singularity
    EPS = 0.01
    alt_c, az_c = view_center
    if alt_c <= -90.0: alt_c = -(90.0 - EPS)
    if alt_c >=  90.0: alt_c =  (90.0 - EPS)

    # The outer circumference is always 90° to match the normalization spec
    theta_max = 90.0

    # Function to "measure" the local radius r_px(θ)
    def screen_radius_px(theta_deg: float) -> float:
        alt_p, az_p = _fwd_offset_altaz(alt_c, az_c, theta_deg, 0.0)
        nx, ny = altaz_to_normalized_xy(alt_p, az_p, (alt_c, az_c))
        # Convert to image coordinates relative to image center (R, R)
        s = normalized_to_screen_xy(nx, ny, geometry)
        sx, sy = s.x(), s.y()
        return max(0.0, math.hypot(sx - cx, sy - cy))

    # Advance angle from 0 to 90° (Δθ is dynamic)
    theta = 0.0
    half = max(1, int(dot_size // 2))
    while True:
        # Current ring radius (px)
        r_px = screen_radius_px(theta)

        # Number of circumferential samples
        circumference = max(1.0, 2.0 * math.pi * r_px)
        n_ang = max(min_ang_samples, int(round(circumference / max(1.0, float(sample_pitch_px)))))

        # Fill the ring
        for i in range(n_ang):
            psi_deg = (360.0 * i) / n_ang  # 0° North, 90° East, clockwise

            # (θ,ψ) -> (alt,az)
            alt, az = _fwd_offset_altaz(alt_c, az_c, theta, psi_deg)

            # (alt,az) -> screen -> image coordinates
            nx, ny = altaz_to_normalized_xy(alt, az, (alt_c, az_c))
            s = normalized_to_screen_xy(nx, ny, geometry)
            # Convert to image coordinates (origin at top-left of the QImage)
            xi = int(round(s.x() - (cx - R)))
            yi = int(round(s.y() - (cy - R)))


            if xi < 0 or xi >= img_w or yi < 0 or yi >= img_h:
                continue

            aa = _alpha_from_alt(alt, alpha, fade_hi=1.0, fade_lo=-1.0)
            if aa <= 0.0:
                continue

            rr, gg, bb = get_sky_color((alt, az), sun_altaz)
            gray = rr*0.299 + gg*0.587 + bb*0.114
            rr, gg, bb = _lerp_color((gray, gray, gray), (rr, gg, bb), saturation)
            rr *= exposure; gg *= exposure; bb *= exposure
            rr = _clamp01(rr); gg = _clamp01(gg); bb = _clamp01(bb)

            ip.fillRect(xi - half, yi - half, 2*half, 2*half, QColor.fromRgbF(rr, gg, bb, aa))

        # Termination condition (ensure the 90° ring is always drawn)
        if theta >= theta_max - 1e-6:
            break

        # ---- Determine Δθ from screen distance: Δθ ≈ ring_step_px / (dr_px/dθ) ----
        # Estimate dr/dθ using finite differences
        probe = min(deriv_probe_deg, theta_max - theta)
        r_next = screen_radius_px(theta + probe)
        dr_dtheta = (r_next - r_px) / max(1e-6, probe)

        if dr_dtheta <= 1e-6:
            dtheta = max_theta_step_deg  # Safe side (almost no change / numerical error at edges)
        else:
            dtheta = float(ring_step_px) / dr_dtheta

        # Limit the angle step
        dtheta = max(min_theta_step_deg, min(max_theta_step_deg, dtheta))

        theta += dtheta

    ip.end()

    # Paste the generated image onto the main painter
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    target_rect = QRect(int(cx - R), int(cy - R), 2*R, 2*R)
    painter.drawImage(target_rect, img)
    painter.restore()