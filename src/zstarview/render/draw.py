import math
from typing import List, Optional, Tuple

import numpy as np
from PIL import Image
from zoneinfo import ZoneInfo

from PySide6.QtCore import QPoint, QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPainterPath, QPen, QPolygonF, QRadialGradient

from ..paths import (
    CELESTIAL_EQUATOR_COLOR,
    DIRECTIONS,
    ECLIPTIC_COLOR,
    FIELD_OF_VIEW_DEG,
    HORIZON_LINE_COLOR,
    TEXT_COLOR,
)
from ..types import ScreenGeometry, CelestialData, ViewerData, CelestialObject
from ..astro import altaz_to_normalized_xy, is_in_fov, calculate_moon_render_data
from ..utils.image import generate_moon_phase_image
from ..utils.qt import pil2qpixmap

DEBUG_ECLIPSE = False


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


def altaz_to_normalized_xy_vectorized(alt: np.ndarray, az: np.ndarray, view_center: Tuple[float, float]) -> Tuple[np.ndarray, np.ndarray]:
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


def normalized_to_screen_xy_vectorized(nx: np.ndarray, ny: np.ndarray, geometry: ScreenGeometry) -> Tuple[np.ndarray, np.ndarray]:
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
    celestial_data: Optional[CelestialData],
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
        celestial_data: A CelestialData object containing information about celestial objects.
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

    if not celestial_data:
        return None

    # Handle stars first (vectorized)
    stars = celestial_data.stars
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
    for body in celestial_data.planets:
        if not body.is_visible:
            continue
        nx, ny = altaz_to_normalized_xy(body.alt, body.az, viewer_data.view_center)
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        dist_sq = (mouse_pos.x() - px) ** 2 + (mouse_pos.y() - py) ** 2
        if dist_sq < min_dist_sq:
            min_dist_sq = dist_sq
            highlighted_object = (body, QPointF(px, py))  # body is already an object/dataclass

    return highlighted_object


def draw_outlined_text(
    painter: QPainter,
    text: str,
    pos: QPointF,
    font: QFont,
    text_color: QColor = QColor(255, 255, 255),
    outline_color: QColor = QColor.fromRgbF(0, 0, 0, 0.3),
    outline_width: float = 3.0,
):
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)

    path = QPainterPath()
    path.addText(pos, font, text)

    pen = QPen(outline_color, outline_width)
    painter.setPen(pen)
    painter.drawPath(path)

    painter.fillPath(path, text_color)

    painter.restore()


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
        if not fragments[-1] or dist(p, fragments[-1][-1]) < 0.2:
            fragments[-1].append(p)
        else:
            fragments.append([p])
    return fragments


def draw_sky_reference_lines(painter: QPainter, geometry: ScreenGeometry, celestial_data: CelestialData) -> None:
    """
    Draw celestial reference lines like the equator, ecliptic, and horizon.

    Args:
        painter: The QPainter to use for drawing.
        geometry: The screen geometry for coordinate conversion.
        celestial_data: The data containing the points for the reference lines.
    """
    point_list_pen_styles: List[Tuple[List[Tuple[float, float]], Tuple[QColor, int, List[int]]]] = [
        (celestial_data.celestial_equator_points, (CELESTIAL_EQUATOR_COLOR, 1, [8, 4])),
        (celestial_data.ecliptic_points, (ECLIPTIC_COLOR, 1, [3, 3])),
        (celestial_data.horizon_points, (HORIZON_LINE_COLOR, 1, [10, 1])),
    ]

    painter.save()
    for points, (color, width, style) in point_list_pen_styles:
        for frag in split_by_gaps(points):
            if len(frag) < 2:
                continue
            pts = [QPointF(*normalized_to_screen_xy(nx, ny, geometry)) for nx, ny in frag]
            poly = QPolygonF(pts)

            bc = QColor.fromRgbF(0, 0, 0, 0.3)
            base = QPen(bc, width + 2, Qt.PenStyle.SolidLine)
            base.setCosmetic(True)
            painter.setPen(base)
            painter.drawPolyline(poly)

            fg = QPen(QColor(*color), width)
            fg.setCosmetic(True)
            fg.setDashPattern(style)
            painter.setPen(fg)
            painter.drawPolyline(poly)
    painter.restore()


def draw_stars(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
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
        celestial_data: The data containing star information (positions, magnitudes, etc.).
        viewer_data: The viewer's data, including the view center.
        star_base_radius: A base radius for scaling star sizes.
    """
    stars = celestial_data.stars

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
            rects = [QRectF(float(cx) - float(s) / 2.0, float(cy) - float(s) / 2.0, float(s), float(s)) for cx, cy, s in zip(sx[m], sy[m], sL[m])]
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
    x, y = normalized_to_screen_xy(nx, ny, geometry)

    s = 7

    painter.setPen(QPen(QColor(*TEXT_COLOR), 1))
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
    moon_img_pil = generate_moon_phase_image(
        img_size, sun_dir_in_moon_frame, view_dir, tint_color=tint_rgba
    )

    if abs(screen_rotation_deg) > 0.1:
        moon_img_pil = moon_img_pil.rotate(
            screen_rotation_deg,
            resample=Image.Resampling.BICUBIC,
            expand=False,
            fillcolor=(0, 0, 0, 0),  # keep transparency outside the rotated bounds
        )

    pixmap = pil2qpixmap(moon_img_pil)  # should handle RGBA → QPixmap with alpha
    target_rect = QRectF(center.x() - img_size / 2, center.y() - img_size / 2, img_size, img_size)

    painter.save()
    painter.setOpacity(opacity)
    # If you ever see halos, try: painter.setCompositionMode(QPainter.CompositionMode_SourceOver)
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
    celestial_data: CelestialData,
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
        celestial_data: The data containing planet information.
        viewer_data: The viewer's data for position calculations.
        enlarge_moon: A boolean indicating whether to draw the moon larger.
        emoji_font: The QFont to use for drawing planet symbols.
    """
    moon_zoom = 5 if enlarge_moon else 1
    moon_body = None
    sun_altaz: Optional[Tuple[float, float]] = None
    moon_altaz: Optional[Tuple[float, float]] = None

    for body in celestial_data.planets:
        if body.name == "sun":
            sun_altaz = (body.alt, body.az)
        elif body.name == "moon":
            moon_altaz = (body.alt, body.az)
            moon_body = body

    text_color = QColor(*TEXT_COLOR)

    for body in celestial_data.planets:
        if not body.is_visible:
            continue

        pos = QPointF(
            *normalized_to_screen_xy(
                *altaz_to_normalized_xy(body.alt, body.az, viewer_data.view_center),
                geometry,
            )
        )

        if body.name == "sun":
            draw_gauge_cross(painter, text_color, pos)

        elif body.name == "moon" and moon_body and sun_altaz and moon_altaz:
            sun_dir_in_moon_frame, screen_rotation_deg = calculate_moon_render_data(sun_altaz, moon_altaz, viewer_data.view_center)
            moon_radius_px = (0.25 / 90.0) * geometry.radius * moon_zoom

            eclipse = body.lunar_eclipse_info
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
            draw_gauge_cross(painter, text_color, pos)

        else:
            draw_outlined_text(painter, body.symbol, pos, emoji_font, text_color)


def draw_direction_labels(painter: QPainter, geometry: ScreenGeometry, view_center: Tuple[float, float], text_font: QFont) -> None:
    """
    Draw compass direction labels (N, S, E, W) on the horizon.

    Args:
        painter: The QPainter for drawing.
        geometry: The screen geometry for coordinate conversion.
        view_center: The current view center to determine which labels are visible.
        text_font: The QFont to use for the labels.
    """
    text_color = QColor(*TEXT_COLOR)
    text_color.setAlphaF(0.7)
    painter.setPen(text_color)
    painter.setFont(text_font)
    alt = 0.0
    for label, az in DIRECTIONS.items():
        if not is_in_fov(alt, az, view_center):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        pos = QPointF(*normalized_to_screen_xy(nx, ny, geometry))
        painter.drawText(pos, label)


def draw_overlay_info(
    painter: QPainter,
    celestial_data: CelestialData,
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
        celestial_data: The data containing the current time.
        viewer_data: The data containing viewer location and view direction.
        vmag_limit: The current visual magnitude limit for stars.
        enlarge_moon: A boolean indicating if the moon is enlarged.
        highlighted_object: The currently highlighted object, if any.
        text_font: The QFont to use for the text.
    """
    text_color = QColor(*TEXT_COLOR)

    line_spacing = QFontMetrics(text_font).lineSpacing()
    line_height = int(line_spacing * 1.2)
    line_x = line_spacing
    line_y = 0

    def print_line(message: str):
        nonlocal line_x, line_y
        line_y += line_height
        draw_outlined_text(painter, message, QPointF(line_x, line_y), text_font, text_color)

    # ---- Local time ----
    utc_time = celestial_data.time
    tz_name = viewer_data.timezone_name
    try:
        local_tz = ZoneInfo(tz_name)
        local_dt = utc_time.to_datetime(timezone=local_tz)
        time_text = local_dt.strftime("%Y-%m-%d %H:%M:%S %Z")
    except Exception:
        time_text = utc_time.to_datetime().strftime("%Y-%m-%d %H:%M:%S UTC")

    print_line(time_text)

    # ---- City, view direction (Alt/Az) ----
    city_name_text = viewer_data.city_name.title()
    print_line(city_name_text)

    alt_deg, az_deg = viewer_data.view_center

    def az_to_compass(az: float) -> str:
        """Converts azimuth in degrees to a compass direction string."""
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
    print_line(view_text)

    print_line(f"Vmag limit {vmag_limit:.1f}")

    if enlarge_moon:
        print_line("Moon size: 5x")

    # ---- Star/planet highlight ----
    if highlighted_object:
        obj, pos = highlighted_object
        painter.setPen(QPen(text_color, 2))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawEllipse(pos, 10, 10)

        # PlanetBody(dataclass) or star(dict)
        name = getattr(obj, "name", "") if hasattr(obj, "name") else obj.get("name", "")
        draw_outlined_text(painter, name, QPointF(pos.x() + 15, pos.y() - 15), text_font, text_color)


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
    center = (int(radius + margin_x), int((height - margin_y * 2) * ud / (ud + dd) + margin_y))
    return ScreenGeometry(center, radius)


def normalized_to_screen_xy(nx: float, ny: float, geometry: ScreenGeometry) -> Tuple[float, float]:
    """
    Convert normalized coordinates to screen coordinates.

    Args:
        nx: Normalized x-coordinate (-1 to 1).
        ny: Normalized y-coordinate (-1 to 1).
        geometry: The screen geometry.

    Returns:
        The corresponding screen coordinates.
    """
    return geometry.center[0] + nx * geometry.radius, geometry.center[1] + ny * geometry.radius
