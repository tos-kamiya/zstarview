import math
from typing import List, Optional, Tuple

import numpy as np
from PIL import Image
from zoneinfo import ZoneInfo

from PySide6.QtCore import QPoint, QPointF, QRect, QRectF, Qt
from PySide6.QtGui import QImage, QColor, QFont, QFontMetrics, QLinearGradient, QPainter, QPainterPath, QPen, QPolygonF, QRadialGradient

from ..paths import (
    BACKGROUND_FIELD_OF_VIEW_DEG,
    CELESTIAL_EQUATOR_COLOR,
    DIRECTIONS,
    ECLIPTIC_COLOR,
    HORIZON_LINE_COLOR,
    TEXT_COLOR,
    STATUS_LINE_COLOR,
)
from ..types import ScreenGeometry, CelestialData, ViewerData, CelestialObject, PlanetBody
from ..astro import altaz_to_normalized_xy, is_in_fov, calculate_moon_render_data
from ..utils.image import generate_moon_phase_image
from ..utils.qt import pil2qpixmap

DEBUG_ECLIPSE = False


def get_text_style(preset: str = "night") -> Tuple[QColor, QColor]:
    """Return (text_color, outline_color) tuned for the selected visual preset."""
    if preset == "white":
        return QColor(18, 29, 48), QColor(245, 250, 255, 210)
    if preset == "day":
        return QColor(22, 33, 52), QColor(238, 245, 255, 200)
    if preset == "black":
        return QColor(246, 249, 255), QColor(2, 2, 3, 236)
    return QColor(*TEXT_COLOR), QColor.fromRgbF(0, 0, 0, 0.3)


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


def effective_star_draw_vmag_limit(base_limit: float, interaction_mode: bool, interaction_cap: float = 7.0) -> float:
    """Return the effective Vmag limit for the current draw pass."""
    if not interaction_mode:
        return float(base_limit)
    return float(min(base_limit, interaction_cap))


def flare_strength_from_vmag(vmag: float) -> float:
    """Return monotonic flare strength [0, 1] for all magnitudes."""
    vmag_bright = -1.5
    vmag_faint = 6.0
    t = (vmag_faint - float(vmag)) / (vmag_faint - vmag_bright)
    t = max(0.0, min(1.0, t))
    return t**1.35


def compute_flare_profile(vmag: float, core_radius_px: float) -> Tuple[float, float]:
    """Compute (core_scale, flare_outer_px) for a star.

    `flare_outer_px` is the additional radial reach from the core radius.
    If `flare_outer_px < 1.0`, flare should not be drawn.
    """
    strength = flare_strength_from_vmag(vmag)
    flare_outer_px = float(core_radius_px) * (0.65 * strength)
    if flare_outer_px < 1.0:
        return 1.0, 0.0
    core_scale = 1.0 / math.sqrt(1.0 + 0.9 * strength)
    return core_scale, flare_outer_px


def planet_disc_style_from_vmag(vmag: Optional[float]) -> Tuple[float, int]:
    """Return (radius_px, alpha) for a planet disc marker."""
    if vmag is None or not math.isfinite(float(vmag)):
        return 3.0, 200
    # Keep explicit clipping here for readability: very bright planets saturate at -1.5.
    vmag_clipped = float(np.clip(float(vmag), -1.5, 6.0))
    # Reuse star-like brightness mapping so brighter planets appear stronger.
    strength = flare_strength_from_vmag(vmag_clipped)
    radius_px = 2.4 + 3.2 * strength
    alpha = int(np.clip(round(125 + 130 * strength), 110, 255))
    return radius_px, alpha


def planet_bloom_profile_from_vmag(vmag: Optional[float], core_radius_px: float) -> Tuple[float, int, int]:
    """Return bloom profile as (radius_px, center_alpha, mid_alpha).

    Planet disk size remains clipped at -1.5 mag via `planet_disc_style_from_vmag`,
    while this bloom profile uses the raw magnitude to express extra brightness
    for very bright planets (e.g. Venus).
    """
    if vmag is None or not math.isfinite(float(vmag)):
        return 0.0, 0, 0

    vm = float(np.clip(float(vmag), -6.0, 6.0))
    # Base term follows existing clipped brightness behavior.
    base = flare_strength_from_vmag(float(np.clip(vm, -1.5, 6.0)))
    # Extra term activates only when brighter than -1.5.
    extra = max(0.0, min(1.0, (-1.5 - vm) / 4.5))
    strength = 0.60 * base + 0.40 * extra

    if strength < 0.12:
        return 0.0, 0, 0

    r_core = max(1.0, float(core_radius_px))
    # Keep bloom close to VR look (halo roughly up to ~2.4x core radius).
    bloom_radius = r_core * (1.20 + 1.20 * strength)
    center_alpha = int(np.clip(round(10 + 58 * strength), 8, 90))
    mid_alpha = int(np.clip(round(4 + 34 * strength), 3, 60))
    return bloom_radius, center_alpha, mid_alpha


def planet_marker_color(name: str) -> QColor:
    """Return display color for a planet marker."""
    palette = {
        "mercury": QColor(190, 190, 182),
        "venus": QColor(245, 226, 176),
        "mars": QColor(232, 126, 96),
        "jupiter": QColor(224, 188, 141),
        "saturn": QColor(226, 214, 154),
        "uranus": QColor(157, 224, 218),
        "neptune": QColor(108, 152, 234),
        "pluto": QColor(194, 166, 132),
    }
    return QColor(palette.get(name, QColor(*TEXT_COLOR)))


def altaz_to_normalized_xy_vectorized(
    alt_deg: np.ndarray,
    az_deg: np.ndarray,
    view_center_altaz_deg: Tuple[float, float],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Vectorized conversion of altitude/azimuth to normalized screen coordinates.

    This function projects spherical coordinates (altitude and azimuth) onto a 2D
    plane using a stereographic projection. The projection is centered on the
    `view_center` coordinates.

    Args:
        alt_deg: A NumPy array of altitude values in degrees.
        az_deg: A NumPy array of azimuth values in degrees.
        view_center_altaz_deg: `(view_alt_deg, view_az_deg)` in degrees.

    Returns:
        A tuple `(nx, ny)` of normalized screen coordinates.
        The normalization is `r_norm = angular_distance_deg / 90`, so values can
        exceed `1.0` when the angular distance is larger than 90 degrees.
    """
    center_alt, center_az = view_center_altaz_deg
    alt1, az1 = np.radians(center_alt), np.radians(center_az)
    alt2, az2 = np.radians(alt_deg), np.radians(az_deg)

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
        # For planets, allow below-horizon display as long as they are in the current FOV.
        if not is_in_fov(body.alt, body.az, viewer_data.view_center):
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


def draw_radial_background(
    painter: QPainter,
    rect: QRectF,
    geometry: ScreenGeometry,
    *,
    preset: str = "night",
) -> None:
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
    fov_middle = 90 + (BACKGROUND_FIELD_OF_VIEW_DEG - 90) / 2
    r90 = float(geometry.radius)
    r_fov = float(geometry.radius * (fov_middle / 90))
    r_max = float(r_fov * 1.4)
    step_px = 0.5

    def pos(r: float) -> float:
        return max(0.0, min(1.0, r / r_max))

    if preset == "white":
        def col(r: float, s: float) -> QColor:
            t = max(0.0, min(1.0, (r - r90) / max(1.0, r_max - r90)))
            gray = int(246 - 54 * t)
            aa = max(0, 248 - (s + int(120 * t)))
            return QColor(gray, gray, gray, aa)
    elif preset == "day":
        def col(r: float, s: float) -> QColor:
            t = max(0.0, min(1.0, (r - r90) / max(1.0, r_max - r90)))
            # Keep day palette bright, but slightly blue-leaning.
            rr = int(230 - 28 * t)
            gg = int(242 - 34 * t)
            bb = int(255 - 34 * t)
            # Keep "day" visibly lighter/more transparent than "white".
            aa = max(0, 170 - (s + int(120 * t)))
            return QColor(rr, gg, bb, aa)
    elif preset == "night":
        def col(r: float, s: float) -> QColor:
            t = max(0.0, min(1.0, (r - r90) / max(1.0, r_max - r90)))
            rr = int(10 - 7 * t)
            gg = int(12 - 9 * t)
            bb = int(16 - 11 * t)
            aa = max(0, 236 - (s + int(95 * t)))
            return QColor(rr, gg, bb, aa)
    elif preset == "black":
        def col(r: float, s: float) -> QColor:
            t = max(0.0, min(1.0, (r - r90) / max(1.0, r_max - r90)))
            gray = int(4 - 3 * t)
            aa = max(0, 255 - (s + int(45 * t)))
            return QColor(gray, gray, gray, aa)
    else:
        def col(r: float, s: float) -> QColor:
            t = max(0.0, min(1.0, (r - r90) / max(1.0, r_max - r90)))
            return QColor(0, 0, 0, max(0, 205 - (s + int(130 * t))))

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

            # Use a tinted outline (same hue) instead of black to avoid dark gaps
            # on bright or highly transparent backgrounds.
            base_color = QColor(*color)
            base_color.setAlpha(70)
            base = QPen(base_color, width + 2, Qt.PenStyle.SolidLine)
            base.setCosmetic(True)
            # Keep the outline only under visible dash segments.
            # This avoids black "gaps" on bright/transparent backgrounds.
            base.setDashPattern(style)
            base.setCapStyle(Qt.PenCapStyle.RoundCap)
            base.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(base)
            painter.drawPolyline(poly)

            fg = QPen(QColor(*color), width)
            fg.setCosmetic(True)
            fg.setDashPattern(style)
            fg.setCapStyle(Qt.PenCapStyle.RoundCap)
            fg.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(fg)
            painter.drawPolyline(poly)
    painter.restore()


def draw_stars(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    star_base_radius: float,
    *,
    visibility_boost: float = 1.0,
    draw_vmag_limit: Optional[float] = None,
    viewport_size: Tuple[int, int] | None = None,
) -> None:
    """
    Draw stars using a numpy canvas that paints uniformly sized rectangles.

    The visual area of each star is derived from its magnitude (a 1st-magnitude
    star maps to a 6×6 rectangle), and the color intensity is reduced for stars
    whose calculated area falls below 1 pixel². The resulting RGB canvas is
    converted into a `QImage` and blended additively to the painter before
    drawing other objects.

    Args:
        painter: The QPainter object for drawing.
        geometry: The screen geometry for coordinate conversions.
        celestial_data: The data containing star information (positions, magnitudes, etc.).
        viewer_data: The viewer's data, including the view center.
        star_base_radius: A base radius for scaling star sizes.
        visibility_boost: Multiplier for the overall star colors.
        draw_vmag_limit: Optional override to skip fainter stars entirely.
        viewport_size: Optional `(width, height)` of the drawing area, used to create the numpy canvas; if omitted, defaults to the painter's clip rect.
    """
    stars = celestial_data.stars
    visibility_boost = float(np.clip(visibility_boost, 0.7, 2.0))

    if draw_vmag_limit is not None:
        draw_mask = stars["vmag"] <= float(draw_vmag_limit)
        if not np.any(draw_mask):
            return
        alt = stars["alt"][draw_mask]
        az = stars["az"][draw_mask]
        bv = stars["bv"][draw_mask]
        size_factor = stars["size_factor"][draw_mask]
        color_factor_base = stars["color_factor_base"][draw_mask]
    else:
        alt = stars["alt"]
        az = stars["az"]
        bv = stars["bv"]
        size_factor = stars["size_factor"]
        color_factor_base = stars["color_factor_base"]

    width_px = None
    height_px = None
    if viewport_size:
        width_px = max(1, int(viewport_size[0]))
        height_px = max(1, int(viewport_size[1]))
    else:
        width_px = painter.viewport().width()
        height_px = painter.viewport().height()

    if width_px <= 0 or height_px <= 0:
        return

    # Convert directions to pixel centers.
    nx, ny = altaz_to_normalized_xy_vectorized(alt, az, viewer_data.view_center)
    x, y = normalized_to_screen_xy_vectorized(nx, ny, geometry)

    bv_clamped = np.nan_to_num(bv, nan=0.45)
    rgb_colors = bv_to_rgb_vectorized(bv_clamped).astype(np.float32)

    max_size = max(12.0, float(max(1.0, star_base_radius)))
    size_float = float(star_base_radius) * size_factor
    size_px = np.clip(np.round(size_float), 1, int(max_size)).astype(int)

    color_factor = np.clip(0.1 + 0.9 * color_factor_base, 0.0, 1.0) * visibility_boost
    color_factor = np.clip(color_factor, 0.0, 1.0)
    star_colors = rgb_colors * color_factor[:, None]

    canvas = np.zeros((height_px, width_px, 3), dtype=np.float32)

    ix = np.round(x).astype(int)
    iy = np.round(y).astype(int)

    half = (size_px // 2).astype(int)
    x0 = ix - half
    y0 = iy - half
    x1 = x0 + size_px
    y1 = y0 + size_px

    x0_clamped = np.clip(x0, 0, width_px)
    y0_clamped = np.clip(y0, 0, height_px)
    x1_clamped = np.clip(x1, 0, width_px)
    y1_clamped = np.clip(y1, 0, height_px)

    valid = (x1_clamped > x0_clamped) & (y1_clamped > y0_clamped) & (size_px > 0)
    indices = np.nonzero(valid)[0]
    for idx in indices:
        canvas[y0_clamped[idx]:y1_clamped[idx], x0_clamped[idx]:x1_clamped[idx], :] += star_colors[idx]

    np.clip(canvas, 0.0, 255.0, out=canvas)
    canvas_uint8 = np.ascontiguousarray(canvas.astype(np.uint8))
    if canvas_uint8.size == 0:
        return

    alpha = (np.any(canvas_uint8 > 0, axis=2)).astype(np.uint8) * 255
    rgba = np.zeros((height_px, width_px, 4), dtype=np.uint8)
    rgba[:, :, :3] = canvas_uint8
    rgba[:, :, 3] = alpha
    image = QImage(rgba.data, width_px, height_px, width_px * 4, QImage.Format_RGBA8888).copy()
    painter.save()
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
    painter.drawImage(0, 0, image)

    painter.restore()
    # Reset composition mode.
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)


def draw_gauge_cross(
    painter: QPainter,
    color: QColor,
    center: QPointF,
    *,
    scale: float = 1.0,
    pen_width: float = 1.0,
) -> None:
    """
    Draws a cross-shaped gauge marker.

    This is used to indicate the position of certain celestial objects, like
    the Sun or the Moon.

    Args:
        painter: The QPainter to use for drawing.
        color: The color of the cross.
        center: The center point (QPointF) of the cross.
    """
    scale = max(0.5, float(scale))
    cross_outer_len = int(round(15 * scale))
    cross_inner_len = max(1, int(round(4 * scale)))
    x, y = center.x(), center.y()
    painter.setPen(QPen(color, float(pen_width)))
    painter.drawLine(QPointF(x - cross_outer_len, y), QPointF(x - cross_inner_len, y))
    painter.drawLine(QPointF(x + cross_inner_len, y), QPointF(x + cross_outer_len, y))
    painter.drawLine(QPointF(x, y - cross_outer_len), QPointF(x, y - cross_inner_len))
    painter.drawLine(QPointF(x, y + cross_inner_len), QPointF(x, y + cross_outer_len))


def draw_zenith_marker(painter: QPainter, geometry: ScreenGeometry, view_center: Tuple[float, float]) -> None:
    """
    Draws markers at zenith and nadir.

    Args:
        painter: The QPainter to use for drawing.
        geometry: The screen geometry for coordinate conversion.
        view_center: The current view center (altitude, azimuth).
    """
    az_ref = view_center[1]
    s = 7
    painter.setPen(QPen(QColor(*TEXT_COLOR), 1))
    for alt in (90.0, -90.0):
        if not is_in_fov(alt, az_ref, view_center):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az_ref, view_center)
        x, y = normalized_to_screen_xy(nx, ny, geometry)
        painter.drawLine(QPointF(x - s, y - s), QPointF(x + s, y + s))
        painter.drawLine(QPointF(x - s, y + s), QPointF(x + s, y - s))


def draw_moon(
    painter: QPainter,
    center: QPointF,
    radius_px: float,
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
        radius_px: The radius of the moon in pixels.
        sun_dir_in_moon_frame: The direction vector of the sun in the moon's reference frame.
        screen_rotation_deg: The rotation angle of the screen in degrees.
        opacity: The opacity of the moon image.
        base_color: An optional QColor to tint the moon, used for eclipses.
    """
    img_size = max(5, int(math.ceil(radius_px * 2.0)))
    view_dir = np.array([0, 0, 1], dtype=float)
    if base_color is not None:
        tint_rgba = (base_color.red(), base_color.green(), base_color.blue(), base_color.alpha())
    else:
        tint_rgba = None
    moon_img_pil = generate_moon_phase_image(img_size, sun_dir_in_moon_frame, view_dir, tint_color=tint_rgba)

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
        painter.drawEllipse(center, radius_px, radius_px)

    painter.restore()


def _collect_sun_moon_context(planets: List[PlanetBody]) -> Tuple[Optional[PlanetBody], Optional[Tuple[float, float]], Optional[Tuple[float, float]]]:
    """Collect moon body and sun/moon alt-az pairs from planet list."""
    moon_body: Optional[PlanetBody] = None
    sun_altaz: Optional[Tuple[float, float]] = None
    moon_altaz: Optional[Tuple[float, float]] = None
    for body in planets:
        if body.name == "sun":
            sun_altaz = (body.alt, body.az)
        elif body.name == "moon":
            moon_altaz = (body.alt, body.az)
            moon_body = body
    return moon_body, sun_altaz, moon_altaz


def _moon_eclipse_overlay_color(body: PlanetBody) -> Optional[QColor]:
    """Return overlay color for lunar eclipse rendering, if needed."""
    eclipse = body.lunar_eclipse_info
    if not eclipse or not eclipse.is_eclipse:
        return None
    if eclipse.eclipse_type == "partial":
        return QColor(30, 0, 0, 60)
    if eclipse.eclipse_type == "total":
        return QColor(40, 10, 10, 180)
    return None


def _draw_moon_planet(
    painter: QPainter,
    pos: QPointF,
    geometry: ScreenGeometry,
    body: PlanetBody,
    viewer_data: ViewerData,
    sun_altaz: Tuple[float, float],
    moon_altaz: Tuple[float, float],
    enlarge_moon: bool,
    cross_color: QColor,
) -> None:
    """Draw moon phase disc (with eclipse tint) and its gauge marker."""
    moon_zoom = 5 if enlarge_moon else 1
    sun_dir_in_moon_frame, screen_rotation_deg = calculate_moon_render_data(sun_altaz, moon_altaz, viewer_data.view_center)
    moon_radius_px = (0.25 / 90.0) * geometry.radius * moon_zoom
    draw_moon(
        painter,
        pos,
        moon_radius_px,
        sun_dir_in_moon_frame=sun_dir_in_moon_frame,
        screen_rotation_deg=screen_rotation_deg,
        opacity=1.0 if not enlarge_moon else 0.7,
        base_color=_moon_eclipse_overlay_color(body),
    )
    draw_gauge_cross(painter, cross_color, pos)


def draw_planet_disc(
    painter: QPainter,
    pos: QPointF,
    color: QColor,
    *,
    radius_px: float,
    alpha: int,
) -> None:
    """Draw a no-flare solid circular marker for planets."""
    r = max(1.5, float(radius_px))
    fill = QColor(color)
    fill.setAlpha(int(np.clip(alpha, 1, 255)))
    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(fill)
    painter.drawEllipse(pos, r, r)
    painter.restore()


def draw_planet_bloom(
    painter: QPainter,
    pos: QPointF,
    color: QColor,
    *,
    core_radius_px: float,
    vmag: Optional[float],
) -> None:
    """Draw a soft additive bloom around a planet marker."""
    bloom_radius, center_alpha, mid_alpha = planet_bloom_profile_from_vmag(vmag, core_radius_px)
    if bloom_radius <= 0.0 or center_alpha <= 0:
        return

    c0 = QColor(color)
    c0.setAlpha(center_alpha)
    c1 = QColor(color)
    c1.setAlpha(mid_alpha)
    c2 = QColor(color)
    c2.setAlpha(0)

    gradient = QRadialGradient(pos, bloom_radius)
    gradient.setColorAt(0.0, c0)
    gradient.setColorAt(0.45, c1)
    gradient.setColorAt(1.0, c2)

    painter.save()
    # Additive blend to emulate bloom-like luminous spill.
    painter.setCompositionMode(QPainter.CompositionMode_Plus)
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(gradient)
    painter.drawEllipse(pos, bloom_radius, bloom_radius)
    painter.restore()


def draw_solar_system_bodies(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    enlarge_moon: bool,
    *,
    text_font: Optional[QFont] = None,
    preset: str = "night",
) -> None:
    """
    Draw major solar system bodies (Sun, Moon, and planets).

    This function iterates through bodies in the sky data, calculates their
    screen positions, and draws them. The Sun is drawn as a gauge cross, the
    Moon is drawn with its phase. Other planets are represented by circular markers.

    Args:
        painter: The QPainter for drawing.
        geometry: The screen geometry for coordinate conversion.
        celestial_data: The data containing solar system body information.
        viewer_data: The viewer's data for position calculations.
        enlarge_moon: A boolean indicating whether to draw the moon larger.
    """
    moon_body, sun_altaz, moon_altaz = _collect_sun_moon_context(celestial_data.planets)

    # Keep body markers in a stable high-contrast color over the sky disc.
    text_color = QColor(*TEXT_COLOR)
    label_text_color, label_outline_color = get_text_style(preset)
    if text_font is not None:
        painter.setFont(text_font)
        label_font = text_font
    else:
        label_font = painter.font() if hasattr(painter, "font") else QFont()

    for body in celestial_data.planets:
        # Draw planets if they are in-view, even below horizon.
        if not is_in_fov(body.alt, body.az, viewer_data.view_center):
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
            _draw_moon_planet(
                painter,
                pos,
                geometry,
                body,
                viewer_data,
                sun_altaz,
                moon_altaz,
                enlarge_moon,
                text_color,
            )

        else:
            radius_px, alpha = planet_disc_style_from_vmag(body.vmag)
            marker_color = planet_marker_color(body.name)
            draw_planet_bloom(painter, pos, marker_color, core_radius_px=radius_px, vmag=body.vmag)
            marker_color.setAlpha(alpha)
            draw_planet_disc(painter, pos, marker_color, radius_px=radius_px, alpha=alpha)
            # Keep planet cross markers shorter than the Moon/Sun gauge marker.
            draw_gauge_cross(painter, text_color, pos, scale=0.55, pen_width=1.0)
            draw_outlined_text(
                painter,
                body.name,
                QPointF(pos.x() + 12.0, pos.y() - 10.0),
                label_font,
                label_text_color,
                label_outline_color,
            )


def draw_direction_labels(
    painter: QPainter,
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    text_font: QFont,
    *,
    preset: str = "night",
) -> None:
    """
    Draw compass direction labels and horizon markers on the horizon.

    Args:
        painter: The QPainter for drawing.
        geometry: The screen geometry for coordinate conversion.
        view_center: The current view center to determine which labels are visible.
        text_font: The QFont to use for the labels.
    """
    # Match direction labels to the horizon line color.
    text_color = QColor(*HORIZON_LINE_COLOR)
    outline_color = QColor.fromRgbF(0, 0, 0, 0.3)
    marker_color = QColor(*HORIZON_LINE_COLOR)
    marker_pen = QPen(marker_color, 1.6)
    marker_pen.setCosmetic(True)
    marker_half_len_px = 6.0  # constant visual length regardless of view altitude
    marker_hit_radius_px = 4.0
    tangent_probe_deg = 0.6
    label_outward_offset_px = 4.0
    painter.setFont(text_font)
    fm = QFontMetrics(text_font)
    alt = 0.0
    for label, az in DIRECTIONS.items():
        if not is_in_fov(alt, az, view_center):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        pos = QPointF(*normalized_to_screen_xy(nx, ny, geometry))

        # Estimate local horizon tangent from azimuth finite difference.
        az_prev = (az - tangent_probe_deg + 360.0) % 360.0
        az_next = (az + tangent_probe_deg) % 360.0
        p_prev = QPointF(*normalized_to_screen_xy(*altaz_to_normalized_xy(alt, az_prev, view_center), geometry))
        p_next = QPointF(*normalized_to_screen_xy(*altaz_to_normalized_xy(alt, az_next, view_center), geometry))
        tx = p_next.x() - p_prev.x()
        ty = p_next.y() - p_prev.y()
        t_norm = math.hypot(tx, ty)
        if t_norm <= 1e-6:
            rx = pos.x() - geometry.center[0]
            ry = pos.y() - geometry.center[1]
            tx, ty = -ry, rx
            t_norm = math.hypot(tx, ty)
        if t_norm <= 1e-6:
            tx, ty, t_norm = 1.0, 0.0, 1.0
        ux, uy = tx / t_norm, ty / t_norm
        # Draw a marker line perpendicular to the local horizon tangent.
        nxp, nyp = -uy, ux
        painter.save()
        painter.setPen(marker_pen)
        painter.drawLine(
            QPointF(pos.x() - nxp * marker_half_len_px, pos.y() - nyp * marker_half_len_px),
            QPointF(pos.x() + nxp * marker_half_len_px, pos.y() + nyp * marker_half_len_px),
        )
        painter.restore()

        label_pos = QPointF(pos)
        # Baseline-relative bounds for draw_outlined_text(path.addText baseline semantics).
        bounds = fm.tightBoundingRect(label)
        label_rect = QRectF(
            label_pos.x() + bounds.x(),
            label_pos.y() + bounds.y(),
            bounds.width(),
            bounds.height(),
        )
        # If label and marker overlap, push label slightly outward from disc center.
        nearest_x = min(max(pos.x(), label_rect.left()), label_rect.right())
        nearest_y = min(max(pos.y(), label_rect.top()), label_rect.bottom())
        dx0 = pos.x() - nearest_x
        dy0 = pos.y() - nearest_y
        overlap = (dx0 * dx0 + dy0 * dy0) <= (marker_hit_radius_px**2)
        if overlap:
            ox = pos.x() - geometry.center[0]
            oy = pos.y() - geometry.center[1]
            norm = math.hypot(ox, oy)
            if norm > 1e-6:
                label_pos = QPointF(
                    label_pos.x() + (ox / norm) * label_outward_offset_px,
                    label_pos.y() + (oy / norm) * label_outward_offset_px,
                )

        draw_outlined_text(
            painter,
            label,
            label_pos,
            text_font,
            text_color,
            outline_color,
            outline_width=2.5,
        )


def draw_overlay_info(
    painter: QPainter,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    vmag_limit: float,
    enlarge_moon: bool,
    highlighted_object: Optional[Tuple[CelestialObject, QPointF]],
    text_font: QFont,
    *,
    preset: str = "night",
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
    text_color, outline_color = get_text_style(preset)

    line_spacing = QFontMetrics(text_font).lineSpacing()
    line_height = int(line_spacing * 1.2)
    line_x = line_spacing
    line_y = 0

    def print_line(message: str):
        nonlocal line_x, line_y
        line_y += line_height
        draw_outlined_text(painter, message, QPointF(line_x, line_y), text_font, text_color, outline_color)

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
        names = tuple(DIRECTIONS.keys())
        sector = 360.0 / len(names)
        idx = int(((az % 360.0) + sector / 2.0) // sector) % len(names)
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
        draw_outlined_text(painter, name, QPointF(pos.x() + 15, pos.y() - 15), text_font, text_color, outline_color)


def draw_status_line_text(
    painter: QPainter,
    message: str,
    status_line_font: QFont,
    viewport_rect: QRect,
    *,
    preset: str = "night",
) -> None:
    """
    Draws a single-line error message at the bottom-left corner, using the same
    outlined text style as draw_overlay_info.

    Args:
        painter: QPainter to draw with.
        message: Error message text (single line).
        text_font: Font used for the overlay text.
        viewport_rect: Target drawing rect (typically the window rect()).
    """
    if not message:
        return

    if preset == "white":
        color = QColor(64, 22, 22)
        outline_color = QColor(255, 245, 245, 220)
    elif preset == "day":
        color = QColor(78, 26, 26)
        outline_color = QColor(250, 242, 242, 215)
    elif preset == "black":
        color = QColor(255, 220, 220)
        outline_color = QColor(2, 2, 3, 236)
    else:
        color = QColor(*STATUS_LINE_COLOR)
        outline_color = QColor.fromRgbF(0, 0, 0, 0.3)

    painter.save()
    painter.setFont(status_line_font)

    fm = painter.fontMetrics()
    margin = fm.lineSpacing()
    baseline_y = viewport_rect.bottom() - margin // 4
    x = margin

    draw_outlined_text(painter, "> " + message, QPointF(x, baseline_y), status_line_font, color, outline_color)
    painter.restore()


def get_screen_geometry(width_px: int, height_px: int, view_alt_deg: float) -> ScreenGeometry:
    """
    Calculate circular viewport geometry.

    Args:
        width_px: The width of the drawing area in pixels.
        height_px: The height of the drawing area in pixels.
        view_alt_deg: Unused (kept for call-site compatibility).

    Returns:
        A ScreenGeometry object with the calculated center and radius.

    Layout rule:
        - Tall/square-ish windows: centered circle that fits inside the window.
        - Wide windows (width >= height): keep the disc top tangent to the top
          margin and maximize radius so the current horizon sits near the bottom.
          Radius is also limited by horizontal space.
    """
    _ = view_alt_deg
    margin_x = 10
    margin_y = 10
    avail_w = max(2, int(width_px) - margin_x * 2)
    avail_h = max(2, int(height_px) - margin_y * 2)
    alt = max(0.0, min(90.0, float(view_alt_deg)))

    if width_px >= height_px:
        # horizon_y ~= margin + radius * (1 + alt/90)
        # Choose radius from height so horizon is near the bottom, but never
        # exceed horizontal fit.
        r_height = int(avail_h / (1.0 + alt / 90.0))
        r_width = avail_w // 2
        radius_px = max(1, min(r_width, r_height))
        center = (int(width_px) // 2, margin_y + radius_px)
    else:
        radius_px = max(1, min(avail_w // 2, avail_h // 2))
        center = (int(width_px) // 2, int(height_px) // 2)
    return ScreenGeometry(center, radius_px)


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
