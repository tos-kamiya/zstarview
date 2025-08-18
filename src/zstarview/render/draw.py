import math
from typing import Callable, Dict, List, Optional, Tuple, Union

import numpy as np
from PIL import Image
from zoneinfo import ZoneInfo

from PySide6.QtCore import QPoint, QPointF, QRect, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF, QRadialGradient, QImage

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

DEBUG_ECLIPSE = True


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
        sx = x[small_star_mask]
        sy = y[small_star_mask]
        sA = area_px[small_star_mask]
        sL = np.sqrt(sA)  # side length from area
        srgb = rgb_colors[small_star_mask]

        # Base alpha kept modest so faint stars don't pop as dots.
        # Then lift α slightly with area (subtle, gamma<1).
        alpha_base = 140  # try 130–160
        alpha_gain = 50  # how much to lift at the upper end (30–60)
        alpha_scale = np.power(np.clip(sA / 4.0, 0.0, 1.0), 0.75)  # area vs. 2x2px reference
        alpha_f = np.clip(alpha_base + alpha_gain * alpha_scale, 60, 210)

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
                QRectF(float(cx) - float(s) / 2.0, float(cy) - float(s) / 2.0, float(s), float(s))
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
    base_color: Optional[QColor] = None,
):
    """Draw the moon with phase based on the sun's direction vector and apply screen rotation."""
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
        painter.setCompositionMode(QPainter.CompositionMode_SourceAtop)
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
):
    moon_zoom = 5 if enlarge_moon else 1
    moon_body = None
    sun_altaz = None
    moon_altaz = None

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

        elif body.name == "moon" and moon_body:
            sun_dir_in_moon_frame, screen_rotation_deg = calculate_moon_render_data(
                sun_altaz, moon_altaz, viewer_data.view_center
            )
            moon_radius_px = (0.25 / 90.0) * geometry.radius * moon_zoom

            eclipse = body.eclipse_info
            base_color = None
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


def _deg2rad(d): return d * math.pi / 180.0
def _clamp01(x): return max(0.0, min(1.0, x))

def _angle_between(alt1_deg, az1_deg, alt2_deg, az2_deg) -> float:
    a1, z1 = _deg2rad(alt1_deg), _deg2rad(az1_deg)
    a2, z2 = _deg2rad(alt2_deg), _deg2rad(az2_deg)
    cos_g = (math.sin(a1)*math.sin(a2) + math.cos(a1)*math.cos(a2)*math.cos(z2 - z1))
    cos_g = max(-1.0, min(1.0, cos_g))
    return math.acos(cos_g)

# 色を混ぜる関数
def _lerp_color(c1, c2, t):
    t = _clamp01(t)
    return (
        c1[0] * (1-t) + c2[0] * t,
        c1[1] * (1-t) + c2[1] * t,
        c1[2] * (1-t) + c2[2] * t,
    )

# --- ▼ 現象模倣モデル ▼ ---

def get_sun_color(sun_alt_deg: float) -> tuple[float, float, float]:
    """
    ステップ1: 太陽の高度に基づいて太陽光の色を決める
    """
    # 色を定義
    zenith_color = (0.3, 0.5, 1.0)  # 天頂にあるときの色 (青)
    horizon_color = (1.0, 0.8, 0.4) # 地平線にあるときの色 (オレンジ)
    night_color = (0.01, 0.02, 0.05) # 夜の色 (濃紺)

    # 太陽高度 -5度(日没) から 90度(天頂) の範囲を 0-1 に正規化
    t = _clamp01((sun_alt_deg + 5.0) / 95.0)

    # 昼の色 (地平線〜天頂)
    day_color = _lerp_color(horizon_color, zenith_color, t)
    
    # 昼と夜の色を混ぜる (太陽が地平線近くで急激に暗くなるのを表現)
    fade = _clamp01((sun_alt_deg) / 10.0) # 0〜10度でフェード
    
    return _lerp_color(night_color, day_color, fade)



# EPS = 1e-3  # 天頂の特異を避けるための余白


# def clamp_alt(alt: float, alt_min: float = -60.0, alt_max: float = 89.5) -> float:
#     # alt_max は 90 より少し小さく
#     if alt < alt_min: 
#         return alt_min
#     if alt > alt_max:
#         return alt_max
#     # ぴったり90°は避ける
#     if abs(alt - 90.0) < EPS:
#         return 90.0 - EPS
#     return alt


# def screen_to_altaz_equidistant(
#     x: int, y: int,
#     geometry: ScreenGeometry,
#     view_center: tuple[float, float],  # (alt_c, az_c)
#     fov_deg: float,
# ) -> tuple[float, float]:
#     """
#     画面上の点(x,y)を、視線中心(view_center)を軸とした(alt, az)に変換する。
#     方位は 北=0°, 東=90°（時計回り）を満たす。天頂付近の不安定性を回避。
#     """
#     cx, cy = geometry.center
#     R = geometry.radius
#     dx = x - cx
#     dy = y - cy
#     rho = math.hypot(dx, dy)
#     if rho == 0 or R <= 0:
#         return view_center

#     # 平面方位 ψ: 北(上)=0°, 時計回り増加
#     psi = math.atan2(dx, -dy)

#     # 距離→角度（方位等距投影）
#     # まれな丸め誤差で rho > R になるのを防ぐ保険
#     rho_ratio = min(rho / R, 1.0)
#     theta = math.radians(rho_ratio * (fov_deg * 0.5))

#     alt_c, az_c = view_center
#     alt_c = clamp_alt(alt_c)
#     lat1 = math.radians(alt_c)
#     lon1 = math.radians(az_c)

#     sin_lat1 = math.sin(lat1)
#     cos_lat1 = math.cos(lat1)
#     sin_theta = math.sin(theta)
#     cos_theta = math.cos(theta)
#     cos_psi = math.cos(psi)
#     sin_psi = math.sin(psi)

#     # φ2（=alt2）
#     sin_lat2 = sin_lat1 * cos_theta + cos_lat1 * sin_theta * cos_psi
#     sin_lat2 = max(-1.0, min(1.0, sin_lat2))
#     lat2 = math.asin(sin_lat2)
#     cos_lat2 = math.cos(lat2)

#     # λ2（=az2）
#     if abs(cos_lat2) < 1e-6:
#         # 天頂/天底で方位は未定義 → 連続性の観点で中心の方位を保持
#         lon2 = lon1
#     else:
#         # 前進式（安定）：Δλ = atan2( sinψ sinθ cosφ1, cosθ − sinφ1 sinφ2 )
#         y_ = sin_psi * sin_theta * cos_lat1
#         x_ = cos_theta - sin_lat1 * sin_lat2
#         lon2 = lon1 + math.atan2(y_, x_)

#     alt = math.degrees(lat2)
#     az  = (math.degrees(lon2) + 360.0) % 360.0
#     return alt, az


def get_sky_color(view_altaz: Tuple[float, float], sun_altaz: Tuple[float, float]) -> tuple[float, float, float]:
    # 日没後は緩やかに暗く（-10°で完全に暗い、0°で昼側へ）
    sun_alt_deg, sun_az_deg = sun_altaz
    view_alt_deg, view_az_deg = view_altaz

    if sun_alt_deg <= -10.0:
        return (0.0, 0.0, 0.0)

    # 基本の太陽色（0..1想定：できればリニアRGB）
    sun_color = get_sun_color(sun_alt_deg)

    # 1) 太陽との角度による明るさ（天頂でも安定）
    gamma = _angle_between(view_alt_deg, view_az_deg, sun_alt_deg, sun_az_deg)  # [0..π]
    cosg = math.cos(gamma)
    brightness = (1.0 + cosg) * 0.5          # 0..1
    brightness = brightness ** 2.0            # 日向を強調

    # 2) 高度に応じたトーン（天頂は濃く、地平は白っぽく）
    t = _clamp01(view_alt_deg / 90.0)         # 0(地平) → 1(天頂)
    zenith_darkness = 0.5 + 0.5 * t           # 0.5..1.0（天頂ほど暗く）
    horizon_whiteness = (1.0 - t) * 0.3       # 0.3..0.0（地平ほど白混ぜ）

    # 3) トワイライト補正（-10..0°を0..1で補間）
    if sun_alt_deg < 0.0:
        twilight = _smoothstep(-10.0, 0.0, sun_alt_deg)  # 0..1
    else:
        twilight = 1.0

    # 合成（簡易：白のミックス + 日向寄与）
    r = sun_color[0] * brightness * zenith_darkness * twilight + horizon_whiteness * twilight
    g = sun_color[1] * brightness * zenith_darkness * twilight + horizon_whiteness * twilight
    b = sun_color[2] * brightness * zenith_darkness * twilight + horizon_whiteness * twilight

    # 4) クリップ（最終安全）
    return (_clamp01(r), _clamp01(g), _clamp01(b))

def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    # Hermite smoothstep
    t = _clamp01((x - edge0) / (edge1 - edge0))
    return t * t * (3.0 - 2.0 * t)


# def draw_sky_color_disc(
#     painter: QPainter,
#     geometry: ScreenGeometry,
#     view_center: Tuple[float, float],  # (alt, az)
#     sun_altaz: Tuple[float, float],  # (alt, az)
#     fov_deg: float,
#     *,
#     exposure: float = 1.0,
#     saturation: float = 1.2,
#     ground_color: Tuple[float, float, float] = (1.0, 0.0, 0.0),
#     alpha: float = 1.0,
# ) -> None:
#     """天球ディスク（center=geometry.center, radius=geometry.radius）を空色で塗る。
#     """
#     R = int(geometry.radius)
#     if R < 2:
#         return

#     cx, cy = geometry.center
#     buf_w = buf_h = R * 2
#     x0 = int(cx - R)
#     y0 = int(cy - R)
#     r2 = R * R

#     img = QImage(buf_w, buf_h, QImage.Format.Format_ARGB32)
#     img.fill(QColor(0, 0, 0, 0))

#     # 走査
#     by_printed = set()  # debug
#     for by in range(0, buf_h):
#         # 実画面Y（丸めの歪み低減）
#         y = y0 + by
#         dy_b = by - R
#         dy2 = dy_b * dy_b
#         for bx in range(0, buf_w):
#             dx_b = bx - R
#             if (dx_b * dx_b + dy2) > r2:
#                 continue  # ディスク外

#             # 実画面X
#             x = x0 + bx

#             # 1) 画面→(alt, az)
#             altaz = screen_to_altaz_equidistant(x, y, geometry, view_center, fov_deg)

#             # 2) 空色（リニアRGB想定）
#             r, g, b = get_sky_color(altaz, sun_altaz)

#             # 3) 彩度調整＆露出（リニア）
#             col = r * 0.299 + g * 0.587 + b * 0.114
#             r, g, b = _lerp_color((col, col, col), (r, g, b), saturation)
#             r *= exposure; g *= exposure; b *= exposure

#             # 地平線下（-5〜0°でフェードイン）
#             alt = altaz[0]
#             if by % 10 == 0 and by not in by_printed:  # debug
#                 by_printed.add(by)  # debug
#                 print(f"{view_center[0]=}, {geometry.radius=}, {by=}, {alt=}")  # debug

#             if abs(alt) <= 0.3:
#                 img.setPixel(bx, by, QColor(0, 255, 0, 255).rgba())
#                 continue

#             if alt < 0.0:
#                 # t = _clamp01((alt + 5.0) / 5.0)
#                 # r, g, b = _lerp_color(ground_color, (r, g, b), t)
#                 r, g, b = ground_color

#             # 4) クリップ＆書き込み
#             r = _clamp01(r); g = _clamp01(g); b = _clamp01(b)
#             img.setPixel(bx, by, QColor.fromRgbF(r, g, b, alpha).rgba())

#     # 描画
#     painter.save()
#     painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
#     target_rect = QRect(cx - R, cy - R, 2 * R, 2 * R)
#     painter.drawImage(target_rect, img)
#     painter.restore()


def _lerp_color(a: Tuple[float,float,float],
                b: Tuple[float,float,float], t: float) -> Tuple[float,float,float]:
    return (a[0] + (b[0]-a[0]) * t,
            a[1] + (b[1]-a[1]) * t,
            a[2] + (b[2]-a[2]) * t)

def _fwd_offset_altaz(alt_c: float, az_c: float,
                      theta_deg: float, psi_deg: float) -> Tuple[float, float]:
    """中心(alt_c,az_c)から角距離theta, 方位psi（北0°/東90°/時計回り）だけ進めた(alt,az)。"""
    φ1  = math.radians(alt_c)
    λ1  = math.radians(az_c)
    θ   = math.radians(theta_deg)
    ψ   = math.radians(psi_deg)

    sinφ1, cosφ1 = math.sin(φ1), math.cos(φ1)
    sinθ, cosθ   = math.sin(θ), math.cos(θ)

    sinφ2 = sinφ1 * cosθ + cosφ1 * sinθ * math.cos(ψ)
    sinφ2 = max(-1.0, min(1.0, sinφ2))
    φ2 = math.asin(sinφ2)

    y = math.sin(ψ) * sinθ * cosφ1
    x = cosθ - sinφ1 * sinφ2
    λ2 = λ1 + math.atan2(y, x)

    alt = math.degrees(φ2)
    az  = (math.degrees(λ2) + 360.0) % 360.0
    return alt, az

def draw_sky_color_disc(
    painter: QPainter,
    geometry,                           # ScreenGeometry(center, radius)
    view_center: Tuple[float, float],   # (alt, az)
    sun_altaz: Tuple[float, float],     # (alt, az)
    *,
    exposure: float = 1.0,
    saturation: float = 1.2,
    ground_color: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    alpha: float = 1.0,
    ring_step_px: int = 1,          # 同心円リングの間隔（おおよそのピクセル）
    sample_pitch_px: float = 1.0,   # 各リング上の目標サンプル間隔（px）
    min_ang_samples: int = 8,       # 各リングの最小サンプル数
    dot_size: int = 1,              # 打点サイズ（1=1px, 2=2x2 ...）
) -> None:
    """
    既存の altaz_to_normalized_xy(alt,az,view_center) と normalized_to_screen_xy(nx,ny,geometry)
    を使って、天球ディスクをリング×方位のサンプルで塗る。
    """
    R = int(geometry.radius)
    if R < 2:
        return

    cx, cy = geometry.center
    img_w = img_h = R * 2
    img = QImage(img_w, img_h, QImage.Format.Format_ARGB32)
    img.fill(QColor(0, 0, 0, 0))  # 透明で開始（下地を生かせる）

    theta_max = 90.0
    alt_c, az_c = view_center

    # リングを「見かけのスクリーン半径」ベースで等間隔にする
    # そのために、リングごとに1点だけ前進写像して r_px を計測 → 円周長からサンプル数を決める
    # スクリーン半径に対する角距離の近似は使わないので、既存投影に完全追従します。
    # 粗密は ring_step_px / sample_pitch_px / dot_size で調整してください。
    rho_px = 0
    while True:
        # そのリングに対応する角距離を線形対応で近似
        theta_deg = theta_max * (rho_px / R)

        # そのリングのスクリーン半径を前進写像で“実測”
        alt_probe, az_probe = _fwd_offset_altaz(alt_c, az_c, theta_deg, 0.0)
        nx_probe, ny_probe = altaz_to_normalized_xy(alt_probe, az_probe, view_center)
        s = normalized_to_screen_xy(nx_probe, ny_probe, geometry)
        sx_probe, sy_probe = s.x(), s.y()
        r_px = math.hypot(sx_probe - cx, sy_probe - cy)

        # 円周長からサンプル数を見積もり
        circumference = max(1.0, 2.0 * math.pi * r_px)
        n_ang = max(min_ang_samples, int(round(circumference / max(0.5, sample_pitch_px))))

        for i in range(n_ang):
            psi_deg = (360.0 * i) / n_ang  # 北0°/東90°/時計回り
            # 1) 中心から(θ,ψ)だけ進めた(alt,az)
            alt, az = _fwd_offset_altaz(alt_c, az_c, theta_deg, psi_deg)

            # 2) 既存前進写像でスクリーン座標へ
            nx, ny = altaz_to_normalized_xy(alt, az, view_center)
            s = normalized_to_screen_xy(nx, ny, geometry)
            sx, sy = s.x(), s.y()
            xi, yi = int(round(sx)), int(round(sy))

            # ディスク内チェック
            dx, dy = xi - cx, yi - cy
            if dx*dx + dy*dy > R*R + 1:
                continue
            if not (cx - R <= xi < cx + R and cy - R <= yi < cy + R):
                continue

            # 3) 色（地面下はベタで ground_color。透明で抜きたいなら continue に変更）
            if alt < 0.0:
                rr, gg, bb = ground_color
                aa = 1.0
            else:
                rr, gg, bb = get_sky_color((alt, az), sun_altaz)
                gray = rr*0.299 + gg*0.587 + bb*0.114
                rr, gg, bb = _lerp_color((gray, gray, gray), (rr, gg, bb), saturation)
                rr *= exposure; gg *= exposure; bb *= exposure
                rr = _clamp01(rr); gg = _clamp01(gg); bb = _clamp01(bb)
                aa = alpha

            # 4) ドットで描く（穴埋め用に dot_size>1 も可）
            half = max(0, dot_size // 2)
            for oy in range(-half, half + 1):
                yy = yi + oy
                if yy < cy - R or yy >= cy + R: 
                    continue
                for ox in range(-half, half + 1):
                    xx = xi + ox
                    if xx < cx - R or xx >= cx + R:
                        continue
                    ddx, ddy = xx - cx, yy - cy
                    if ddx*ddx + ddy*ddy > R*R + 1:
                        continue
                    img.setPixel(xx - (cx - R), yy - (cy - R),
                                 QColor.fromRgbF(rr, gg, bb, aa).rgba())

        rho_px += max(1, ring_step_px)
        if rho_px > R:
            break

    # 貼り付け
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    target_rect = QRect(cx - R, cy - R, 2*R, 2*R)
    painter.drawImage(target_rect, img)
    painter.restore()
