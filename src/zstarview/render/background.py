import math
from typing import Callable
from zoneinfo import ZoneInfo

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import (
    QColor,
    QPainter,
    QPainterPath,
    QPen,
    QPolygonF,
    QRadialGradient,
)

from ..astro import altaz_to_normalized_xy
from ..paths import (
    BACKGROUND_FIELD_OF_VIEW_DEG1,
    BACKGROUND_FIELD_OF_VIEW_DEG2,
    DIRECTIONS,
    GUI_BUTTON_SIZE,
    ThemeStyle,
)
from ..types import CelestialData, ScreenGeometry, ViewerData
from ..utils.location_display import build_location_info_lines

FRAMELESS_WINDOW_BORDER_WIDTH = 24.0
ALT_RING_HIGHLIGHT_ALPHA = 12
ALT_RING_PEN_WIDTH = 2.0
ALT_RING_DIMALT_BRIGHTNESS_THRESHOLD = 64.0
ALT_RING_DIMALT_DARKEN_FACTOR = 0.87
ALT_RING_DIMALT_LIGHTEN_MIX = 0.05
ATLAS_HORIZON_TINT_RGBA = (245, 168, 82, 255)
ATLAS_NIGHT_TINT_RGBA = (48, 52, 58, 255)
ATLAS_DAY_TINT_RGBA = (150, 200, 235, 255)


def atlas_background_tint_rgba(sun_alt_deg: float | None) -> tuple[int, int, int, int] | None:
    """Return a subtle Atlas background tint derived from the Sun altitude."""
    if sun_alt_deg is None or not math.isfinite(float(sun_alt_deg)):
        return None
    sun_alt = float(sun_alt_deg)
    if sun_alt <= 0.0:
        lower = ATLAS_NIGHT_TINT_RGBA
        upper = ATLAS_HORIZON_TINT_RGBA
        weight = (max(-6.0, sun_alt) + 6.0) / 6.0
    else:
        lower = ATLAS_HORIZON_TINT_RGBA
        upper = ATLAS_DAY_TINT_RGBA
        weight = min(6.0, sun_alt) / 6.0
    return tuple(
        int(round(start * (1.0 - weight) + end * weight))
        for start, end in zip(lower, upper)
    )


def dimalt_ring_brightness_score_from_rgba(
    red: int,
    green: int,
    blue: int,
    alpha: int,
) -> float:
    """Return an 8-bit brightness score weighted by alpha."""
    luma = ((77 * int(red)) + (150 * int(green)) + (29 * int(blue))) >> 8
    return float(luma) * float(int(alpha)) / 255.0


def dimalt_ring_pen_color_from_color(color: QColor) -> QColor:
    """Return a dimalt stroke color derived from a sampled QColor."""
    score = dimalt_ring_brightness_score_from_rgba(
        color.red(),
        color.green(),
        color.blue(),
        color.alpha(),
    )
    if score > ALT_RING_DIMALT_BRIGHTNESS_THRESHOLD:
        factor = max(0.0, min(1.0, ALT_RING_DIMALT_DARKEN_FACTOR))
        red = int(round(float(color.red()) * factor))
        green = int(round(float(color.green()) * factor))
        blue = int(round(float(color.blue()) * factor))
    else:
        mix = max(0.0, min(1.0, ALT_RING_DIMALT_LIGHTEN_MIX))
        red = int(round(float(color.red()) + (255.0 - float(color.red())) * mix))
        green = int(round(float(color.green()) + (255.0 - float(color.green())) * mix))
        blue = int(round(float(color.blue()) + (255.0 - float(color.blue())) * mix))
    return QColor(red, green, blue, color.alpha())


def format_overlay_info_lines(
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    vmag_limit: float,
    *,
    include_vmag_limit: bool = False,
) -> list[str]:
    """Return the static observation overlay text lines."""

    alt_deg, az_deg = viewer_data.view_center

    def az_to_compass(az: float) -> str:
        names = tuple(DIRECTIONS.keys())
        sector = 360.0 / len(names)
        idx = int(((az % 360.0) + sector / 2.0) // sector) % len(names)
        return names[idx]

    compass = az_to_compass(az_deg)
    deg = "\N{DEGREE SIGN}"

    utc_time = celestial_data.time
    tz_name = viewer_data.timezone_name
    lines: list[str] = []
    try:
        local_tz = ZoneInfo(tz_name)
        local_dt = utc_time.to_datetime(timezone=local_tz)
        lines.append(local_dt.strftime("%Y-%m-%d %H:%M:%S %Z"))
    except Exception:
        lines.append(utc_time.to_datetime().strftime("%Y-%m-%d %H:%M:%S UTC"))

    lines = build_location_info_lines(
        viewer_data.city_name,
        viewer_data.lat_deg,
        viewer_data.lon_deg,
        ground_elevation_m=viewer_data.ground_elevation_m,
        location_height_m=viewer_data.location_height_m,
        height_add_m=viewer_data.height_add_m,
    ) + lines
    lines.append(f"Alt {alt_deg:.0f}{deg}  Az {az_deg:.0f}{deg} ({compass})")
    if include_vmag_limit:
        lines.append(f"Vmag limit {vmag_limit:.1f}")
    return lines


def draw_radial_background(
    painter: QPainter,
    rect: QRectF,
    geometry: ScreenGeometry,
    *,
    theme: ThemeStyle,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float = BACKGROUND_FIELD_OF_VIEW_DEG2,
    opaque: bool = False,
    altaz_rings_mode: str = "dimalt",
    view_center: tuple[float, float] = (0.0, 0.0),
    terrain_profile_altaz: list[tuple[float, float]] | None = None,
) -> None:
    """Draw a radial sky gradient background."""
    assert geometry.radius >= 10
    fov_outer = max(float(BACKGROUND_FIELD_OF_VIEW_DEG1), float(content_fov_deg))
    r_content = float(geometry.radius * (fov_outer / max(1.0e-6, float(edge_fov_deg))))
    cx = float(geometry.center[0])
    cy = float(geometry.center[1])
    corners = (
        (float(rect.left()), float(rect.top())),
        (float(rect.right()), float(rect.top())),
        (float(rect.left()), float(rect.bottom())),
        (float(rect.right()), float(rect.bottom())),
    )
    r_window = max(math.hypot(x - cx, y - cy) for x, y in corners)
    r_max = float(max(r_content + 1.0, r_window))

    def pos(r: float) -> float:
        return max(0.0, min(1.0, r / r_max))

    bg = theme.window_background
    if opaque:
        bg_alpha = 255
    else:
        bg_alpha = None

    if bg.flat_background:
        fill_color = QColor(*bg.inner_rgba)
        if bg_alpha is not None:
            fill_color.setAlpha(bg_alpha)
        painter.save()
        painter.fillRect(rect, fill_color)
        painter.restore()
        return

    def col(r: float, s: float) -> QColor:
        t = max(0.0, min(1.0, r / max(1.0, r_max)))
        rr = int(bg.base_rgb[0] - bg.delta_rgb[0] * t)
        gg = int(bg.base_rgb[1] - bg.delta_rgb[1] * t)
        bb = int(bg.base_rgb[2] - bg.delta_rgb[2] * t)
        aa = int(bg.outer_alpha * (1.0 - s) + bg.edge_alpha * s)
        if bg_alpha is not None:
            aa = bg_alpha
        return QColor(rr, gg, bb, aa)

    c = geometry.center
    g = QRadialGradient(QPointF(c[0], c[1]), r_max)
    inner_color = QColor(*bg.inner_rgba)
    if bg_alpha is not None:
        inner_color.setAlpha(bg_alpha)
    boundary_color = col(r_content, 0.3)
    edge_color = col(r_max, 1.0)
    g.setColorAt(0.0, inner_color)
    g.setColorAt(pos(r_content), inner_color)
    g.setColorAt(pos(r_content + 1.0), boundary_color)
    g.setColorAt(1.0, edge_color)

    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(g)
    painter.drawRect(rect)
    if altaz_rings_mode == "dimalt":
        background_sample_color = sample_background_disc_edge_color(
            rect,
            geometry,
            theme=theme,
            edge_fov_deg=edge_fov_deg,
            content_fov_deg=content_fov_deg,
            opaque=opaque,
        )
        draw_altitude_ring_overlay(
            painter,
            rect,
            geometry,
            view_center=view_center,
            theme=theme,
            edge_fov_deg=edge_fov_deg,
            content_fov_deg=content_fov_deg,
            ring_color=dimalt_ring_pen_color_from_color(background_sample_color),
        )
    painter.restore()


def draw_altitude_ring_overlay(
    painter: QPainter,
    rect: QRectF,
    geometry: ScreenGeometry,
    *,
    view_center: tuple[float, float],
    theme: ThemeStyle,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float = BACKGROUND_FIELD_OF_VIEW_DEG2,
    ring_color: QColor | None = None,
    ring_color_for_alt_deg: Callable[[float], QColor | None] | None = None,
) -> None:
    """Overlay subtle altitude rings clipped to the sky disc."""
    _apply_background_altitude_ring_highlights(
        painter,
        rect,
        geometry,
        view_center=view_center,
        theme=theme,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
        ring_color=ring_color,
        ring_color_for_alt_deg=ring_color_for_alt_deg,
    )


def sample_background_disc_edge_color(
    rect: QRectF,
    geometry: ScreenGeometry,
    *,
    theme: ThemeStyle,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float = BACKGROUND_FIELD_OF_VIEW_DEG2,
    opaque: bool = False,
) -> QColor:
    """Return the background color used nearest the sky-disc boundary."""
    bg = theme.window_background
    if opaque:
        bg_alpha = 255
    else:
        bg_alpha = None

    if bg.flat_background:
        fill_color = QColor(*bg.inner_rgba)
        if bg_alpha is not None:
            fill_color.setAlpha(bg_alpha)
        return fill_color

    fov_outer = max(float(BACKGROUND_FIELD_OF_VIEW_DEG1), float(content_fov_deg))
    r_content = float(geometry.radius * (fov_outer / max(1.0e-6, float(edge_fov_deg))))
    cx = float(geometry.center[0])
    cy = float(geometry.center[1])
    corners = (
        (float(rect.left()), float(rect.top())),
        (float(rect.right()), float(rect.top())),
        (float(rect.left()), float(rect.bottom())),
        (float(rect.right()), float(rect.bottom())),
    )
    r_window = max(math.hypot(x - cx, y - cy) for x, y in corners)
    r_max = float(max(r_content + 1.0, r_window))
    t = max(0.0, min(1.0, r_content / max(1.0, r_max)))
    rr = int(bg.base_rgb[0] - bg.delta_rgb[0] * t)
    gg = int(bg.base_rgb[1] - bg.delta_rgb[1] * t)
    bb = int(bg.base_rgb[2] - bg.delta_rgb[2] * t)
    aa = int(bg.outer_alpha * 0.7 + bg.edge_alpha * 0.3)
    if bg_alpha is not None:
        aa = bg_alpha
    return QColor(rr, gg, bb, aa)


def _apply_background_altitude_ring_highlights(
    painter: QPainter,
    rect: QRectF,
    geometry: ScreenGeometry,
    *,
    view_center: tuple[float, float],
    theme: ThemeStyle,
    edge_fov_deg: float,
    content_fov_deg: float,
    ring_color: QColor | None,
    ring_color_for_alt_deg: Callable[[float], QColor | None] | None,
) -> None:
    """Overlay subtle Alt-ring highlights clipped to the sky disc."""
    if geometry.radius < 1:
        return

    cx = float(geometry.center[0])
    cy = float(geometry.center[1])
    radius = float(geometry.radius)
    disc_fov = max(float(BACKGROUND_FIELD_OF_VIEW_DEG1), float(content_fov_deg))
    disc_radius = radius * (disc_fov / max(1.0e-6, float(edge_fov_deg)))
    clip_path = QPainterPath()
    clip_path.addEllipse(QPointF(cx, cy), float(disc_radius), float(disc_radius))

    painter.save()
    try:
        painter.setClipPath(clip_path)
        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)
        for alt_deg in range(30, 90, 30):
            base_color = (
                ring_color_for_alt_deg(float(alt_deg))
                if ring_color_for_alt_deg is not None
                else ring_color
            )
            ring_pen_color = QColor(255, 255, 255) if base_color is None else QColor(base_color)
            ring_pen = QPen(ring_pen_color, ALT_RING_PEN_WIDTH)
            ring_pen.setCosmetic(True)
            ring_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            ring_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(ring_pen)
            points = QPolygonF()
            az_deg = 0.0
            while az_deg <= 360.0 + 1.0e-6:
                nx, ny = altaz_to_normalized_xy(
                    float(alt_deg),
                    float(az_deg),
                    view_center,
                    edge_fov_deg=edge_fov_deg,
                )
                points.append(QPointF(cx + (nx * radius), cy + (ny * radius)))
                az_deg += 1.0
            painter.drawPolyline(points)
    finally:
        painter.restore()


def draw_window_border(
    painter: QPainter,
    rect: QRectF,
    *,
    theme: ThemeStyle,
    draw_menu_button: bool = True,
) -> None:
    """Draw the menu and resize affordances for the custom window chrome."""
    border_width = FRAMELESS_WINDOW_BORDER_WIDTH
    max_border_width = 0.25 * min(float(rect.width()), float(rect.height()))
    border_width = min(border_width, max_border_width)
    if border_width <= 0.0:
        return

    chrome = theme.window_chrome
    chrome_fill_color = QColor(*chrome.menu_fill_rgba)
    chrome_line_color = QColor(*chrome.menu_icon_rgba)

    top = float(rect.top())
    right = float(rect.right())
    menu_size = float(GUI_BUTTON_SIZE)
    menu_left_edge = right - menu_size + 1.0
    menu_top_edge = top

    painter.save()
    border_width = float(GUI_BUTTON_SIZE)
    border_color = QColor(*theme.window_background.border_rgba)
    border_w = float(border_width)
    left = float(rect.left())
    top = float(rect.top())
    width = float(rect.width())
    height = float(rect.height())
    inner_width = max(0.0, width - (2.0 * border_w))
    painter.fillRect(QRectF(left + border_w, top, inner_width, border_w), border_color)
    painter.fillRect(
        QRectF(left + border_w, top + height - border_w, inner_width, border_w),
        border_color,
    )
    painter.fillRect(QRectF(left, top, border_w, height), border_color)
    painter.fillRect(QRectF(left + width - border_w, top, border_w, height), border_color)
    if draw_menu_button:
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(chrome_fill_color)
        painter.drawRect(QRectF(menu_left_edge, menu_top_edge, menu_size, menu_size))
        menu_icon_pen = QPen(chrome_line_color, 2.0)
        menu_icon_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(menu_icon_pen)
        menu_left = menu_left_edge + 9.0
        menu_right = menu_left_edge + menu_size - 9.0
        for y in (menu_top_edge + 9.0, menu_top_edge + 14.0, menu_top_edge + 19.0):
            painter.drawLine(QPointF(menu_left, y), QPointF(menu_right, y))
    painter.restore()
