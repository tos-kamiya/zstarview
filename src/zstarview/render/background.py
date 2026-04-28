import math
from zoneinfo import ZoneInfo

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QRadialGradient

from ..paths import (
    BACKGROUND_FIELD_OF_VIEW_DEG1,
    BACKGROUND_FIELD_OF_VIEW_DEG2,
    DIRECTIONS,
    GUI_BUTTON_SIZE,
    THEME_STYLES_BY_PRESET,
)
from ..types import CelestialData, ScreenGeometry, ViewerData


FRAMELESS_WINDOW_BORDER_WIDTH = 24.0


def format_overlay_info_lines(
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    vmag_limit: float,
    *,
    include_vmag_limit: bool = False,
) -> list[str]:
    """Return the static observation overlay text lines."""

    def format_height_m(value_m: float) -> str:
        rounded = round(float(value_m))
        if abs(float(value_m) - rounded) < 0.05:
            return f"{int(rounded)} m"
        return f"{float(value_m):.1f} m"

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

    location_parts = [viewer_data.city_name]
    if viewer_data.location_height_label and viewer_data.location_height_m is not None:
        location_parts.append(f"{viewer_data.location_height_label} {format_height_m(viewer_data.location_height_m)}")
    if viewer_data.show_observer_height:
        location_parts.append(f"Observer height {format_height_m(viewer_data.observer_height_m)}")
    lines = location_parts + lines
    lines.append(f"Alt {alt_deg:.0f}{deg}  Az {az_deg:.0f}{deg} ({compass})")
    if include_vmag_limit:
        lines.append(f"Vmag limit {vmag_limit:.1f}")
    return lines


def draw_radial_background(
    painter: QPainter,
    rect: QRectF,
    geometry: ScreenGeometry,
    *,
    preset: str = "night",
    edge_fov_deg: float = 90.0,
    content_fov_deg: float = BACKGROUND_FIELD_OF_VIEW_DEG2,
    opaque: bool = False,
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

    theme = THEME_STYLES_BY_PRESET.get(preset, THEME_STYLES_BY_PRESET["black"])
    bg = theme.window_background
    if opaque:
        bg_alpha = 255
    else:
        bg_alpha = None

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
    painter.restore()


def draw_window_border(
    painter: QPainter,
    rect: QRectF,
    *,
    preset: str = "night",
) -> None:
    """Draw the menu and resize affordances for the custom window chrome."""
    border_width = FRAMELESS_WINDOW_BORDER_WIDTH
    max_border_width = 0.25 * min(float(rect.width()), float(rect.height()))
    border_width = min(border_width, max_border_width)
    if border_width <= 0.0:
        return

    theme = THEME_STYLES_BY_PRESET.get(preset, THEME_STYLES_BY_PRESET["night"])
    chrome = theme.window_chrome
    chrome_fill_color = QColor(*chrome.menu_fill_rgba)
    chrome_line_color = QColor(*chrome.menu_icon_rgba)

    top = float(rect.top())
    right = float(rect.right())
    bottom = float(rect.bottom())
    menu_size = float(GUI_BUTTON_SIZE)
    menu_left_edge = right - menu_size + 1.0
    menu_top_edge = top

    painter.save()
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
    painter.setBrush(Qt.BrushStyle.NoBrush)
    grip_inset = 0.0
    inner_right = right - grip_inset
    grip_size = 33.6
    inner_bottom = bottom - grip_inset
    grip_line_inset = 6.0
    grip_pen = QPen(chrome_line_color, 2.0)
    grip_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    painter.setPen(grip_pen)
    painter.drawLine(
        QPointF(inner_right - grip_line_inset, inner_bottom - grip_size + grip_line_inset),
        QPointF(inner_right - grip_size + grip_line_inset, inner_bottom - grip_line_inset),
    )
    if preset in {"white", "day"}:
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
    painter.restore()
