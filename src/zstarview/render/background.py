import math
from zoneinfo import ZoneInfo

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF, QRadialGradient

from ..paths import (
    BACKGROUND_FIELD_OF_VIEW_DEG1,
    BACKGROUND_FIELD_OF_VIEW_DEG2,
    DIRECTIONS,
    GUI_BUTTON_SIZE,
    GUI_MENU_TEXT_COLOR,
)
from ..types import CelestialData, ScreenGeometry, ViewerData


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
    content_fov_deg: float = BACKGROUND_FIELD_OF_VIEW_DEG2,
) -> None:
    """Draw a radial sky gradient background."""
    assert geometry.radius >= 10
    fov_outer = max(float(BACKGROUND_FIELD_OF_VIEW_DEG1), float(content_fov_deg))
    r_content = float(geometry.radius * (fov_outer / 90.0))
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

    col_params = {
        "white": (242, 46, 245, 48, 250, 50, 255, 200),
        "black": (12, 9, 12, 9, 12, 9, 255, 200),
        "day": (230, 28, 242, 34, 255, 34, 200, 60),
        "night": (10, 7, 12, 9, 16, 11, 200, 60),
    }
    param = col_params.get(preset, None) or col_params["black"]

    def col(r: float, s: float) -> QColor:
        t = max(0.0, min(1.0, r / max(1.0, r_max)))
        rr = int(param[0] - param[1] * t)
        gg = int(param[2] - param[3] * t)
        bb = int(param[4] - param[5] * t)
        aa = int(param[6] * (1.0 - s) + param[7] * s)
        return QColor(rr, gg, bb, aa)

    c = geometry.center
    g = QRadialGradient(QPointF(c[0], c[1]), r_max)
    inner_color = QColor(4, 4, 4, 255) if preset in ("white", "day", "night", "black") else col(0.0, 0.0)
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
    """Draw a broad but subtle border around the window edges."""
    border_width = 16.0
    max_border_width = 0.25 * min(float(rect.width()), float(rect.height()))
    border_width = min(border_width, max_border_width)
    if border_width <= 0.0:
        return

    border_colors = {
        "white": (254, 254, 255, 112),
        "black": (34, 34, 36, 128),
        "day": (250, 252, 255, 35),
        "night": (30, 34, 40, 45),
    }
    rr, gg, bb, aa = border_colors.get(preset, border_colors["night"])
    border_color = QColor(rr, gg, bb, aa)

    left = float(rect.left())
    top = float(rect.top())
    right = float(rect.right())
    bottom = float(rect.bottom())
    menu_size = float(GUI_BUTTON_SIZE)
    menu_left_edge = right - menu_size + 1.0
    menu_top_edge = top

    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(border_color)
    painter.drawRect(QRectF(left, top, border_width, rect.height()))
    painter.drawRect(
        QRectF(
            right - border_width,
            top + menu_size,
            border_width,
            max(0.0, rect.height() - menu_size),
        )
    )
    painter.drawRect(
        QRectF(
            left + border_width,
            top,
            max(0.0, menu_left_edge - (left + border_width)),
            border_width,
        )
    )
    painter.drawRect(
        QRectF(
            left + border_width,
            bottom - border_width,
            max(0.0, rect.width() - 2.0 * border_width),
            border_width,
        )
    )
    painter.drawRect(QRectF(menu_left_edge, menu_top_edge, menu_size, menu_size))
    menu_icon_color = QColor(*GUI_MENU_TEXT_COLOR)
    menu_icon_pen = QPen(menu_icon_color, 2.0)
    menu_icon_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    painter.setPen(menu_icon_pen)
    menu_left = menu_left_edge + 9.0
    menu_right = menu_left_edge + menu_size - 9.0
    for y in (menu_top_edge + 9.0, menu_top_edge + 14.0, menu_top_edge + 19.0):
        painter.drawLine(QPointF(menu_left, y), QPointF(menu_right, y))
    painter.setPen(Qt.PenStyle.NoPen)
    inner_right = right - border_width
    grip_size = 14.0
    inner_bottom = bottom - border_width
    painter.drawPolygon(
        QPolygonF(
            [
                QPointF(inner_right, inner_bottom),
                QPointF(inner_right - grip_size, inner_bottom),
                QPointF(inner_right, inner_bottom - grip_size),
            ]
        )
    )
    painter.restore()
