from __future__ import annotations

from datetime import timezone

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPen

from ..paths import ThemeStyle
from ..tropical_cyclones.models import TropicalCycloneSnapshot


def _project_lon_lat_to_rect(
    lon_deg: float,
    lat_deg: float,
    rect: QRectF,
) -> QPointF:
    lon = (float(lon_deg) + 180.0) / 360.0
    lat = (90.0 - float(lat_deg)) / 180.0
    x = rect.left() + rect.width() * max(0.0, min(1.0, lon))
    y = rect.top() + rect.height() * max(0.0, min(1.0, lat))
    return QPointF(x, y)


def _fmt_point(snapshot: TropicalCycloneSnapshot) -> str:
    point = snapshot.observed_position
    return f"{point.lat_deg:.1f}, {point.lon_deg:.1f}"


def draw_tropical_cyclone_overlay(
    painter: QPainter,
    viewport_rect: QRectF,
    snapshot: TropicalCycloneSnapshot | None,
    *,
    theme: ThemeStyle,
    text_font: QFont,
    enabled: bool = True,
) -> None:
    if not enabled or snapshot is None:
        return

    margin = 14.0
    card_w = min(320.0, max(230.0, viewport_rect.width() * 0.38))
    card_h = min(190.0, max(140.0, viewport_rect.height() * 0.26))
    card = QRectF(
        viewport_rect.right() - margin - card_w,
        viewport_rect.bottom() - margin - card_h,
        card_w,
        card_h,
    )

    painter.save()
    painter.setRenderHint(QPainter.Antialiasing, True)

    bg = theme.window_chrome.menu_fill_rgba
    border = theme.window_background.border_rgba
    painter.setPen(QPen(QColor(*border), 1.2))
    painter.setBrush(QColor(bg[0], bg[1], bg[2], min(220, max(90, bg[3]))))
    painter.drawRoundedRect(card, 10.0, 10.0)

    padding = 10.0
    title_rect = QRectF(
        card.left() + padding,
        card.top() + padding,
        card.width() - padding * 2.0,
        22.0,
    )
    font = QFont(text_font)
    font.setPointSizeF(max(8.0, float(text_font.pointSizeF()) * 0.92))
    painter.setFont(font)
    painter.setPen(QColor(*theme.text.foreground_rgb[:3]))
    title = snapshot.storm_name
    if snapshot.basin:
        title = f"{title} ({snapshot.basin})"
    painter.drawText(title_rect, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, title)

    subtitle_font = QFont(text_font)
    subtitle_font.setPointSizeF(max(7.0, float(text_font.pointSizeF()) * 0.78))
    painter.setFont(subtitle_font)
    painter.setPen(QColor(*theme.status_text.foreground_rgb[:3]))
    subtitle_y = title_rect.bottom() + 2.0
    if snapshot.advdate_utc is not None:
        advdate_text = snapshot.advdate_utc.astimezone(timezone.utc).strftime("%Y-%m-%d %H:%MZ")
    else:
        advdate_text = "?"
    painter.drawText(
        QRectF(card.left() + padding, subtitle_y, card.width() - padding * 2.0, 18.0),
        Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
        f"ADVDATE {advdate_text}  |  {snapshot.observed_position.lat_deg:.1f}, {snapshot.observed_position.lon_deg:.1f}",
    )

    map_top = subtitle_y + 20.0
    map_rect = QRectF(
        card.left() + padding,
        map_top,
        card.width() - padding * 2.0,
        card.bottom() - padding - map_top,
    )
    painter.setPen(QPen(QColor(255, 255, 255, 32), 1.0))
    painter.setBrush(QColor(0, 0, 0, 26))
    painter.drawRoundedRect(map_rect, 8.0, 8.0)

    grid_pen = QPen(QColor(255, 255, 255, 28), 1.0)
    painter.setPen(grid_pen)
    for lon in (-180, -120, -60, 0, 60, 120, 180):
        p1 = _project_lon_lat_to_rect(lon, -90, map_rect)
        p2 = _project_lon_lat_to_rect(lon, 90, map_rect)
        painter.drawLine(p1, p2)
    for lat in (-60, -30, 0, 30, 60):
        p1 = _project_lon_lat_to_rect(-180, lat, map_rect)
        p2 = _project_lon_lat_to_rect(180, lat, map_rect)
        painter.drawLine(p1, p2)

    forecast = list(snapshot.forecast_positions)
    if forecast:
        line_pen = QPen(QColor(255, 170, 70, 210), 2.0)
        painter.setPen(line_pen)
        prev_point = _project_lon_lat_to_rect(
            snapshot.observed_position.lon_deg,
            snapshot.observed_position.lat_deg,
            map_rect,
        )
        for point in forecast:
            next_point = _project_lon_lat_to_rect(point.lon_deg, point.lat_deg, map_rect)
            painter.drawLine(prev_point, next_point)
            prev_point = next_point

    current_point = _project_lon_lat_to_rect(
        snapshot.observed_position.lon_deg,
        snapshot.observed_position.lat_deg,
        map_rect,
    )
    painter.setBrush(QColor(255, 70, 70, 230))
    painter.setPen(QPen(QColor(255, 210, 210, 240), 1.5))
    painter.drawEllipse(current_point, 4.5, 4.5)

    painter.setFont(subtitle_font)
    forecast_label_pen = QPen(QColor(255, 235, 200, 210), 1.0)
    painter.setPen(forecast_label_pen)
    for idx, point in enumerate(forecast[:5], start=1):
        marker = _project_lon_lat_to_rect(point.lon_deg, point.lat_deg, map_rect)
        painter.setBrush(QColor(255, 170, 70, 220))
        painter.drawEllipse(marker, 3.0, 3.0)
        label = point.label or (f"+{point.tau_hr}h" if point.tau_hr is not None else str(idx))
        painter.drawText(
            QRectF(marker.x() + 5.0, marker.y() - 8.0, 60.0, 16.0),
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
            label,
        )

    painter.restore()
