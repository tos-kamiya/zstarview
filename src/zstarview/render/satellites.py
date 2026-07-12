import astropy.time
from PySide6.QtCore import QPoint, QPointF
from PySide6.QtGui import QColor, QFont, QPainter

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import THEME_STYLES_BY_PRESET, ThemeStyle
from ..satellite_constants import SATELLITE_OVERLAY_MARKER_MAX_ALPHA
from ..satellites import project_satellite_records
from ..satellites.types import SatelliteOverlayPoint
from ..types import ScreenGeometry, ViewerData
from . import text as render_text
from .geometry import normalized_to_screen_xy
from .guides import draw_gauge_cross

_SATELLITE_HOVER_MIN_RADIUS_PX = 12.0
_SATELLITE_HOVER_RADIUS_SCALE = 20.0
_SIMPLIFIED_SATELLITE_LABEL_ALPHA = 0.7


def _satellite_hover_radius_px(point: SatelliteOverlayPoint) -> float:
    return max(
        _SATELLITE_HOVER_MIN_RADIUS_PX,
        float(point.marker_scale) * _SATELLITE_HOVER_RADIUS_SCALE,
    )


def find_highlighted_satellite(
    satellite_records_by_group: object | None = None,
    mouse_pos: QPoint | QPointF | None = None,
    geometry: ScreenGeometry | None = None,
    *,
    viewer_data: ViewerData | None = None,
    time_obj: astropy.time.Time | None = None,
) -> tuple[SatelliteOverlayPoint, QPointF] | None:
    if viewer_data is None or time_obj is None:
        return None
    view_center = viewer_data.view_center
    edge_fov_deg = float(viewer_data.edge_fov_deg)
    content_fov_deg = float(viewer_data.content_fov_deg)
    if geometry is None:
        return None
    satellite_points = project_satellite_records(
        satellite_records_by_group or {},
        observer_lat=float(viewer_data.lat_deg),
        observer_lon=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m),
        time_obj=time_obj,
    )
    if not satellite_points or mouse_pos is None:
        return None

    mouse_x = float(mouse_pos.x())
    mouse_y = float(mouse_pos.y())
    best_point: tuple[SatelliteOverlayPoint, QPointF] | None = None
    best_dist_sq = float("inf")
    for point in satellite_points:
        alt = float(point.alt_deg)
        az = float(point.az_deg)
        if not is_in_fov(alt, az, view_center, fov_deg=content_fov_deg):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center, edge_fov_deg=edge_fov_deg)
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        dist_sq = (mouse_x - float(px)) ** 2 + (mouse_y - float(py)) ** 2
        hover_radius = _satellite_hover_radius_px(point)
        if dist_sq <= hover_radius * hover_radius and dist_sq < best_dist_sq:
            best_dist_sq = dist_sq
            best_point = (point, QPointF(float(px), float(py)))
    return best_point


def draw_satellite_overlay(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer_data: ViewerData | None = None,
    satellite_records_by_group: object | None = None,
    *,
    time_obj: astropy.time.Time | None = None,
    opacity: float = 1.0,
    highlighted_satellite: SatelliteOverlayPoint | None = None,
    marker_scale: float = 1.0,
    draw_simplified_labels: bool = False,
    text_font: QFont | None = None,
    theme: ThemeStyle = THEME_STYLES_BY_PRESET["night"],
) -> None:
    if viewer_data is None or time_obj is None:
        return
    view_center = viewer_data.view_center
    edge_fov_deg = float(viewer_data.edge_fov_deg)
    content_fov_deg = float(viewer_data.content_fov_deg)
    satellite_style = theme.overlays.satellite
    layer_opacity = max(
        0.0,
        min(1.0, float(opacity) * float(satellite_style.alpha_scale)),
    )
    if layer_opacity <= 0.0:
        return
    satellite_points = project_satellite_records(
        satellite_records_by_group or {},
        observer_lat=float(viewer_data.lat_deg),
        observer_lon=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m),
        time_obj=time_obj,
    )
    if not satellite_points:
        return

    painter.save()
    width_scale = max(
        1.0,
        float(marker_scale) * float(satellite_style.width_scale),
    )
    marker_color = QColor(
        *satellite_style.rgb,
        max(
            0, min(255, int(round(SATELLITE_OVERLAY_MARKER_MAX_ALPHA * layer_opacity)))
        ),
    )
    label_font = text_font
    if label_font is None:
        try:
            label_font = painter.font()
        except Exception:
            label_font = QFont()
    label_style = None
    if draw_simplified_labels:
        label_style = render_text.resolve_overlay_label_text_style(
            theme,
            satellite_style,
            label_font,
            opacity=layer_opacity * _SIMPLIFIED_SATELLITE_LABEL_ALPHA,
        )
    for point in satellite_points:
        alt = float(point.alt_deg)
        az = float(point.az_deg)
        if not is_in_fov(alt, az, view_center, fov_deg=content_fov_deg):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center, edge_fov_deg=edge_fov_deg)
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        pos = QPointF(float(px), float(py))
        cross_scale = float(point.marker_scale) * width_scale
        cross_pen_width = float(satellite_style.marker_width) + (
            1.0 if point is highlighted_satellite else 0.0
        )
        if satellite_style.outline_rgba is not None:
            outline_color = QColor(*satellite_style.outline_rgba)
            outline_color.setAlpha(
                max(0, min(255, int(round(outline_color.alpha() * layer_opacity))))
            )
            draw_gauge_cross(
                painter,
                outline_color,
                pos,
                scale=cross_scale,
                pen_width=cross_pen_width + 2.0,
            )
        draw_gauge_cross(
            painter,
            marker_color,
            pos,
            scale=cross_scale,
            pen_width=cross_pen_width,
        )
        if label_style is not None:
            satellite_name = str(point.satellite_name).strip()
            if satellite_name:
                text_bounds = render_text._text_bounds_at_baseline(
                    satellite_name,
                    label_font,
                    QPointF(0.0, 0.0),
                )
                label_pos = QPointF(
                    float(pos.x()) - float(text_bounds.left()),
                    float(pos.y()) - float(text_bounds.bottom()),
                )
                render_text.draw_outlined_text(
                    painter,
                    satellite_name,
                    label_pos,
                    style=label_style,
                )
    painter.restore()


def draw_satellite_highlight_overlay(
    painter: QPainter,
    highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None,
    *,
    opacity: float = 1.0,
    marker_scale: float = 1.0,
    theme: ThemeStyle = THEME_STYLES_BY_PRESET["night"],
) -> None:
    if highlighted_satellite is None:
        return
    point, pos = highlighted_satellite
    satellite_style = theme.overlays.satellite
    layer_opacity = max(
        0.0,
        min(1.0, float(opacity) * float(satellite_style.alpha_scale)),
    )
    if layer_opacity <= 0.0:
        return
    marker_color = QColor(
        *satellite_style.rgb,
        max(
            0, min(255, int(round(SATELLITE_OVERLAY_MARKER_MAX_ALPHA * layer_opacity)))
        ),
    )
    cross_scale = float(point.marker_scale) * max(
        1.0,
        float(marker_scale) * float(satellite_style.width_scale),
    )
    cross_pen_width = float(satellite_style.marker_width) + 1.0
    if satellite_style.outline_rgba is not None:
        outline_color = QColor(*satellite_style.outline_rgba)
        outline_color.setAlpha(
            max(0, min(255, int(round(outline_color.alpha() * layer_opacity))))
        )
        draw_gauge_cross(
            painter,
            outline_color,
            pos,
            scale=cross_scale,
            pen_width=cross_pen_width + 2.0,
        )
    draw_gauge_cross(
        painter,
        marker_color,
        pos,
        scale=cross_scale,
        pen_width=cross_pen_width,
    )
