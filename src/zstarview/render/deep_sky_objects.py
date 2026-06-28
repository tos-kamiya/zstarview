import math
import re
from typing import Optional, Tuple

import numpy as np
from PySide6.QtCore import QPoint, QPointF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy
from ..paths import ThemeStyle, _rgb_from_hsv
from ..types import CelestialData, CelestialObject, ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy

_DSO_SHAPE_SIZE_GAIN = 1.0
_DSO_HOVER_SIZE_GAIN = 3.0
_DSO_SHAPE_MIN_MAJOR_ARCMIN = 15.0
_DSO_CATALOG_LIKE_NAME_RE = re.compile(r"^(M\d+|NGC\d+|IC\d+|MEL\d+|MWSC\d+)$", re.IGNORECASE)
DSO_LABEL_RGB = (122, 173, 240)
DSO_LABEL_TEXT_RGB = _rgb_from_hsv(214.1, 49.2, 86.0)


def _is_named_dso(name: object, obj_id: object) -> bool:
    name_s = "" if name is None else str(name).strip()
    if not name_s:
        return False
    if _DSO_CATALOG_LIKE_NAME_RE.match(name_s):
        return False
    id_s = "" if obj_id is None else str(obj_id).strip()
    return name_s.upper().replace(" ", "") != id_s.upper().replace(" ", "")


def _dso_ellipse_polygon(
    *,
    alt_deg: float,
    az_deg: float,
    major_arcmin: float,
    minor_arcmin: float,
    pa_deg: float,
    viewer_data: ViewerData,
    geometry: ScreenGeometry,
    gain: float,
    samples: int = 48,
) -> QPolygonF:
    a_deg = 0.5 * (float(major_arcmin) / 60.0) * float(gain)
    b_deg = 0.5 * (float(minor_arcmin) / 60.0) * float(gain)
    pa_rad = math.radians(float(pa_deg) if math.isfinite(pa_deg) else 0.0)
    c = math.cos(pa_rad)
    s = math.sin(pa_rad)
    cos_alt = max(0.15, math.cos(math.radians(float(alt_deg))))

    pts: list[QPointF] = []
    n = max(16, int(samples))
    for i in range(n):
        t = 2.0 * math.pi * (i / n)
        local_major = a_deg * math.cos(t)
        local_minor = b_deg * math.sin(t)
        d_north = (local_major * c) + (local_minor * -s)
        d_east = (local_major * s) + (local_minor * c)

        alt_i = float(alt_deg) + d_north
        alt_i = max(-89.8, min(89.8, alt_i))
        az_i = (float(az_deg) + (d_east / cos_alt)) % 360.0
        nx, ny = altaz_to_normalized_xy(
            alt_i,
            az_i,
            viewer_data.view_center,
            edge_fov_deg=float(viewer_data.edge_fov_deg),
        )
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        pts.append(QPointF(px, py))
    return QPolygonF(pts)


def find_highlighted_dso(
    celestial_data: Optional[CelestialData],
    viewer_data: ViewerData,
    mouse_pos: QPoint,
    geometry: ScreenGeometry,
) -> Optional[Tuple[CelestialObject, QPointF]]:
    if not celestial_data:
        return None

    dso = celestial_data.deep_sky_objects
    if dso["alt"].size == 0:
        return None

    names = np.asarray(dso["name"], dtype=object)
    ids = np.asarray(dso["id"], dtype=object)
    major = np.asarray(dso["major_arcmin"], dtype=float)
    minor = np.asarray(dso["minor_arcmin"], dtype=float)
    pa = np.asarray(dso["pa_deg"], dtype=float)
    named_mask = np.array([_is_named_dso(n, i) for n, i in zip(names, ids)], dtype=bool)
    shape_mask = np.isfinite(major) & np.isfinite(minor) & (major >= _DSO_SHAPE_MIN_MAJOR_ARCMIN) & (minor > 0.0)
    valid = np.nonzero(named_mask & shape_mask)[0]
    if valid.size == 0:
        return None

    mouse_pf = QPointF(float(mouse_pos.x()), float(mouse_pos.y()))
    best: Optional[Tuple[CelestialObject, QPointF]] = None
    best_dist = float("inf")
    near_best: Optional[Tuple[CelestialObject, QPointF]] = None
    near_best_dist = 30.0**2
    for idx in valid:
        alt = float(dso["alt"][idx])
        az = float(dso["az"][idx])
        nx, ny = altaz_to_normalized_xy(
            alt,
            az,
            viewer_data.view_center,
            edge_fov_deg=float(viewer_data.edge_fov_deg),
        )
        x, y = normalized_to_screen_xy(nx, ny, geometry)
        dist_sq = (mouse_pf.x() - x) ** 2 + (mouse_pf.y() - y) ** 2
        poly = _dso_ellipse_polygon(
            alt_deg=alt,
            az_deg=az,
            major_arcmin=float(major[idx]),
            minor_arcmin=float(minor[idx]),
            pa_deg=float(pa[idx]) if np.isfinite(pa[idx]) else 0.0,
            viewer_data=viewer_data,
            geometry=geometry,
            gain=_DSO_SHAPE_SIZE_GAIN,
            samples=48,
        )
        candidate: CelestialObject = {
            "name": dso["name"][idx],
            "alt": alt,
            "az": az,
            "id": dso["id"][idx],
            "type": dso["type"][idx],
            "major_arcmin": major[idx],
            "minor_arcmin": minor[idx],
            "pa_deg": pa[idx],
            "_kind": "dso",
        }
        if poly.containsPoint(mouse_pf, Qt.FillRule.OddEvenFill) and dist_sq < best_dist:
            best_dist = dist_sq
            best = (candidate, QPointF(x, y))
        if dist_sq < near_best_dist:
            near_best_dist = dist_sq
            near_best = (candidate, QPointF(x, y))

    return best if best is not None else near_best


def draw_deep_sky_shapes(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    celestial_data: CelestialData,
    *,
    opacity_scale: float = 1.0,
) -> None:
    dso = celestial_data.deep_sky_objects
    if dso["alt"].size == 0:
        return

    major = np.asarray(dso["major_arcmin"], dtype=float)
    minor = np.asarray(dso["minor_arcmin"], dtype=float)
    finite_shape = np.isfinite(major) & np.isfinite(minor) & (major >= _DSO_SHAPE_MIN_MAJOR_ARCMIN) & (minor > 0.0)
    if not np.any(finite_shape):
        return

    alpha_scale = max(0.0, min(1.0, float(opacity_scale)))
    if alpha_scale <= 0.0:
        return

    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)

    indices = np.nonzero(finite_shape)[0]
    for idx in indices:
        alt = float(dso["alt"][idx])
        az = float(dso["az"][idx])
        major_deg = float(major[idx]) / 60.0
        alpha = int(
            np.clip(
                round((95.0 - 14.0 * math.sqrt(max(0.0, major_deg))) * alpha_scale),
                0,
                255,
            )
        )
        painter.setBrush(QColor(DSO_LABEL_RGB[0], DSO_LABEL_RGB[1], DSO_LABEL_RGB[2], alpha))
        poly = _dso_ellipse_polygon(
            alt_deg=alt,
            az_deg=az,
            major_arcmin=float(major[idx]),
            minor_arcmin=float(minor[idx]),
            pa_deg=float(dso["pa_deg"][idx]) if np.isfinite(dso["pa_deg"][idx]) else 0.0,
            viewer_data=viewer_data,
            geometry=geometry,
            gain=_DSO_SHAPE_SIZE_GAIN,
        )
        painter.drawPolygon(poly)

    painter.restore()


def draw_dso_hover_info(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    highlighted_dso: Optional[Tuple[CelestialObject, QPointF]],
    text_font: QFont,
    *,
    theme: ThemeStyle,
) -> None:
    del text_font
    if not highlighted_dso:
        return
    obj, _ = highlighted_dso
    major_arcmin = float(obj.get("major_arcmin", 0.0))
    minor_arcmin = float(obj.get("minor_arcmin", 0.0))
    pa_deg = float(obj.get("pa_deg", 0.0))
    alt = float(obj.get("alt", 0.0))
    az = float(obj.get("az", 0.0))
    hover_fill_alpha = 70 if theme.label_outline_suppressed else 62
    hover_pen_alpha = 180 if theme.label_outline_suppressed else 190
    hover_pen = QColor(DSO_LABEL_RGB[0], DSO_LABEL_RGB[1], DSO_LABEL_RGB[2], hover_pen_alpha)
    hover_fill = QColor(DSO_LABEL_RGB[0], DSO_LABEL_RGB[1], DSO_LABEL_RGB[2], hover_fill_alpha)
    base_poly = _dso_ellipse_polygon(
        alt_deg=alt,
        az_deg=az,
        major_arcmin=major_arcmin,
        minor_arcmin=minor_arcmin,
        pa_deg=pa_deg if math.isfinite(pa_deg) else 0.0,
        viewer_data=viewer_data,
        geometry=geometry,
        gain=_DSO_SHAPE_SIZE_GAIN,
        samples=60,
    )
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(hover_fill)
    painter.drawPolygon(base_poly)

    hover_poly = _dso_ellipse_polygon(
        alt_deg=alt,
        az_deg=az,
        major_arcmin=major_arcmin,
        minor_arcmin=minor_arcmin,
        pa_deg=pa_deg if math.isfinite(pa_deg) else 0.0,
        viewer_data=viewer_data,
        geometry=geometry,
        gain=_DSO_HOVER_SIZE_GAIN,
        samples=60,
    )
    painter.setPen(QPen(hover_pen, 2.2))
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawPolygon(hover_poly)
