import math
from typing import Any

import numpy as np
from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPainterPath, QPen, QRadialGradient

from ..astro import altaz_to_normalized_xy, calculate_moon_render_data, is_in_fov
from ..paths import ThemeStyle
from ..types import (
    CelestialData,
    CelestialObject,
    PlanetBody,
    ScreenGeometry,
    ViewerData,
)
from ..utils.image import generate_moon_phase_rgba
from .geometry import normalized_to_screen_xy
from .guides import draw_gauge_cross
from .photometry import (
    body_label_text,
    planet_bloom_profile_from_vmag,
    planet_disc_style_from_vmag,
    planet_marker_color,
)
from .qt_image import np_rgba_to_qimage
from .text import (
    LABEL_COLOR_WHITE_BLEND_AMOUNT,
    _rect_overlap_count,
    _text_bounds_at_baseline,
    blend_color_toward_white,
    draw_outlined_text,
    recolor_text_style,
    resolve_label_text_style,
)

_NIGHT_ANNOTATION_RGB = (180, 180, 180)
_ATLAS_PLANET_OUTLINE_RGBA = (24, 24, 24, 150)
_ATLAS_PLANET_OUTLINE_MARGIN_PX = 1.0
_SCENIC_PLANET_OUTLINE_RGBA = (0, 0, 0, 115)
_SCENIC_PLANET_OUTLINE_MARGIN_PX = 1.0


def _content_fov_deg_from_viewer(viewer_data: ViewerData) -> float:
    return float(viewer_data.content_fov_deg)


def _solar_system_annotation_rgb(
    theme: ThemeStyle,
    *,
    instrument_presentation: bool = False,
) -> tuple[int, int, int]:
    if instrument_presentation:
        return tuple(int(value) for value in theme.text.foreground_rgb[:3])
    base_rgb = theme.window_background.base_rgb
    luminance = (
        (77 * int(base_rgb[0])) + (150 * int(base_rgb[1])) + (29 * int(base_rgb[2]))
    ) / 256.0
    if luminance >= 128.0:
        return _NIGHT_ANNOTATION_RGB
    return theme.text.foreground_rgb


def _solar_system_label_style(
    theme: ThemeStyle,
    label_font: QFont,
    label_rgb: tuple[int, int, int],
    *,
    instrument_presentation: bool = False,
) -> Any:
    style = resolve_label_text_style(theme, label_font)
    label_color = (
        QColor(*theme.text.foreground_rgb[:3])
        if instrument_presentation
        else blend_color_toward_white(
            QColor(*label_rgb),
            amount=LABEL_COLOR_WHITE_BLEND_AMOUNT,
        )
    )
    return recolor_text_style(
        style,
        (label_color.red(), label_color.green(), label_color.blue()),
    )


def draw_moon(
    painter: QPainter,
    center: QPointF,
    radius_px: float,
    sun_dir_in_moon_frame: np.ndarray,
    screen_rotation_deg: float,
    opacity: float = 1.0,
    base_color: QColor | None = None,
) -> None:
    img_size = max(5, math.ceil(radius_px * 2.0))
    view_dir = np.array([0, 0, 1], dtype=float)
    if base_color is not None:
        tint_rgba = (base_color.red(), base_color.green(), base_color.blue(), base_color.alpha())
    else:
        tint_rgba = None
    moon_rgba = generate_moon_phase_rgba(img_size, sun_dir_in_moon_frame, view_dir, tint_color=tint_rgba)
    moon_image = np_rgba_to_qimage(moon_rgba)
    source_rect = QRectF(0, 0, img_size, img_size)

    painter.save()
    painter.setOpacity(opacity)
    painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    painter.translate(center)
    if abs(screen_rotation_deg) > 0.1:
        painter.rotate(screen_rotation_deg)
    target_rect = QRectF(-img_size / 2, -img_size / 2, img_size, img_size)
    painter.drawImage(target_rect, moon_image, source_rect)

    if base_color is not None:
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceAtop)
        painter.setBrush(base_color)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawEllipse(QPointF(0.0, 0.0), radius_px, radius_px)

    painter.restore()


def _collect_sun_moon_context(planets: list[PlanetBody]) -> tuple[PlanetBody | None, tuple[float, float] | None, tuple[float, float] | None]:
    moon_body: PlanetBody | None = None
    sun_altaz: tuple[float, float] | None = None
    moon_altaz: tuple[float, float] | None = None
    for body in planets:
        if body.name == "sun":
            sun_altaz = (body.alt, body.az)
        elif body.name == "moon":
            moon_altaz = (body.alt, body.az)
            moon_body = body
    return moon_body, sun_altaz, moon_altaz


def _moon_eclipse_overlay_color(body: PlanetBody) -> QColor | None:
    eclipse = body.lunar_eclipse_info
    if not eclipse or not eclipse.is_eclipse:
        return None
    if eclipse.eclipse_type == "partial":
        return QColor(30, 0, 0, 60)
    if eclipse.eclipse_type == "total":
        return QColor(40, 10, 10, 180)
    return None


def draw_moon_outline(
    painter: QPainter,
    center: QPointF,
    radius_px: float,
    color: QColor,
) -> None:
    outline_radius = max(1.5, float(radius_px))
    pen_width = max(1.25, min(3.0, outline_radius * 0.08))
    pen = QPen(color, pen_width)
    pen.setCosmetic(True)
    painter.save()
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.setPen(pen)
    painter.drawEllipse(center, outline_radius, outline_radius)
    painter.restore()


def draw_moon_phase_outline(
    painter: QPainter,
    center: QPointF,
    radius_px: float,
    sun_dir_in_moon_frame: np.ndarray,
    screen_rotation_deg: float,
    color: QColor,
) -> None:
    """Draw a compact, outline-only representation of the current moon phase."""
    sun_dir = np.asarray(sun_dir_in_moon_frame, dtype=float)
    norm = float(np.linalg.norm(sun_dir))
    if norm <= 1e-9:
        draw_moon_outline(painter, center, radius_px, color)
        return
    sun_dir /= norm

    outline_radius = max(1.5, float(radius_px))
    pen_width = max(1.25, min(3.0, outline_radius * 0.08))
    pen = QPen(color, pen_width)
    pen.setCosmetic(True)

    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setPen(pen)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.translate(center)
    if abs(screen_rotation_deg) > 0.1:
        painter.rotate(screen_rotation_deg)
    projected_sun = sun_dir[:2]
    projected_sun_norm = float(np.linalg.norm(projected_sun))
    if projected_sun_norm <= 1e-6:
        # A front-lit Moon has no visible dark-side limb.  A back-lit Moon
        # has no illuminated limb to draw; the cross still identifies it.
        if sun_dir[2] > 0.0:
            painter.drawEllipse(QPointF(0.0, 0.0), outline_radius, outline_radius)
    else:
        outer_path = QPainterPath()
        drawing = False
        for index in range(193):
            theta = 2.0 * math.pi * index / 192.0
            limb_x = math.cos(theta)
            limb_y = math.sin(theta)
            illuminated = float(limb_x * projected_sun[0] + limb_y * projected_sun[1]) > 0.0
            if not illuminated:
                drawing = False
                continue
            point = QPointF(limb_x * outline_radius, limb_y * outline_radius)
            if not drawing:
                outer_path.moveTo(point)
                drawing = True
            else:
                outer_path.lineTo(point)
        painter.drawPath(outer_path)

    # The terminator is a great circle perpendicular to the Sun direction.
    # Project only its front half; this produces the inner arc of the compact
    # crescent/gibbous marker while retaining the correct waxing/waning side.
    view_dir = np.array([0.0, 0.0, 1.0])
    projected_tangent = np.cross(view_dir, sun_dir)
    tangent_norm = float(np.linalg.norm(projected_tangent))
    if tangent_norm > 1e-6 and abs(float(sun_dir[2])) < 0.9999:
        tangent = projected_tangent / tangent_norm
        other = np.cross(sun_dir, tangent)
        path = QPainterPath()
        drawing = False
        for index in range(97):
            theta = math.pi * index / 96.0
            point_3d = tangent * math.cos(theta) + other * math.sin(theta)
            if point_3d[2] <= 1e-5:
                drawing = False
                continue
            point = QPointF(
                float(point_3d[0] * outline_radius),
                float(point_3d[1] * outline_radius),
            )
            if not drawing:
                path.moveTo(point)
                drawing = True
            else:
                path.lineTo(point)
        painter.drawPath(path)
    painter.restore()


def _draw_moon_planet(
    painter: QPainter,
    pos: QPointF,
    geometry: ScreenGeometry,
    body: PlanetBody,
    viewer_data: ViewerData,
    sun_altaz: tuple[float, float],
    moon_altaz: tuple[float, float],
    enlarge_moon: bool,
    outline_bright_bodies: bool,
    cross_color: QColor,
    marker_scale: float,
    instrument_presentation: bool = False,
    draw_cross: bool = True,
) -> None:
    moon_zoom = 5 if enlarge_moon else 1
    marker_scale = max(1.0, float(marker_scale))
    sun_dir_in_moon_frame, screen_rotation_deg = calculate_moon_render_data(
        sun_altaz,
        moon_altaz,
        viewer_data.view_center,
        edge_fov_deg=float(viewer_data.edge_fov_deg),
    )
    base_moon_radius_px = max((0.25 / float(viewer_data.edge_fov_deg)) * geometry.radius, 2.5)
    moon_radius_px = base_moon_radius_px * moon_zoom * marker_scale
    use_outline = outline_bright_bodies and not enlarge_moon and not instrument_presentation
    if enlarge_moon and draw_cross:
        draw_gauge_cross(painter, cross_color, pos, scale=marker_scale, pen_width=marker_scale)
    if use_outline:
        outline_color = _moon_eclipse_overlay_color(body)
        if outline_color is None:
            outline_color = QColor(220, 220, 220, 220)
        draw_moon_phase_outline(
            painter,
            pos,
            moon_radius_px,
            sun_dir_in_moon_frame=sun_dir_in_moon_frame,
            screen_rotation_deg=screen_rotation_deg,
            color=outline_color,
        )
    else:
        draw_moon(
            painter,
            pos,
            moon_radius_px,
            sun_dir_in_moon_frame=sun_dir_in_moon_frame,
            screen_rotation_deg=screen_rotation_deg,
            opacity=1.0 if instrument_presentation or not enlarge_moon else 0.7,
            base_color=_moon_eclipse_overlay_color(body),
        )
    if not enlarge_moon and draw_cross:
        draw_gauge_cross(painter, cross_color, pos, scale=marker_scale, pen_width=marker_scale)


def _marker_intersects_viewport(painter: QPainter, pos: QPointF, radius_px: float) -> bool:
    viewport_getter = getattr(painter, "viewport", None)
    if not callable(viewport_getter):
        return True
    viewport = viewport_getter()
    if viewport.isNull():
        return True
    r = max(1.0, float(radius_px))
    marker_bounds = QRectF(pos.x() - r, pos.y() - r, r * 2.0, r * 2.0)
    return marker_bounds.intersects(QRectF(viewport))


def draw_planet_disc(
    painter: QPainter,
    pos: QPointF,
    color: QColor,
    *,
    radius_px: float,
    alpha: int,
) -> None:
    r = max(1.5, float(radius_px))
    fill = QColor(color)
    fill.setAlpha(int(np.clip(alpha, 1, 255)))
    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(fill)
    painter.drawEllipse(pos, r, r)
    painter.restore()


def draw_planet_outline(
    painter: QPainter,
    pos: QPointF,
    color: QColor,
    *,
    radius_px: float,
) -> None:
    r = max(1.5, float(radius_px))
    pen_width = max(1.25, min(3.0, r * 0.08))
    pen = QPen(QColor(color))
    pen.setWidthF(pen_width)
    pen.setCosmetic(True)
    painter.save()
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.setPen(pen)
    painter.drawEllipse(pos, r, r)
    painter.restore()


def draw_planet_bloom(
    painter: QPainter,
    pos: QPointF,
    color: QColor,
    *,
    core_radius_px: float,
    vmag: float | None,
) -> None:
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
    painter.setCompositionMode(QPainter.CompositionMode_Plus)
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(gradient)
    painter.drawEllipse(pos, bloom_radius, bloom_radius)
    painter.restore()


def draw_solar_system_bodies(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    celestial_data: CelestialData,
    enlarge_moon: bool,
    outline_bright_bodies: bool = False,
    *,
    text_font: QFont | None = None,
    label_candidates: list[dict[str, Any]] | None = None,
    label_reservations: list[QRectF] | None = None,
    draw_markers: bool = True,
    draw_labels: bool = True,
    theme: ThemeStyle,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float | None = None,
    marker_scale: float = 1.0,
    instrument_presentation: bool = False,
    dark_contrast_enabled: bool = False,
    planet_bodies: list[PlanetBody] | None = None,
    suppress_moon_marker: bool = False,
) -> None:
    bodies = celestial_data.planets if planet_bodies is None else planet_bodies
    moon_body, sun_altaz, moon_altaz = _collect_sun_moon_context(bodies)
    effective_fov_deg = _content_fov_deg_from_viewer(viewer_data) if content_fov_deg is None else float(content_fov_deg)
    marker_scale = max(1.0, float(marker_scale))

    annotation_rgb = _solar_system_annotation_rgb(
        theme,
        instrument_presentation=instrument_presentation,
    )
    text_color = QColor(*annotation_rgb)
    if text_font is not None:
        painter.setFont(text_font)
        label_font = text_font
    else:
        label_font = painter.font() if hasattr(painter, "font") else QFont()

    for body in bodies:
        if not is_in_fov(body.alt, body.az, viewer_data.view_center, fov_deg=effective_fov_deg):
            continue

        pos = QPointF(
            *normalized_to_screen_xy(
                *altaz_to_normalized_xy(
                    body.alt,
                    body.az,
                    viewer_data.view_center,
                    edge_fov_deg=float(edge_fov_deg),
                ),
                geometry,
            )
        )
        marker_visible = True
        if body.name == "moon":
            base_moon_radius_px = max((0.25 / float(edge_fov_deg)) * geometry.radius, 2.5)
            moon_zoom = 5 if enlarge_moon else 1
            marker_visible = _marker_intersects_viewport(
                painter,
                pos,
                base_moon_radius_px * moon_zoom * marker_scale,
            )
        else:
            radius_px, _alpha = planet_disc_style_from_vmag(body.vmag)
            scaled_radius_px = radius_px * marker_scale
            bloom_radius, _center_alpha, _mid_alpha = planet_bloom_profile_from_vmag(
                body.vmag,
                scaled_radius_px,
            )
            marker_visible = _marker_intersects_viewport(
                painter,
                pos,
                max(scaled_radius_px, bloom_radius),
            )

        if draw_markers:
            if body.name == "sun":
                if instrument_presentation:
                    draw_gauge_cross(
                        painter,
                        text_color,
                        pos,
                        scale=marker_scale,
                        pen_width=marker_scale,
                    )
                else:
                    draw_gauge_cross(
                        painter,
                        text_color,
                        pos,
                        scale=marker_scale,
                        pen_width=marker_scale,
                    )
            elif (
                body.name == "moon"
                and not suppress_moon_marker
                and moon_body
                and sun_altaz
                and moon_altaz
            ):
                _draw_moon_planet(
                    painter,
                    pos,
                    geometry,
                    body,
                    viewer_data,
                    sun_altaz,
                    moon_altaz,
                    enlarge_moon,
                    outline_bright_bodies,
                    text_color,
                    marker_scale,
                    instrument_presentation,
                )
            else:
                radius_px, alpha = planet_disc_style_from_vmag(body.vmag)
                radius_px = radius_px * marker_scale
                marker_color = planet_marker_color(body.name)
                if (
                    dark_contrast_enabled
                    and not instrument_presentation
                    and body.vmag is not None
                    and math.isfinite(float(body.vmag))
                    and float(body.vmag) <= 4.0
                ):
                    draw_planet_outline(
                        painter,
                        pos,
                        QColor(*_SCENIC_PLANET_OUTLINE_RGBA),
                        radius_px=radius_px + _SCENIC_PLANET_OUTLINE_MARGIN_PX,
                    )
                if outline_bright_bodies and not instrument_presentation:
                    draw_planet_outline(painter, pos, marker_color, radius_px=radius_px)
                    draw_gauge_cross(
                        painter,
                        text_color,
                        pos,
                        scale=0.55 * marker_scale,
                        pen_width=marker_scale,
                    )
                else:
                    if instrument_presentation:
                        draw_planet_outline(
                            painter,
                            pos,
                            QColor(*_ATLAS_PLANET_OUTLINE_RGBA),
                            radius_px=radius_px + _ATLAS_PLANET_OUTLINE_MARGIN_PX,
                        )
                    draw_planet_bloom(
                        painter,
                        pos,
                        marker_color,
                        core_radius_px=radius_px,
                        vmag=body.vmag,
                    )
                    marker_color.setAlpha(alpha)
                    draw_planet_disc(
                        painter,
                        pos,
                        marker_color,
                        radius_px=radius_px,
                        alpha=255 if instrument_presentation else alpha,
                    )
                    draw_gauge_cross(
                        painter,
                        text_color,
                        pos,
                        scale=0.55 * marker_scale,
                        pen_width=marker_scale,
                    )

        if draw_labels and marker_visible:
            label_text = body_label_text(body.name)
            label_pos = QPointF(pos.x() + 12.0, pos.y() - 10.0)
            label_rgb = (
                planet_marker_color("moon").getRgb()[:3]
                if body.name == "sun"
                else planet_marker_color(body.name).getRgb()[:3]
            )
            label_style = _solar_system_label_style(
                theme,
                label_font,
                label_rgb,
                instrument_presentation=instrument_presentation,
            )
            if label_candidates is not None:
                label_candidates.append(
                    {
                        "text": label_text,
                        "pos": label_pos,
                        "style": label_style,
                        "priority": 40,
                        "hide_on_overlap": True,
                    }
                )
                continue
            if label_reservations is not None:
                rect = _text_bounds_at_baseline(label_text, label_font, label_pos)
                if _rect_overlap_count(rect, label_reservations) > 0:
                    continue
                label_reservations.append(rect)
            draw_outlined_text(
                painter,
                label_text,
                label_pos,
                label_font,
                style=label_style,
            )


def draw_hovered_moon_overlay(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    celestial_data: CelestialData,
    highlighted_object: tuple[CelestialObject, QPointF] | None,
    *,
    marker_scale: float = 1.0,
    outline_bright_bodies: bool = False,
    theme: ThemeStyle,
) -> None:
    if highlighted_object is None:
        return
    obj, pos = highlighted_object
    obj_name = getattr(obj, "name", None) if hasattr(obj, "name") else obj.get("name")
    if str(obj_name).strip().lower() != "moon":
        return
    moon_body, sun_altaz, moon_altaz = _collect_sun_moon_context(celestial_data.planets)
    if moon_body is None or sun_altaz is None or moon_altaz is None:
        return
    text_color = QColor(*theme.text.foreground_rgb)
    _draw_moon_planet(
        painter,
        pos,
        geometry,
        moon_body,
        viewer_data,
        sun_altaz,
        moon_altaz,
        True,
        outline_bright_bodies,
        text_color,
        marker_scale,
        draw_cross=True,
    )
