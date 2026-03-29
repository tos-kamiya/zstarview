from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Any

import numpy as np
from PySide6.QtCore import QPoint, QRect, QRectF, Qt
from PySide6.QtGui import QFont, QImage, QPainter

from ..types import CelestialData, ViewerData
from . import asterisms as render_asterisms
from . import background as render_background
from . import draw as render_draw
from . import deep_sky_objects as render_deep_sky_objects
from . import guides as render_guides
from . import overlay_info as render_overlay_info
from . import aircraft as render_aircraft
from . import satellites as render_satellites
from . import solar_system as render_solar_system
from . import terrain as render_terrain
from . import text as render_text

ORIENTATION_INTERACTION_STAR_VMAG_LIMIT = 4.0


def compute_star_render_surface_size(
    width_px: int,
    height_px: int,
    disc_width_px: int,
    expected_width_px: int,
) -> tuple[int, int]:
    """Return internal star-render surface size."""
    w = max(1, int(width_px))
    h = max(1, int(height_px))
    disc_w = max(1, int(disc_width_px))
    base = max(1, int(expected_width_px))
    if disc_w <= base:
        return (w, h)
    rendered_disc_w = float(base) * math.sqrt(float(disc_w) / float(base))
    scale = rendered_disc_w / float(disc_w)
    return (
        max(1, int(round(w * scale))),
        max(1, int(round(h * scale))),
    )


def compute_star_render_upscale_factor(
    disc_width_px: int,
    expected_width_px: int,
) -> float:
    """Return the on-screen enlargement factor caused by star-layer downsampling."""
    rendered_w, _ = compute_star_render_surface_size(
        disc_width_px,
        disc_width_px,
        disc_width_px,
        expected_width_px,
    )
    return float(disc_width_px) / float(max(1, rendered_w))


@dataclass(frozen=True)
class RenderSceneData:
    viewer: ViewerData
    celestial_data: CelestialData
    sky_disc_image: QImage | None
    cloud_image: np.ndarray | None
    cloud_missing_mask: np.ndarray | None
    cloud_stripe_density: Any | None
    terrain_horizon_profile: list[tuple[float, float]] | None
    urban_outlines: Any | None
    satellite_overlay_points: Any | None
    aircraft_overlay_points: Any | None


@dataclass(frozen=True)
class RenderStyle:
    visual_preset: str
    text_font: QFont
    status_line_font: QFont
    show_background_gradient: bool
    show_overlay_info: bool
    show_dso: bool
    show_asterisms: bool
    show_guidelines: bool
    enlarge_moon: bool
    star_base_radius: float
    star_visibility_boost: float
    vmag_limit: float
    cloud_disc_alpha: float
    satellite_opacity: float
    terrain_horizon_opacity: float
    urban_outline_opacity: float
    show_urban_outline_layer: bool
    aircraft_opacity: float
    star_render_expected_width: int


@dataclass(frozen=True)
class RenderHudState:
    mouse_pos: QPoint | None
    viewport_interaction_mode: bool
    viewport_interaction_stars: Any | None
    status_message: str | None


def _content_fov_deg(scene: RenderSceneData) -> float:
    return float(scene.viewer.content_fov_deg)


def _window_size(viewport_rect: QRect) -> tuple[int, int]:
    return (int(viewport_rect.width()), int(viewport_rect.height()))


def render_base_scene_into_painter(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    style: RenderStyle,
    hud: RenderHudState,
    compositor: Any,
) -> None:
    clear_background_layer(painter, viewport_rect)
    draw_background_layer(
        painter,
        geometry=geometry,
        viewport_rect=viewport_rect,
        scene=scene,
        style=style,
    )
    sky_cloud_style = replace(style, cloud_disc_alpha=0.0) if hud.viewport_interaction_mode else style
    draw_sky_cloud_layers(
        painter,
        geometry=geometry,
        scene=scene,
        style=sky_cloud_style,
        compositor=compositor,
    )
    draw_guide_layer(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
    )
    if hud.viewport_interaction_mode:
        draw_viewport_interaction_layers(
            painter,
            geometry=geometry,
            viewport_rect=viewport_rect,
            scene=scene,
            style=style,
            hud=hud,
        )
        return

    label_reservations: list[QRectF] = []
    label_candidates: list[dict[str, Any]] = []
    draw_terrain_layers(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
        highlighted_object=None,
        label_reservations=label_reservations,
        label_candidates=label_candidates,
    )
    draw_star_layer(
        painter,
        geometry=geometry,
        viewport_rect=viewport_rect,
        scene=scene,
        style=style,
    )
    draw_planet_layer(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
        enlarge_moon=bool(style.enlarge_moon),
        label_candidates=label_candidates,
    )
    draw_satellite_layer(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
        label_candidates=label_candidates,
    )
    draw_aircraft_layer(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
        label_candidates=label_candidates,
    )
    draw_label_layer(
        painter,
        style=style,
        label_candidates=label_candidates,
    )
def render_hud_overlay_into_painter(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    style: RenderStyle,
    hud: RenderHudState,
    highlighted_object: Any | None,
    highlighted_dso: Any | None,
) -> None:
    if hud.viewport_interaction_mode:
        draw_status_line(
            painter,
            viewport_rect=viewport_rect,
            style=style,
            hud=hud,
        )
        return

    draw_hover_overlay_layer(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
        highlighted_object=highlighted_object,
        highlighted_dso=highlighted_dso,
    )
    draw_overlay_layer(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
        mouse_pos=hud.mouse_pos,
        highlighted_object=None,
        highlighted_dso=None,
        enlarge_moon=bool(style.enlarge_moon),
        label_reservations=[],
        label_candidates=[],
    )
    draw_status_line(
        painter,
        viewport_rect=viewport_rect,
        style=style,
        hud=hud,
    )


def draw_viewport_interaction_layers(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    style: RenderStyle,
    hud: RenderHudState,
) -> None:
    line_width_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    interaction_celestial_data = scene.celestial_data
    if hud.viewport_interaction_stars is not None:
        interaction_celestial_data = replace(
            scene.celestial_data,
            stars=hud.viewport_interaction_stars,
        )
    if style.show_guidelines:
        render_guides.draw_sky_reference_lines(
            painter,
            geometry,
            scene.celestial_data,
            scene.viewer,
        )
    draw_star_layer(
        painter,
        geometry=geometry,
        viewport_rect=viewport_rect,
        scene=replace(scene, celestial_data=interaction_celestial_data),
        style=style,
        draw_vmag_limit=ORIENTATION_INTERACTION_STAR_VMAG_LIMIT,
    )
    render_terrain.draw_terrain_horizon_line(
        painter,
        geometry,
        scene.terrain_horizon_profile,
        scene.viewer.view_center,
        opacity=style.terrain_horizon_opacity,
        line_width_scale=line_width_scale,
        content_fov_deg=_content_fov_deg(scene),
    )


def clear_background_layer(painter: QPainter, viewport_rect: QRect) -> None:
    painter.save()
    painter.setCompositionMode(QPainter.CompositionMode_Clear)
    painter.fillRect(viewport_rect, Qt.transparent)
    painter.setCompositionMode(QPainter.CompositionMode_SourceOver)
    painter.restore()


def draw_background_layer(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    style: RenderStyle,
) -> None:
    if not style.show_background_gradient:
        return
    render_background.draw_radial_background(
        painter,
        QRectF(viewport_rect),
        geometry,
        preset=style.visual_preset,
        content_fov_deg=_content_fov_deg(scene),
    )
    render_background.draw_window_frame(
        painter,
        QRectF(viewport_rect),
        preset=style.visual_preset,
    )


def draw_guide_layer(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
) -> None:
    """Draw guide annotations that should float above sky/cloud but below scene overlays."""
    if not style.show_guidelines:
        return
    content_fov_deg = _content_fov_deg(scene)
    render_guides.draw_direction_labels(
        painter,
        geometry,
        scene.viewer.view_center,
        style.text_font,
        None,
        preset=style.visual_preset,
        content_fov_deg=content_fov_deg,
    )
    render_guides.draw_zenith_marker(
        painter,
        geometry,
        scene.viewer.view_center,
        content_fov_deg=content_fov_deg,
    )


def draw_sky_cloud_layers(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    compositor: Any,
) -> None:
    compositor.draw(
        painter,
        geometry,
        scene.sky_disc_image,
        scene.cloud_image,
        cloud_alpha=style.cloud_disc_alpha,
        view_center=scene.viewer.view_center,
        observer_lat_deg=scene.viewer.location[0],
        stripe_density=scene.cloud_stripe_density,
        missing_mask=scene.cloud_missing_mask,
        terrain_profile_altaz=(
            scene.terrain_horizon_profile if style.terrain_horizon_opacity > 0.0 else None
        ),
        content_fov_deg=_content_fov_deg(scene),
    )


def draw_terrain_layers(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    highlighted_object: Any | None,
    label_reservations: list[QRectF],
    label_candidates: list[dict[str, Any]],
) -> None:
    content_fov_deg = _content_fov_deg(scene)
    line_width_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    if style.show_dso:
        render_deep_sky_objects.draw_deep_sky_shapes(
            painter,
            geometry,
            scene.celestial_data,
            scene.viewer,
            preset=style.visual_preset,
        )
    if style.show_asterisms:
        render_asterisms.draw_asterisms(
            painter,
            geometry,
            scene.celestial_data,
            scene.viewer,
            highlighted_object,
            style.text_font,
            label_reservations,
            label_candidates=label_candidates,
            preset=style.visual_preset,
            line_width_scale=line_width_scale,
            content_fov_deg=content_fov_deg,
            draw_base=True,
            draw_highlight=False,
        )
    if style.show_guidelines:
        render_guides.draw_sky_reference_lines(
            painter,
            geometry,
            scene.celestial_data,
            scene.viewer,
        )
    render_terrain.draw_terrain_horizon_line(
        painter,
        geometry,
        scene.terrain_horizon_profile,
        scene.viewer.view_center,
        opacity=style.terrain_horizon_opacity,
        line_width_scale=line_width_scale,
        content_fov_deg=content_fov_deg,
    )
    draw_urban_outline_layer(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
    )


def draw_dso_hover_layer(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    highlighted_dso: Any | None,
) -> None:
    if not style.show_dso:
        return
    render_deep_sky_objects.draw_dso_hover_info(
        painter,
        geometry,
        scene.viewer,
        highlighted_dso,
        style.text_font,
        preset=style.visual_preset,
    )


def draw_urban_outline_layer(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
) -> None:
    if not style.show_urban_outline_layer:
        return
    line_width_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    render_terrain.draw_urban_outlines(
        painter,
        geometry,
        scene.urban_outlines,
        scene.viewer.view_center,
        opacity=style.urban_outline_opacity,
        line_width_scale=line_width_scale,
        content_fov_deg=_content_fov_deg(scene),
    )


def draw_aircraft_layer(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    label_candidates: list[dict[str, Any]],
) -> None:
    line_width_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    render_aircraft.draw_aircraft_overlay(
        painter,
        geometry,
        scene.aircraft_overlay_points,
        scene.viewer.view_center,
        opacity=style.aircraft_opacity,
        line_width_scale=line_width_scale,
        label_candidates=label_candidates,
        preset=style.visual_preset,
        content_fov_deg=_content_fov_deg(scene),
    )


def draw_star_layer(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    style: RenderStyle,
    draw_vmag_limit: float | None = None,
) -> None:
    draw_data = scene.celestial_data
    win_w, win_h = _window_size(viewport_rect)
    low_w, low_h = compute_star_render_surface_size(
        win_w,
        win_h,
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    content_fov_deg = _content_fov_deg(scene)
    if low_w == win_w and low_h == win_h:
        render_draw.draw_stars(
            painter,
            geometry,
            draw_data,
            scene.viewer,
            style.star_base_radius,
            visibility_boost=style.star_visibility_boost,
            draw_vmag_limit=draw_vmag_limit if draw_vmag_limit is not None else style.vmag_limit,
            viewport_size=(win_w, win_h),
            content_fov_deg=content_fov_deg,
        )
        return

    low_img = QImage(low_w, low_h, QImage.Format.Format_ARGB32_Premultiplied)
    low_img.fill(Qt.GlobalColor.transparent)
    low_painter = QPainter(low_img)
    low_painter.setRenderHint(QPainter.Antialiasing, False)
    low_painter.setRenderHint(QPainter.SmoothPixmapTransform, False)
    sx = low_w / max(1.0, float(win_w))
    sy = low_h / max(1.0, float(win_h))
    low_geometry = render_draw.ScreenGeometry(
        center=(
            int(round(geometry.center[0] * sx)),
            int(round(geometry.center[1] * sy)),
        ),
        radius=max(1, int(round(geometry.radius * min(sx, sy)))),
    )
    render_draw.draw_stars(
        low_painter,
        low_geometry,
        draw_data,
        scene.viewer,
        style.star_base_radius,
        visibility_boost=style.star_visibility_boost,
        draw_vmag_limit=draw_vmag_limit if draw_vmag_limit is not None else style.vmag_limit,
        viewport_size=(low_w, low_h),
        content_fov_deg=content_fov_deg,
    )
    low_painter.end()

    painter.save()
    painter.setRenderHint(QPainter.SmoothPixmapTransform, False)
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
    painter.drawImage(viewport_rect, low_img)
    painter.restore()


def draw_planet_layer(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    enlarge_moon: bool,
    label_candidates: list[dict[str, Any]],
) -> None:
    render_solar_system.draw_solar_system_bodies(
        painter,
        geometry,
        scene.celestial_data,
        scene.viewer,
        enlarge_moon,
        text_font=style.text_font,
        label_candidates=label_candidates,
        draw_labels=True,
        preset=style.visual_preset,
        content_fov_deg=_content_fov_deg(scene),
    )


def draw_overlay_layer(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    mouse_pos: QPoint | None,
    highlighted_object: Any | None,
    highlighted_dso: Any | None,
    enlarge_moon: bool,
    label_reservations: list[QRectF],
    label_candidates: list[dict[str, Any]],
) -> None:
    if not style.show_overlay_info:
        return
    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        scene.celestial_data,
        scene.viewer,
        style.vmag_limit,
        enlarge_moon,
        highlighted_dso,
        highlighted_object,
        style.text_font,
        label_candidates=label_candidates,
        label_reservations=label_reservations,
        mouse_pos=mouse_pos,
        preset=style.visual_preset,
        draw_static_info=True,
        draw_hover_info=False,
        get_text_style_func=render_text.get_text_style,
        draw_outlined_text_func=render_text.draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
    )


def draw_satellite_layer(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    label_candidates: list[dict[str, Any]],
) -> None:
    render_satellites.draw_satellite_overlay(
        painter,
        geometry,
        scene.satellite_overlay_points,
        scene.viewer.view_center,
        opacity=style.satellite_opacity,
        label_candidates=label_candidates,
        preset=style.visual_preset,
        content_fov_deg=_content_fov_deg(scene),
    )


def draw_hover_overlay_layer(
    painter: QPainter,
    *,
    geometry: render_draw.ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    highlighted_object: Any | None,
    highlighted_dso: Any | None,
) -> None:
    line_width_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    if style.show_asterisms and highlighted_object is not None:
        render_asterisms.draw_asterisms(
            painter,
            geometry,
            scene.celestial_data,
            scene.viewer,
            highlighted_object,
            style.text_font,
            preset=style.visual_preset,
            line_width_scale=line_width_scale,
            content_fov_deg=_content_fov_deg(scene),
            draw_base=False,
            draw_highlight=True,
        )
    render_solar_system.draw_hovered_moon_overlay(
        painter,
        geometry,
        scene.celestial_data,
        scene.viewer,
        highlighted_object,
    )
    draw_dso_hover_layer(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
        highlighted_dso=highlighted_dso,
    )
    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        scene.celestial_data,
        scene.viewer,
        style.vmag_limit,
        False,
        highlighted_dso,
        highlighted_object,
        style.text_font,
        preset=style.visual_preset,
        draw_static_info=False,
        draw_hover_info=True,
        get_text_style_func=render_text.get_text_style,
        draw_outlined_text_func=render_text.draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
    )


def draw_label_layer(
    painter: QPainter,
    *,
    style: RenderStyle,
    label_candidates: list[dict[str, Any]],
) -> None:
    render_text._draw_label_candidates(painter, label_candidates, style.text_font)


def draw_status_line(
    painter: QPainter,
    *,
    viewport_rect: QRect,
    style: RenderStyle,
    hud: RenderHudState,
) -> None:
    if not hud.status_message:
        return
    render_text._draw_status_line_text(
        painter=painter,
        message=hud.status_message,
        status_line_font=style.status_line_font,
        viewport_rect=viewport_rect,
        preset=style.visual_preset,
    )
