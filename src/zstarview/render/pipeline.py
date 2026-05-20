from __future__ import annotations

import math
from dataclasses import dataclass, replace
from datetime import datetime
from typing import Any

import astropy.time
import numpy as np
from PySide6.QtCore import QPoint, QPointF, QRect, QRectF, Qt
from PySide6.QtGui import QFont, QImage, QPainter

from ..aircraft.types import AircraftSnapshot
from ..gui.composite import CloudAmountField, SkyCompositorCache
from ..night_lights import NightLightGlowProfile
from ..paths import THEME_STYLES_BY_PRESET, ThemeStyle
from ..satellites.types import SatelliteOmmRecord, SatelliteOverlayPoint
from ..search.models import SearchJumpTarget
from ..types import (
    CelestialData,
    CelestialObject,
    ScreenGeometry,
    StarsTable,
    UrbanOutlinePolyline,
    ViewerData,
)
from ..water_overlay import WaterOverlayPoint
from . import aircraft as render_aircraft
from . import asterisms as render_asterisms
from . import background as render_background
from . import deep_sky_objects as render_deep_sky_objects
from . import guides as render_guides
from . import overlay_info as render_overlay_info
from . import satellites as render_satellites
from . import search_overlay as render_search_overlay
from . import solar_system as render_solar_system
from . import stars as render_stars
from . import terrain as render_terrain
from . import text as render_text

ORIENTATION_INTERACTION_STAR_VMAG_LIMIT = 4.0


def _sun_alt_deg(celestial_data: CelestialData) -> float | None:
    for body in celestial_data.planets:
        if body.name == "sun":
            return float(body.alt)
    return None


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
class FrameContext:
    viewer: ViewerData
    time_obj: astropy.time.Time | None
    geometry: ScreenGeometry
    viewport_rect: QRect


@dataclass(frozen=True)
class RenderSceneData:
    viewer: ViewerData
    celestial_data: CelestialData
    sky_disc_image: QImage | None
    cloud_image: np.ndarray | None
    cloud_missing_mask: np.ndarray | None
    cloud_amount_field: CloudAmountField | None
    terrain_horizon_profile: list[tuple[float, float]] | None
    terrain_horizon_profile_distances_m: list[float] | None
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None
    terrain_secondary_ridges_distances_m_layers: list[list[float]] | None
    urban_outlines: list[UrbanOutlinePolyline] | None
    satellite_element_epoch_utc: datetime | None = None
    satellite_records_by_group: dict[str, list[SatelliteOmmRecord]] | None = None
    aircraft_snapshots: list[AircraftSnapshot] | None = None
    time_obj: astropy.time.Time | None = None
    night_light_glow_profile: NightLightGlowProfile | None = None
    water_overlay_dots: list[WaterOverlayPoint] | None = None


@dataclass(frozen=True)
class RenderStyle:
    visual_preset: str
    text_font: QFont
    status_line_font: QFont
    show_background_gradient: bool
    show_custom_window_frame: bool
    show_observation_info: bool
    show_dso: bool
    show_asterisms: bool
    show_guidelines: bool
    enlarge_moon: bool
    bright_bodies_mode: str
    star_base_radius: float
    star_visibility_boost: float
    asterism_visibility_boost: float
    earth_guide_visibility_boost: float
    vmag_limit: float
    sky_disc_altaz_rings: str
    sky_disc_altaz_rings_hover: str
    cloud_disc_alpha: float
    satellite_opacity: float
    terrain_horizon_opacity: float
    earth_guide_opacity: float
    night_light_opacity: float = 0.02
    urban_outline_opacity: float = 0.2
    show_urban_outline_layer: bool = True
    water_overlay_opacity: float = 0.12
    aircraft_opacity: float = 0.5
    star_render_expected_width: int = 600
    theme: ThemeStyle = THEME_STYLES_BY_PRESET["night"]


def _should_draw_water_overlay(scene: RenderSceneData, style: RenderStyle) -> bool:
    return float(style.terrain_horizon_opacity) > 0.0 and scene.terrain_horizon_profile is not None


def _terrain_horizon_water_overlay_dots(scene: RenderSceneData) -> list[WaterOverlayPoint] | None:
    water_dots = scene.water_overlay_dots
    return list(water_dots) if water_dots else None


@dataclass(frozen=True)
class RenderHudState:
    mouse_pos: QPoint | None
    overlay_info_bottom_left: bool
    viewport_interaction_mode: bool
    viewport_interaction_stars: StarsTable | None
    status_message: str | None


def _content_fov_deg(scene: RenderSceneData) -> float:
    return float(scene.viewer.content_fov_deg)


def _bright_bodies_mode(style: RenderStyle) -> str:
    return str(style.bright_bodies_mode)


def _window_size(viewport_rect: QRect) -> tuple[int, int]:
    return (int(viewport_rect.width()), int(viewport_rect.height()))


def _resolve_frame_context(
    *,
    frame: FrameContext | None,
    scene: RenderSceneData,
    geometry: ScreenGeometry | None,
    viewport_rect: QRect | None,
) -> FrameContext:
    if frame is not None:
        return frame
    if geometry is None or viewport_rect is None:
        raise TypeError("frame or geometry/viewport_rect must be provided")
    return FrameContext(
        viewer=scene.viewer,
        time_obj=scene.time_obj,
        geometry=geometry,
        viewport_rect=viewport_rect,
    )


def render_base_scene_into_painter(
    painter: QPainter,
    *,
    frame: FrameContext,
    scene: RenderSceneData,
    style: RenderStyle,
    hud: RenderHudState,
    compositor: SkyCompositorCache,
    draw_fast_overlays: bool = True,
    label_candidates: list[dict[str, Any]] | None = None,
    draw_labels: bool = True,
    draw_direction_labels: bool = True,
) -> None:
    win_w, win_h = _window_size(frame.viewport_rect)
    star_surface_size = compute_star_render_surface_size(
        win_w,
        win_h,
        frame.geometry.radius * 2,
        style.star_render_expected_width,
    )
    _clear_background_layer(painter, frame.viewport_rect)
    _draw_background_layer(
        painter,
        geometry=frame.geometry,
        viewport_rect=frame.viewport_rect,
        scene=scene,
        style=style,
        draw_menu_button=not hud.viewport_interaction_mode,
    )
    sky_cloud_style = (
        replace(style, cloud_disc_alpha=0.0) if hud.viewport_interaction_mode else style
    )
    sky_cloud_scene = (
        replace(scene, sky_disc_image=None) if hud.viewport_interaction_mode else scene
    )
    _draw_sky_cloud_layers(
        painter,
        geometry=frame.geometry,
        scene=sky_cloud_scene,
        style=sky_cloud_style,
        compositor=compositor,
        star_render_surface_size=star_surface_size,
        fast_mode=hud.viewport_interaction_mode,
    )
    _draw_guide_layer(
        painter,
        geometry=frame.geometry,
        viewport_rect=frame.viewport_rect,
        scene=scene,
        style=style,
        draw_direction_labels=draw_direction_labels,
    )
    if hud.viewport_interaction_mode:
        _draw_viewport_interaction_layers(
            painter,
            geometry=frame.geometry,
            viewport_rect=frame.viewport_rect,
            scene=scene,
            style=style,
            hud=hud,
        )
        return

    label_reservations: list[QRectF] = []
    local_label_candidates = label_candidates if label_candidates is not None else []
    _draw_terrain_layers(
        painter,
        geometry=frame.geometry,
        scene=scene,
        style=style,
        highlighted_object=None,
        label_reservations=label_reservations,
        label_candidates=local_label_candidates,
    )
    _draw_star_layer(
        painter,
        geometry=frame.geometry,
        viewport_rect=frame.viewport_rect,
        scene=scene,
        style=style,
        star_render_surface_size=star_surface_size,
    )
    _draw_planet_layer(
        painter,
        geometry=frame.geometry,
        scene=scene,
        style=style,
        enlarge_moon=bool(style.enlarge_moon),
        outline_bright_bodies=_bright_bodies_mode(style) == "outline",
        label_candidates=local_label_candidates,
    )
    if draw_fast_overlays:
        _draw_satellite_layer(
            painter,
            geometry=frame.geometry,
            scene=scene,
            style=style,
            highlighted_satellite=None,
            label_candidates=local_label_candidates,
        )
        _draw_aircraft_layer(
            painter,
            geometry=frame.geometry,
            scene=scene,
            style=style,
            label_candidates=local_label_candidates,
        )
    if draw_labels:
        render_text._draw_label_candidates(painter, local_label_candidates, style.text_font)


def render_fast_overlay_layers_into_painter(
    painter: QPainter,
    *,
    frame: FrameContext,
    scene: RenderSceneData,
    style: RenderStyle,
    highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None = None,
    label_candidates: list[dict[str, Any]] | None = None,
    draw_labels: bool = True,
) -> None:
    """Draw dynamic satellite/aircraft overlays and their labels."""
    if style.satellite_opacity <= 0.0 and style.aircraft_opacity <= 0.0:
        return
    local_label_candidates = label_candidates if label_candidates is not None else []
    _draw_satellite_layer(
        painter,
        geometry=frame.geometry,
        scene=scene,
        style=style,
        highlighted_satellite=highlighted_satellite,
        label_candidates=local_label_candidates,
    )
    _draw_aircraft_layer(
        painter,
        geometry=frame.geometry,
        scene=scene,
        style=style,
        label_candidates=local_label_candidates,
    )
    if draw_labels:
        render_text._draw_label_candidates(painter, local_label_candidates, style.text_font)


def render_hud_overlay_into_painter(
    painter: QPainter,
    *,
    frame: FrameContext,
    scene: RenderSceneData,
    style: RenderStyle,
    hud: RenderHudState,
    highlighted_object: tuple[CelestialObject, QPointF] | None,
    highlighted_dso: tuple[CelestialObject, QPointF] | None,
    highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None = None,
    label_candidates: list[dict[str, Any]] | None = None,
    search_overlay_target: SearchJumpTarget | None = None,
) -> None:
    if hud.viewport_interaction_mode:
        _draw_status_line(
            painter,
            viewport_rect=frame.viewport_rect,
            style=style,
            hud=hud,
        )
        return

    _draw_hover_overlay_layer(
        painter,
        geometry=frame.geometry,
        viewport_rect=frame.viewport_rect,
        scene=scene,
        style=style,
        mouse_pos=hud.mouse_pos,
        highlighted_object=highlighted_object,
        highlighted_dso=highlighted_dso,
        highlighted_satellite=highlighted_satellite,
        label_candidates=label_candidates,
    )
    if search_overlay_target is not None:
        render_search_overlay.draw_search_target_overlay(
            painter,
            frame.geometry,
            search_overlay_target,
            viewer_data=scene.viewer,
            text_font=style.text_font,
            draw_marker=True,
            draw_label=True,
            marker_scale=compute_star_render_upscale_factor(
                frame.geometry.radius * 2,
                style.star_render_expected_width,
            ),
            label_candidates=label_candidates,
            theme=style.theme,
        )
    if label_candidates:
        render_text._draw_label_candidates(painter, label_candidates, style.text_font)
    _draw_overlay_layer(
        painter,
        geometry=frame.geometry,
        viewport_rect=frame.viewport_rect,
        scene=scene,
        style=style,
        mouse_pos=hud.mouse_pos,
        overlay_info_bottom_left=hud.overlay_info_bottom_left,
        highlighted_object=None,
        highlighted_dso=None,
        enlarge_moon=bool(style.enlarge_moon),
        label_reservations=[],
        label_candidates=label_candidates,
    )
    _draw_status_line(
        painter,
        viewport_rect=frame.viewport_rect,
        style=style,
        hud=hud,
    )


def _draw_viewport_interaction_layers(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
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
    _draw_star_layer(
        painter,
        geometry=geometry,
        viewport_rect=viewport_rect,
        scene=replace(scene, celestial_data=interaction_celestial_data),
        style=style,
        star_render_surface_size=None,
        draw_vmag_limit=ORIENTATION_INTERACTION_STAR_VMAG_LIMIT,
        fast_mode=True,
    )
    render_terrain.draw_terrain_horizon_fast(
        painter,
        geometry,
        scene.terrain_horizon_profile,
        scene.terrain_horizon_profile_distances_m,
        scene.viewer.view_center,
        opacity=style.terrain_horizon_opacity,
        line_width_scale=line_width_scale,
        edge_fov_deg=float(scene.viewer.edge_fov_deg),
        content_fov_deg=_content_fov_deg(scene),
    )
    if _should_draw_water_overlay(scene, style):
        water_dots = _terrain_horizon_water_overlay_dots(scene)
        render_terrain.draw_water_overlay_dots(
            painter,
            geometry,
            water_dots,
            scene.viewer.view_center,
            opacity=style.water_overlay_opacity,
            line_width_scale=line_width_scale,
            edge_fov_deg=float(scene.viewer.edge_fov_deg),
            content_fov_deg=_content_fov_deg(scene),
        )


def _clear_background_layer(painter: QPainter, viewport_rect: QRect) -> None:
    painter.save()
    painter.setCompositionMode(QPainter.CompositionMode_Clear)
    painter.fillRect(viewport_rect, Qt.transparent)
    painter.setCompositionMode(QPainter.CompositionMode_SourceOver)
    painter.restore()


def _draw_background_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    style: RenderStyle,
    draw_menu_button: bool = True,
) -> None:
    if not style.show_background_gradient:
        return
    render_background.draw_radial_background(
        painter,
        QRectF(viewport_rect),
        geometry,
        theme=style.theme,
        edge_fov_deg=float(scene.viewer.edge_fov_deg),
        content_fov_deg=_content_fov_deg(scene),
        opaque=not style.show_custom_window_frame,
        altaz_rings_mode=style.sky_disc_altaz_rings,
        view_center=scene.viewer.view_center,
        terrain_profile_altaz=scene.terrain_horizon_profile,
    )
    if style.show_custom_window_frame:
        render_background.draw_window_border(
            painter,
            QRectF(viewport_rect),
            theme=style.theme,
            draw_menu_button=draw_menu_button,
        )


def _draw_guide_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    style: RenderStyle,
    draw_direction_labels: bool = True,
) -> None:
    """Draw guide annotations that should float above sky/cloud but below scene overlays."""
    if not style.show_guidelines:
        return
    if style.sky_disc_altaz_rings == "altaz":
        render_guides.draw_direction_grid_overlay(
            painter,
            geometry,
            _window_size(viewport_rect),
            scene.viewer,
        )
    if draw_direction_labels:
        render_guides.draw_direction_labels(
            painter,
            geometry,
            scene.viewer,
            style.text_font,
            None,
            theme=style.theme,
        )
    render_guides.draw_zenith_marker(
        painter,
        geometry,
        scene.viewer,
        theme=style.theme,
    )


def _draw_sky_cloud_layers(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    compositor: SkyCompositorCache,
    star_render_surface_size: tuple[int, int],
    fast_mode: bool = False,
) -> None:
    compositor.draw(
        painter,
        geometry,
        scene.sky_disc_image,
        scene.cloud_image,
        cloud_alpha=style.cloud_disc_alpha,
        density_reference_size=star_render_surface_size,
        view_center=scene.viewer.view_center,
        edge_fov_deg=float(scene.viewer.edge_fov_deg),
        observer_lat_deg=scene.viewer.location[0],
        observer_lon_deg=scene.viewer.location[1],
        observer_height_m=scene.viewer.observer_height_m,
        cloud_amount_field=scene.cloud_amount_field,
        missing_mask=scene.cloud_missing_mask,
        show_guidelines=style.show_guidelines,
        terrain_profile_altaz=(
            scene.terrain_horizon_profile
            if style.terrain_horizon_opacity > 0.0
            else None
        ),
        terrain_profile_distances_m=(
            scene.terrain_horizon_profile_distances_m
            if style.terrain_horizon_opacity > 0.0
            else None
        ),
        terrain_secondary_ridges_altaz_layers=(
            scene.terrain_secondary_ridges_altaz_layers
            if style.terrain_horizon_opacity > 0.0
            else None
        ),
        terrain_secondary_ridges_distances_m_layers=(
            scene.terrain_secondary_ridges_distances_m_layers
            if style.terrain_horizon_opacity > 0.0
            else None
        ),
        night_light_glow_profile=scene.night_light_glow_profile,
        earth_guide_opacity=style.earth_guide_opacity,
        earth_guide_visibility_boost=style.earth_guide_visibility_boost,
        night_light_opacity=style.night_light_opacity,
        night_light_sun_alt_deg=_sun_alt_deg(scene.celestial_data),
        ground_reset_rgba=style.theme.window_background.inner_rgba,
        theme=style.theme,
        content_fov_deg=_content_fov_deg(scene),
        fast_mode=bool(fast_mode),
        sky_disc_altaz_rings=str(style.sky_disc_altaz_rings),
    )
def _draw_terrain_layers(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    highlighted_object: tuple[CelestialObject, QPointF] | None,
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
            theme=style.theme,
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
            theme=style.theme,
            line_width_scale=line_width_scale,
            base_line_width_scale=line_width_scale * float(style.asterism_visibility_boost),
            base_line_alpha_scale=float(style.asterism_visibility_boost),
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
    render_terrain.draw_terrain_secondary_ridges(
        painter,
        geometry,
        scene.viewer,
        scene.terrain_secondary_ridges_altaz_layers,
        scene.terrain_secondary_ridges_distances_m_layers,
        opacity=max(0.0, float(style.terrain_horizon_opacity) * 0.72),
        line_width_scale=line_width_scale,
    )
    if _should_draw_water_overlay(scene, style):
        water_dots = _terrain_horizon_water_overlay_dots(scene)
        render_terrain.draw_water_overlay_dots(
            painter,
            geometry,
            water_dots,
            scene.viewer.view_center,
            opacity=style.water_overlay_opacity,
            line_width_scale=line_width_scale,
            edge_fov_deg=float(scene.viewer.edge_fov_deg),
            content_fov_deg=content_fov_deg,
        )
    _draw_urban_outline_layer(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
    )


def _draw_dso_hover_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    highlighted_dso: tuple[CelestialObject, QPointF] | None,
) -> None:
    if not style.show_dso:
        return
    render_deep_sky_objects.draw_dso_hover_info(
        painter,
        geometry,
        scene.viewer,
        highlighted_dso,
        style.text_font,
        theme=style.theme,
    )


def _draw_urban_outline_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
) -> None:
    if not style.show_urban_outline_layer:
        return
    render_terrain.draw_urban_outlines(
        painter,
        geometry,
        scene.urban_outlines,
        scene.viewer.view_center,
        opacity=style.urban_outline_opacity,
        line_width_scale=1.0,
        edge_fov_deg=float(scene.viewer.edge_fov_deg),
        content_fov_deg=_content_fov_deg(scene),
    )


def _draw_aircraft_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
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
        scene.aircraft_snapshots,
        viewer_data=scene.viewer,
        time_obj=scene.time_obj,
        opacity=style.aircraft_opacity,
        line_width_scale=line_width_scale,
        label_candidates=label_candidates,
        theme=style.theme,
    )


def _draw_star_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    style: RenderStyle,
    star_render_surface_size: tuple[int, int] | None = None,
    draw_vmag_limit: float | None = None,
    fast_mode: bool = False,
) -> None:
    draw_data = scene.celestial_data
    win_w, win_h = _window_size(viewport_rect)
    outline_render_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    outline_bright_bodies = _bright_bodies_mode(style) == "outline"
    low_w, low_h = (
        (win_w, win_h)
        if outline_bright_bodies
        else compute_star_render_surface_size(
            win_w,
            win_h,
            geometry.radius * 2,
            style.star_render_expected_width,
        )
        if star_render_surface_size is None
        else (max(1, int(star_render_surface_size[0])), max(1, int(star_render_surface_size[1])))
    )
    content_fov_deg = _content_fov_deg(scene)
    if low_w == win_w and low_h == win_h:
        draw_stars = (
            render_stars.draw_stars_fast if fast_mode else render_stars.draw_stars_normal
        )
        draw_stars(
            painter,
            geometry,
            draw_data,
            scene.viewer,
            style.star_base_radius,
            visibility_boost=style.star_visibility_boost,
            outline_bright_bodies=outline_bright_bodies,
            outline_render_scale=outline_render_scale,
            draw_vmag_limit=draw_vmag_limit
            if draw_vmag_limit is not None
            else style.vmag_limit,
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
    low_geometry = ScreenGeometry(
        center=(
            int(round(geometry.center[0] * sx)),
            int(round(geometry.center[1] * sy)),
        ),
        radius=max(1, int(round(geometry.radius * min(sx, sy)))),
    )
    draw_stars = render_stars.draw_stars_fast if fast_mode else render_stars.draw_stars_normal
    draw_stars(
        low_painter,
        low_geometry,
        draw_data,
        scene.viewer,
        style.star_base_radius,
        visibility_boost=style.star_visibility_boost,
        outline_bright_bodies=outline_bright_bodies,
        outline_render_scale=outline_render_scale,
        draw_vmag_limit=draw_vmag_limit
        if draw_vmag_limit is not None
        else style.vmag_limit,
        viewport_size=(low_w, low_h),
        content_fov_deg=content_fov_deg,
    )
    low_painter.end()

    painter.save()
    painter.setRenderHint(QPainter.SmoothPixmapTransform, False)
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
    painter.drawImage(viewport_rect, low_img)
    painter.restore()


def _draw_planet_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    enlarge_moon: bool,
    outline_bright_bodies: bool = False,
    label_candidates: list[dict[str, Any]],
) -> None:
    marker_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    render_solar_system.draw_solar_system_bodies(
        painter,
        geometry,
        scene.celestial_data,
        scene.viewer,
        enlarge_moon,
        outline_bright_bodies=outline_bright_bodies,
        text_font=style.text_font,
        label_candidates=label_candidates,
        draw_labels=True,
        theme=style.theme,
        edge_fov_deg=float(scene.viewer.edge_fov_deg),
        content_fov_deg=_content_fov_deg(scene),
        marker_scale=marker_scale,
    )


def _draw_overlay_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    style: RenderStyle,
    mouse_pos: QPoint | None,
    overlay_info_bottom_left: bool,
    highlighted_object: tuple[CelestialObject, QPointF] | None,
    highlighted_dso: tuple[CelestialObject, QPointF] | None,
    enlarge_moon: bool,
    label_reservations: list[QRectF],
    label_candidates: list[dict[str, Any]],
) -> None:
    if not style.show_observation_info:
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
        viewport_rect=viewport_rect,
        mouse_pos=mouse_pos,
        bottom_left=overlay_info_bottom_left,
        theme=style.theme,
        draw_static_info=True,
        draw_hover_info=False,
        draw_outlined_text_func=render_text.draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
    )


def _draw_satellite_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None,
    label_candidates: list[dict[str, Any]],
) -> None:
    render_satellites.draw_satellite_overlay(
        painter,
        geometry,
        scene.satellite_records_by_group,
        viewer_data=scene.viewer,
        time_obj=scene.time_obj,
        opacity=style.satellite_opacity,
        highlighted_satellite=(
            highlighted_satellite[0] if highlighted_satellite is not None else None
        ),
        label_candidates=label_candidates,
        theme=style.theme,
        marker_scale=compute_star_render_upscale_factor(
            geometry.radius * 2,
            style.star_render_expected_width,
        ),
    )


def _draw_hover_overlay_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    style: RenderStyle,
    mouse_pos: QPoint | None = None,
    highlighted_object: tuple[CelestialObject, QPointF] | None,
    highlighted_dso: tuple[CelestialObject, QPointF] | None,
    highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None = None,
    label_candidates: list[dict[str, Any]] | None = None,
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
            theme=style.theme,
            line_width_scale=line_width_scale,
            content_fov_deg=_content_fov_deg(scene),
            draw_base=False,
            draw_highlight=True,
            label_candidates=label_candidates,
        )
    render_solar_system.draw_hovered_moon_overlay(
        painter,
        geometry,
        scene.celestial_data,
        scene.viewer,
        highlighted_object,
        marker_scale=line_width_scale,
        outline_bright_bodies=_bright_bodies_mode(style) == "outline",
        theme=style.theme,
    )
    _draw_dso_hover_layer(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
        highlighted_dso=highlighted_dso,
    )
    direction_hover = None
    if style.show_guidelines and mouse_pos is not None:
        direction_hover = render_guides.resolve_direction_marker_hover(
            geometry,
            scene.viewer,
            mouse_pos,
        )
    if direction_hover is not None:
        if style.sky_disc_altaz_rings_hover == "dimalt":
            dimalt_sample_color = render_background.sample_background_disc_edge_color(
                QRectF(viewport_rect),
                geometry,
                theme=style.theme,
                edge_fov_deg=float(scene.viewer.edge_fov_deg),
                content_fov_deg=_content_fov_deg(scene),
                opaque=not style.show_custom_window_frame,
            )
            render_background.draw_altitude_ring_overlay(
                painter,
                QRectF(viewport_rect),
                geometry,
                view_center=scene.viewer.view_center,
                theme=style.theme,
                edge_fov_deg=float(scene.viewer.edge_fov_deg),
                content_fov_deg=_content_fov_deg(scene),
                ring_color=render_background.dimalt_ring_pen_color_from_color(
                    dimalt_sample_color
                ),
            )
        elif style.sky_disc_altaz_rings_hover == "altaz":
            render_guides.draw_direction_grid_overlay(
                painter,
                geometry,
                _window_size(viewport_rect),
                scene.viewer,
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
        highlighted_satellite,
        label_candidates=label_candidates,
        theme=style.theme,
        draw_static_info=False,
        draw_hover_info=True,
        draw_outlined_text_func=render_text.draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
    )

def _draw_status_line(
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
        theme=style.theme,
    )
