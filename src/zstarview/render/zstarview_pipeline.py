"""Scenic presentation used by the regular ``zstarview`` viewer."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

from PySide6.QtCore import QRect, QRectF
from PySide6.QtGui import QPainter

from ..gui.composite import SkyCompositorCache
from ..types import ScreenGeometry
from . import pipeline as shared
from .render_types import FrameContext, RenderHudState, RenderSceneData, RenderStyle

ORIENTATION_INTERACTION_STAR_VMAG_LIMIT = 4.0


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
    """Render the regular scenic base scene."""
    win_w, win_h = int(frame.viewport_rect.width()), int(frame.viewport_rect.height())
    star_surface_size = shared.compute_star_render_surface_size(
        win_w,
        win_h,
        frame.geometry.radius * 2,
        style.star_render_expected_width,
    )
    shared._clear_background_layer(painter, frame.viewport_rect)
    shared._draw_background_layer(
        painter,
        geometry=frame.geometry,
        viewport_rect=frame.viewport_rect,
        scene=scene,
        style=style,
        draw_menu_button=not hud.viewport_interaction_mode,
    )
    sky_cloud_style = (
        replace(
            style,
            cloud_disc_alpha=0.0,
            earth_guide_opacity=0.0,
        )
        if shared._simplified_view_active(hud)
        else (
            replace(style, cloud_disc_alpha=0.0)
            if hud.viewport_interaction_mode
            else style
        )
    )
    shared._draw_sky_cloud_layers(
        painter,
        geometry=frame.geometry,
        scene=scene,
        style=sky_cloud_style,
        compositor=compositor,
        star_render_surface_size=star_surface_size,
        simplified_view_active=shared._simplified_view_active(hud),
        fast_mode=hud.viewport_interaction_mode,
    )
    shared._draw_guide_layer(
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
    shared._draw_terrain_layers(
        painter,
        geometry=frame.geometry,
        scene=scene,
        style=style,
        fast_mode=not draw_fast_overlays,
        simplified_view_active=shared._simplified_view_active(hud),
        highlighted_object=None,
        label_reservations=label_reservations,
        label_candidates=local_label_candidates,
    )
    shared._draw_star_layer(
        painter,
        geometry=frame.geometry,
        viewport_rect=frame.viewport_rect,
        scene=scene,
        style=style,
        star_render_surface_size=star_surface_size,
    )
    shared._draw_planet_layer(
        painter,
        geometry=frame.geometry,
        scene=scene,
        style=style,
        enlarge_moon=bool(style.enlarge_moon),
        outline_bright_bodies=str(style.bright_bodies_mode) == "outline",
        label_candidates=local_label_candidates,
    )
    if draw_fast_overlays:
        shared._draw_satellite_layer(
            painter,
            geometry=frame.geometry,
            scene=scene,
            style=style,
            highlighted_satellite=None,
            draw_simplified_labels=shared._simplified_view_labels_visible(hud),
        )
        shared._draw_aircraft_layer(
            painter,
            geometry=frame.geometry,
            scene=scene,
            style=style,
            label_candidates=local_label_candidates,
        )
    if draw_labels:
        shared.render_text._draw_label_candidates(
            painter,
            local_label_candidates,
            style.text_font,
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
    """Render the reduced scenic frame used while the viewport is moving."""
    line_width_scale = shared.compute_star_render_upscale_factor(
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
        shared.render_guides.draw_sky_reference_lines(
            painter,
            geometry,
            scene.viewer,
            scene.celestial_data,
            theme=style.theme,
        )
    shared._draw_star_layer(
        painter,
        geometry=geometry,
        viewport_rect=viewport_rect,
        scene=replace(scene, celestial_data=interaction_celestial_data),
        style=style,
        star_render_surface_size=None,
        draw_vmag_limit=ORIENTATION_INTERACTION_STAR_VMAG_LIMIT,
        fast_mode=True,
    )
    shared._draw_planet_layer(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
        enlarge_moon=bool(style.enlarge_moon),
        outline_bright_bodies=str(style.bright_bodies_mode) == "outline",
        label_candidates=[],
        draw_labels=False,
    )
    shared.render_terrain._draw_terrain_profile_layer(
        painter,
        geometry,
        scene.viewer,
        scene.terrain_horizon_profile,
        scene.terrain_horizon_profile_distances_m,
        spec=shared.render_terrain.TerrainHorizonRenderSpec(
            opacity=style.terrain_horizon_opacity,
            base_width=shared.render_terrain.TERRAIN_HORIZON_FAST_WIDTH,
            far_base_width=shared.render_terrain.TERRAIN_HORIZON_FAR_BASE_WIDTH,
            fg_alpha=shared.render_terrain.terrain_horizon_line_alpha(
                style.terrain_horizon_opacity
                * style.theme.overlays.terrain_horizon.alpha_scale
            ),
            line_width_scale=(
                line_width_scale * style.theme.overlays.terrain_horizon.width_scale
            ),
            color_rgb=style.theme.overlays.terrain_horizon.rgb,
            fast_mode=True,
            distance_widths=True,
        ),
        is_in_fov_func=shared.render_terrain.is_in_fov,
        altaz_to_normalized_xy_func=shared.render_terrain.altaz_to_normalized_xy,
        normalized_to_screen_xy_func=shared.render_terrain.normalized_to_screen_xy,
        split_by_gaps_func=shared.render_terrain.split_by_gaps,
    )
    if shared._should_draw_water_overlay(scene, style):
        water_dots = list(scene.water_overlay_dots) if scene.water_overlay_dots else None
        shared.render_terrain.draw_water_overlay_dots(
            painter,
            geometry,
            scene.viewer,
            water_dots,
            opacity=style.water_overlay_opacity,
            line_width_scale=line_width_scale,
            layer_style=style.theme.overlays.water,
        )
