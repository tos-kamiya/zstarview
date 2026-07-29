"""Scenic presentation used by the regular ``zstarview`` viewer."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np
from PySide6.QtCore import QPointF, QRect, QRectF, Qt
from PySide6.QtGui import QPainter

from ..gui.composite import SkyCompositorCache
from ..types import CelestialObject, ScreenGeometry
from . import instrument_background as render_instrument_background
from . import molecular_cloud_overlay as render_molecular_cloud_overlay
from . import pipeline as shared
from . import sky_disc as render_sky_disc
from .render_types import FrameContext, RenderHudState, RenderSceneData, RenderStyle
from .star_interpolation import (
    STAR_INTERPOLATION_COVERAGE,
    build_star_interpolation_homography,
    should_interpolate_stars,
)

ORIENTATION_INTERACTION_STAR_VMAG_LIMIT = 4.0
TIME_OF_DAY_MARKER_SKY_ALT_DEG = 20.0


def _star_interpolation_matrix(
    *, frame: FrameContext, scene: RenderSceneData
) -> np.ndarray | None:
    if not should_interpolate_stars(frame.sky_update_interval):
        return None
    snapshot_time = scene.celestial_data.star_time or scene.celestial_data.time
    current_time = frame.time_obj
    if snapshot_time is None or current_time is None:
        return None
    elapsed_seconds = float(current_time.unix - snapshot_time.unix)
    if abs(elapsed_seconds) <= 1.0e-9:
        return None
    half_interval_seconds = max(0.0, float(frame.sky_update_interval)) / 2.0
    elapsed_seconds = max(
        -half_interval_seconds,
        min(half_interval_seconds, elapsed_seconds),
    )
    elapsed_seconds *= STAR_INTERPOLATION_COVERAGE
    return build_star_interpolation_homography(
        width_px=int(frame.viewport_rect.width()),
        height_px=int(frame.viewport_rect.height()),
        geometry_center=(
            float(frame.geometry.center[0]),
            float(frame.geometry.center[1]),
        ),
        geometry_radius=float(frame.geometry.radius),
        view_center_altaz_deg=(
            float(frame.viewer.view_center[0]),
            float(frame.viewer.view_center[1]),
        ),
        observer_lat_deg=float(frame.viewer.location[0]),
        edge_fov_deg=float(frame.viewer.edge_fov_deg),
        elapsed_seconds=elapsed_seconds,
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
    draw_stars: bool = True,
    draw_planets: bool = True,
    draw_asterisms: bool = True,
) -> None:
    """Render the regular scenic base scene."""
    win_w, win_h = int(frame.viewport_rect.width()), int(frame.viewport_rect.height())
    star_surface_size = shared.compute_star_render_surface_size(
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
    _draw_sky_cloud_layers(
        painter,
        geometry=frame.geometry,
        scene=scene,
        style=sky_cloud_style,
        compositor=compositor,
        star_render_surface_size=star_surface_size,
        simplified_view_active=shared._simplified_view_active(hud),
        fast_mode=hud.viewport_interaction_mode,
        draw_sky_disc=not hud.viewport_interaction_mode,
        time_obj=frame.time_obj,
    )
    sun_altaz = shared._sun_altaz(scene.celestial_data)
    if sun_altaz is not None:
        render_instrument_background.draw_instrument_time_of_day_marker(
            painter,
            frame.viewport_rect,
            sun_alt_deg=sun_altaz[0],
            bottom_left=(
                hud.time_of_day_marker_bottom_left
                if hud.time_of_day_marker_bottom_left is not None
                else not hud.overlay_info_bottom_left
            ),
            tint_rgba=render_sky_disc.sky_color_at_direction(
                TIME_OF_DAY_MARKER_SKY_ALT_DEG,
                sun_altaz[1],
                sun_altaz,
                alpha=0.6,
                exposure=1.0,
                observer_height_m=float(scene.viewer.observer_height_m),
            ),
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
    _draw_terrain_layers(
        painter,
        geometry=frame.geometry,
        scene=scene,
        style=style,
        fast_mode=not draw_fast_overlays,
        simplified_view_active=shared._simplified_view_active(hud),
        highlighted_object=None,
        label_reservations=label_reservations,
        label_candidates=local_label_candidates,
        draw_asterisms=draw_asterisms,
    )
    if draw_stars:
        shared._draw_star_layer(
            painter,
            geometry=frame.geometry,
            viewport_rect=frame.viewport_rect,
            scene=scene,
            style=style,
            star_render_surface_size=star_surface_size,
            separate_bright_stars=True,
            star_interpolation_matrix=_star_interpolation_matrix(frame=frame, scene=scene),
        )
    if draw_planets:
        shared._draw_planet_layer(
            painter,
            geometry=frame.geometry,
            scene=scene,
            style=style,
            enlarge_moon=bool(style.enlarge_moon),
            outline_bright_bodies=str(style.bright_bodies_mode) == "outline",
            dark_contrast_enabled=float(getattr(style, "sky_disc_alpha", 0.0)) > 0.0,
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
    shared.render_background.draw_radial_background(
        painter,
        QRectF(viewport_rect),
        geometry,
        theme=style.theme,
        edge_fov_deg=float(scene.viewer.edge_fov_deg),
        content_fov_deg=float(scene.viewer.content_fov_deg),
        opaque=not style.show_custom_window_frame,
        altaz_rings_mode=style.sky_disc_altaz_rings,
        view_center=scene.viewer.view_center,
        terrain_profile_altaz=scene.terrain_horizon_profile,
    )
    if style.show_custom_window_frame:
        shared.render_background.draw_window_border(
            painter,
            QRectF(viewport_rect),
            theme=style.theme,
            draw_menu_button=draw_menu_button,
        )


def _clear_background_layer(painter: QPainter, viewport_rect: QRect) -> None:
    painter.save()
    painter.setCompositionMode(QPainter.CompositionMode_Clear)
    painter.fillRect(viewport_rect, Qt.transparent)
    painter.setCompositionMode(QPainter.CompositionMode_SourceOver)
    painter.restore()


def _draw_sky_cloud_layers(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    compositor: SkyCompositorCache,
    star_render_surface_size: tuple[int, int],
    simplified_view_active: bool = False,
    fast_mode: bool = False,
    draw_sky_disc: bool = True,
    time_obj: Any | None = None,
) -> None:
    effective_night_light_opacity = (
        0.0 if simplified_view_active else float(style.night_light_opacity)
    )
    compositor.draw(
        painter,
        geometry,
        scene.sky_disc_image,
        cloud_alpha=style.cloud_disc_alpha,
        density_reference_size=star_render_surface_size,
        view_center=scene.viewer.view_center,
        edge_fov_deg=float(scene.viewer.edge_fov_deg),
        observer_lat_deg=scene.viewer.location[0],
        observer_lon_deg=scene.viewer.location[1],
        observer_height_m=scene.viewer.observer_height_m,
        cloud_altaz_grid=scene.cloud_altaz_grid,
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
        night_light_glow_profile=(
            None if simplified_view_active else scene.night_light_glow_profile
        ),
        earth_guide_opacity=style.earth_guide_opacity,
        earth_guide_visibility_boost=style.earth_guide_visibility_boost,
        night_light_opacity=effective_night_light_opacity,
        ridge_glow_opacity=(
            0.0 if simplified_view_active else float(style.ridge_glow_opacity)
        ),
        night_light_sun_alt_deg=shared._sun_alt_deg(scene.celestial_data),
        molecular_cloud_overlay=render_molecular_cloud_overlay.render_molecular_cloud_overlay(
            width=int(
                painter.viewport().width()
                if callable(getattr(painter, "viewport", None))
                else geometry.radius * 2
            ),
            height=int(
                painter.viewport().height()
                if callable(getattr(painter, "viewport", None))
                else geometry.radius * 2
            ),
            geometry=geometry,
            view_center=scene.viewer.view_center,
            edge_fov_deg=float(scene.viewer.edge_fov_deg),
            content_fov_deg=float(scene.viewer.content_fov_deg),
            sun_alt_deg=shared._sun_alt_deg(scene.celestial_data),
            time_obj=time_obj,
            observer_lat_deg=scene.viewer.location[0],
            observer_lon_deg=scene.viewer.location[1],
            opacity=float(style.akari_ir_bands_opacity),
        ),
        ground_reset_rgba=shared._ground_reset_rgba_for_theme(style.theme),
        theme=style.theme,
        content_fov_deg=float(scene.viewer.content_fov_deg),
        fast_mode=bool(fast_mode),
        draw_sky_disc=bool(draw_sky_disc),
        sky_disc_altaz_rings=str(style.sky_disc_altaz_rings),
    )


def _draw_terrain_layers(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    fast_mode: bool = False,
    simplified_view_active: bool = False,
    highlighted_object: tuple[CelestialObject, QPointF] | None,
    label_reservations: list[QRectF],
    label_candidates: list[dict[str, Any]],
    draw_asterisms: bool = True,
) -> None:
    content_fov_deg = float(scene.viewer.content_fov_deg)
    line_width_scale = shared.compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    simplified_view_content_alpha_scale = 0.4 if simplified_view_active else 1.0
    if style.show_dso:
        shared.render_deep_sky_objects.draw_deep_sky_shapes(
            painter,
            geometry,
            scene.viewer,
            scene.celestial_data,
            opacity_scale=simplified_view_content_alpha_scale,
        )
    if draw_asterisms and style.show_asterisms:
        shared.render_asterisms.draw_asterisms(
            painter,
            geometry,
            scene.viewer,
            scene.celestial_data,
            highlighted_object,
            style.text_font,
            label_reservations,
            label_candidates=label_candidates,
            theme=style.theme,
            line_width_scale=line_width_scale,
            base_line_width_scale=line_width_scale * float(style.asterism_visibility_boost),
            base_line_alpha_scale=float(style.asterism_visibility_boost)
            * simplified_view_content_alpha_scale,
            content_fov_deg=content_fov_deg,
            draw_base=True,
            draw_highlight=False,
        )
    if style.show_guidelines:
        shared.render_guides.draw_sky_reference_lines(
            painter,
            geometry,
            scene.viewer,
            scene.celestial_data,
            theme=style.theme,
        )
    if simplified_view_active:
        shared._draw_main_terrain_profile_layer(
            painter,
            geometry=geometry,
            scene=scene,
            style=style,
            line_width_scale=line_width_scale,
            fast_mode=True,
        )
    if not simplified_view_active:
        shared.render_terrain.draw_terrain_secondary_ridges(
            painter,
            geometry,
            scene.viewer,
            scene.terrain_secondary_ridges_altaz_layers,
            scene.terrain_secondary_ridges_distances_m_layers,
            opacity=max(0.0, float(style.terrain_horizon_opacity) * 0.72),
            line_width_scale=line_width_scale,
            layer_style=style.theme.overlays.terrain_horizon,
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
                fast_mode=fast_mode,
            )
            shared.render_terrain.draw_water_overlay_polylines(
                painter,
                geometry,
                scene.viewer,
                list(scene.water_overlay_polylines)
                if scene.water_overlay_polylines
                else None,
                opacity=style.water_overlay_opacity * 0.85,
                line_width_scale=line_width_scale,
                layer_style=style.theme.overlays.water,
            )
        shared._draw_urban_outline_layer(
            painter,
            geometry=geometry,
            scene=scene,
            style=style,
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
        dark_contrast_enabled=float(getattr(style, "sky_disc_alpha", 0.0)) > 0.0,
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
        shared.render_terrain.draw_water_overlay_polylines(
            painter,
            geometry,
            scene.viewer,
            list(scene.water_overlay_polylines)
            if scene.water_overlay_polylines
            else None,
            opacity=style.water_overlay_opacity * 0.85,
            line_width_scale=line_width_scale,
            layer_style=style.theme.overlays.water,
        )
