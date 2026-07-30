"""Instrument presentation used by ``zstarview-atlas``."""

from __future__ import annotations

from typing import Any

from PySide6.QtCore import QRect
from PySide6.QtGui import QPainter

from ..gui.composite import SkyCompositorCache
from ..types import ScreenGeometry, ViewProjection
from . import deep_sky_objects as render_deep_sky_objects
from . import earth_guide as render_earth_guide
from . import guides as render_guides
from . import instrument_background as render_instrument_background
from . import pipeline as shared
from . import terrain as render_terrain
from .render_types import FrameContext, RenderHudState, RenderSceneData, RenderStyle


class InstrumentSkyPresentation:
    """Stable positional presentation used by zstarview-atlas."""

    def render_base_scene_into_painter(
        self,
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
        local_label_candidates = label_candidates if label_candidates is not None else []
        simplified_view_active = shared._simplified_view_active(hud)
        simplified_view_labels_visible = shared._simplified_view_labels_visible(hud)
        render_instrument_background.draw_instrument_background(
            painter,
            frame.viewport_rect,
            theme=style.theme,
        )
        if not simplified_view_active:
            render_terrain.draw_ground_tint(
                painter,
                frame.geometry,
                scene.viewer,
                scene.terrain_horizon_profile,
                opacity=style.ground_tint_opacity,
            )
        _draw_instrument_guide_layer(
            painter,
            geometry=frame.geometry,
            viewport_rect=frame.viewport_rect,
            scene=scene,
            style=style,
            draw_direction_labels=draw_direction_labels,
        )
        _draw_instrument_context_layers(
            painter,
            geometry=frame.geometry,
            scene=scene,
            style=style,
            label_candidates=local_label_candidates,
            simplified_view_active=simplified_view_active,
            time_obj=frame.time_obj,
        )
        _draw_instrument_cloud_layer(
            painter,
            geometry=frame.geometry,
            viewport_rect=frame.viewport_rect,
            scene=scene,
            style=style,
            compositor=compositor,
        )
        render_instrument_background.draw_instrument_time_of_day_marker(
            painter,
            frame.viewport_rect,
            sun_alt_deg=shared._sun_alt_deg(scene.celestial_data),
            bottom_left=(
                hud.time_of_day_marker_bottom_left
                if hud.time_of_day_marker_bottom_left is not None
                else not hud.overlay_info_bottom_left
            ),
        )
        if draw_stars:
            shared._draw_star_layer(
                painter,
                geometry=frame.geometry,
                viewport_rect=frame.viewport_rect,
                scene=scene,
                style=style,
                star_render_surface_size=None,
                fast_mode=False,
            )
        if simplified_view_labels_visible:
            shared._draw_simplified_named_star_labels(
                painter,
                geometry=frame.geometry,
                viewport_rect=frame.viewport_rect,
                scene=scene,
                style=style,
                highlighted_object=None,
            )
        if draw_planets:
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
                draw_simplified_labels=False,
            )
            shared._draw_aircraft_layer(
                painter,
                geometry=frame.geometry,
                scene=scene,
                style=style,
                label_candidates=local_label_candidates,
            )
        if draw_labels and (not simplified_view_active or simplified_view_labels_visible):
            shared.render_text._draw_label_candidates(
                painter,
                local_label_candidates,
                style.text_font,
            )


def _draw_instrument_guide_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    style: RenderStyle,
    draw_direction_labels: bool = True,
) -> None:
    if not style.show_guidelines:
        return
    render_guides.draw_direction_grid_overlay(
        painter,
        geometry,
        scene.viewer,
        (int(viewport_rect.width()), int(viewport_rect.height())),
        theme=style.theme,
    )
    render_guides.draw_sky_reference_lines(
        painter,
        geometry,
        scene.viewer,
        scene.celestial_data,
        theme=style.theme,
    )
    shared._draw_guide_layer(
        painter,
        geometry=geometry,
        viewport_rect=viewport_rect,
        scene=scene,
        style=style,
        draw_direction_labels=draw_direction_labels,
    )


def _draw_instrument_context_layers(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    label_candidates: list[dict[str, Any]],
    simplified_view_active: bool = False,
    time_obj: object | None = None,
) -> None:
    line_width_scale = shared.compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    if style.show_dso:
        render_deep_sky_objects.draw_deep_sky_shapes(
            painter,
            geometry,
            scene.viewer,
            scene.celestial_data,
            theme=style.theme,
        )
    if style.show_asterisms:
        shared.render_asterisms.draw_asterisms(
            painter,
            geometry,
            scene.viewer,
            scene.celestial_data,
            None,
            style.text_font,
            [],
            label_candidates=label_candidates,
            theme=style.theme,
            line_width_scale=line_width_scale,
            base_line_width_scale=line_width_scale,
            base_line_alpha_scale=0.55,
            content_fov_deg=float(scene.viewer.content_fov_deg),
            draw_base=True,
            draw_highlight=False,
        )
    shared._draw_main_terrain_profile_layer(
        painter,
        geometry=geometry,
        scene=scene,
        style=style,
        line_width_scale=line_width_scale,
        fast_mode=False,
    )
    if shared._should_draw_water_overlay(scene, style):
        render_terrain.draw_water_overlay_dots(
            painter,
            geometry,
            scene.viewer,
            list(scene.water_overlay_dots) if scene.water_overlay_dots else None,
            opacity=style.water_overlay_opacity,
            line_width_scale=line_width_scale,
            layer_style=style.theme.overlays.water,
            fast_mode=False,
        )
        render_terrain.draw_water_overlay_polylines(
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
    if not simplified_view_active:
        shared._draw_urban_outline_layer(
            painter,
            geometry=geometry,
            scene=scene,
            style=style,
            time_obj=time_obj,
        )
    if style.earth_guide_opacity > 0.0:
        render_earth_guide.draw_earth_guide(
            painter,
            geometry=geometry,
            viewer_data=scene.viewer,
            terrain_profile_altaz=scene.terrain_horizon_profile,
            earth_guide_opacity=style.earth_guide_opacity,
            visibility_boost=style.earth_guide_visibility_boost,
            single_line=True,
            layer_style=style.theme.overlays.earth_guide,
        )


def _draw_instrument_cloud_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    style: RenderStyle,
    compositor: SkyCompositorCache,
) -> None:
    """Draw Atlas clouds directly over its flat background."""
    grid = scene.cloud_altaz_grid
    if grid is None or float(style.cloud_disc_alpha) <= 0.0:
        return

    projection = ViewProjection(
        view_center=tuple(float(value) for value in scene.viewer.view_center),
        edge_fov_deg=float(scene.viewer.edge_fov_deg),
        content_fov_deg=float(scene.viewer.content_fov_deg),
    )
    cloud_image, missing_image = compositor.render_atlas_cloud_layer(
        width=int(viewport_rect.width()),
        height=int(viewport_rect.height()),
        geometry=geometry,
        projection=projection,
        grid=grid,
        missing_mask=scene.cloud_missing_mask,
        target_stripes=compositor.cloud_target_stripes,
        width_factor=compositor.cloud_stripe_width_factor,
        opacity=float(style.cloud_disc_alpha),
        style=style.theme.overlays.cloud,
    )
    origin = viewport_rect.topLeft()
    if missing_image is not None:
        painter.drawImage(origin, missing_image)
    painter.drawImage(origin, cloud_image)
