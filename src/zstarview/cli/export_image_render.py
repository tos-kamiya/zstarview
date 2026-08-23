"""Render an export-image frame from prepared scene inputs."""

from __future__ import annotations

from PySide6.QtCore import QPoint, QRect
from PySide6.QtGui import QFont, QFontDatabase, QImage, QPainter

from ..gui.composite import SkyCompositorCache
from ..gui.window_inputs import (
    PreparedWindowCatalogs,
    SkyWindowRuntimeOptions,
    SkyWindowUserOptions,
)
from ..paths import (
    CLOUD_MISSING_TINT_RGBA,
    OVERLAY_FONT_SIZE_DEFAULT,
    STATUS_LINE_FONT_SIZE,
    TEXT_FONT_PATH,
    ThemeStyle,
)
from ..render import geometry as render_geometry
from ..render import guides as render_guides
from ..render import text as render_text
from ..render.pipeline import (
    FrameContext,
    RenderHudState,
    RenderSceneData,
    RenderStyle,
    compute_star_render_upscale_factor,
)
from ..render.search_overlay import draw_search_target_overlay
from ..search.models import SearchJumpTarget
from .export_image_support import DEFAULT_EXPORT_IMAGE_SKY_UPDATE_INTERVAL, host


def _load_fonts(
    overlay_font_size: float = float(OVERLAY_FONT_SIZE_DEFAULT),
) -> tuple[QFont, QFont]:
    text_font_id = QFontDatabase.addApplicationFont(TEXT_FONT_PATH)
    text_font_family = QFontDatabase.applicationFontFamilies(text_font_id)[0]
    text_font = QFont(text_font_family)
    text_font.setPointSizeF(float(overlay_font_size))
    status_line_font = QFont(text_font_family)
    status_line_font.setPointSizeF(float(STATUS_LINE_FONT_SIZE))
    return (text_font, status_line_font)

def _build_compositor(
    runtime_options: SkyWindowRuntimeOptions, user_options: SkyWindowUserOptions
) -> SkyCompositorCache:
    target_stripes, width_factor = runtime_options.cloud_stripe_style
    missing_tint_alpha = int(round(255.0 * runtime_options.cloud_missing_tint_opacity))
    missing_tint_rgba = (
        int(CLOUD_MISSING_TINT_RGBA[0]),
        int(CLOUD_MISSING_TINT_RGBA[1]),
        int(CLOUD_MISSING_TINT_RGBA[2]),
        missing_tint_alpha,
    )
    return SkyCompositorCache(
        cloud_target_stripes=int(target_stripes),
        cloud_stripe_width_factor=float(width_factor),
        cloud_stripe_mode=runtime_options.cloud_stripe_mode,
        missing_tint_rgba=missing_tint_rgba,
        ground_tint_opacity=user_options.ground_tint_opacity,
    )

def _render_image(
    *,
    image_size: tuple[int, int],
    scene: RenderSceneData,
    style: RenderStyle,
    compositor: SkyCompositorCache,
    draw_direction_grid: bool = False,
    search_overlay_target: SearchJumpTarget | None = None,
) -> QImage:
    width, height = image_size
    image = QImage(width, height, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    painter.setRenderHint(QPainter.Antialiasing)
    painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
    try:
        geometry = render_geometry.get_screen_geometry(
            width,
            height,
            scene.viewer.view_alt_deg,
            edge_fov_deg=scene.viewer.edge_fov_deg,
            content_fov_deg=scene.viewer.content_fov_deg,
        )
        frame = FrameContext(
            viewer=scene.viewer,
            time_obj=scene.time_obj,
            geometry=geometry,
            viewport_rect=QRect(0, 0, width, height),
            sky_update_interval=DEFAULT_EXPORT_IMAGE_SKY_UPDATE_INTERVAL,
        )
        label_candidates: list[dict[str, object]] = []
        host().render_base_scene_into_painter(
            painter,
            frame=frame,
            scene=scene,
            style=style,
            hud=RenderHudState(
                mouse_pos=QPoint(),
                overlay_info_bottom_left=False,
                time_of_day_marker_bottom_left=False,
                viewport_interaction_mode=False,
                viewport_interaction_stars=None,
                status_message=(
                    "Forecast: Open-Meteo"
                    if float(style.precipitation_opacity) > 0.0
                    else None
                ),
            ),
            compositor=compositor,
            label_candidates=label_candidates,
            draw_labels=False,
        )
        if draw_direction_grid:
            render_guides.draw_direction_grid_overlay(
                painter,
                frame.geometry,
                frame.viewer,
                (width, height),
            )
        if search_overlay_target is not None:
            draw_search_target_overlay(
                painter,
                frame.geometry,
                search_overlay_target,
                viewer_data=frame.viewer,
                theme=style.theme,
                text_font=style.text_font,
                draw_marker=True,
                marker_scale=compute_star_render_upscale_factor(
                    frame.geometry.radius * 2,
                    style.star_render_expected_width,
                ),
                label_candidates=label_candidates,
            )
        if label_candidates:
            render_text._draw_label_candidates(
                painter, label_candidates, style.text_font
            )
    finally:
        painter.end()
    return image

def _build_render_style(
    *,
    text_font: QFont,
    status_line_font: QFont,
    catalogs: PreparedWindowCatalogs,
    user_options: SkyWindowUserOptions,
    runtime_options: SkyWindowRuntimeOptions,
    theme: ThemeStyle,
) -> RenderStyle:
    show_dso = catalogs.dso_catalog_np is not None
    if user_options.show_dso_initial is not None:
        show_dso = (
            bool(user_options.show_dso_initial) and catalogs.dso_catalog_np is not None
        )
    show_asterisms = (
        True
        if user_options.show_asterisms_initial is None
        else bool(user_options.show_asterisms_initial)
    )
    show_guidelines = (
        True
        if user_options.show_guidelines_initial is None
        else bool(user_options.show_guidelines_initial)
    )
    return RenderStyle(
        theme=theme,
        visual_preset=user_options.visual_preset,
        text_font=text_font,
        status_line_font=status_line_font,
        show_background_gradient=False,
        show_custom_window_frame=False,
        show_observation_info=False,
        show_dso=show_dso,
        show_asterisms=show_asterisms,
        show_guidelines=show_guidelines,
        moon_style=str(user_options.moon_style),
        moon_scale=int(user_options.moon_scale),
        bright_bodies_mode=str(user_options.bright_bodies_mode),
        star_base_radius=float(user_options.star_base_radius),
        star_visibility_boost=float(user_options.star_visibility_boost),
        sky_disc_alpha=float(user_options.sky_disc_alpha),
        asterism_visibility_boost=float(user_options.asterism_visibility_boost),
        earth_guide_visibility_boost=float(user_options.earth_guide_visibility_boost),
        vmag_limit=float(user_options.vmag_limit),
        sky_disc_altaz_rings=str(user_options.sky_disc_altaz_rings),
        sky_disc_altaz_rings_hover=str(user_options.sky_disc_altaz_rings_hover),
        cloud_disc_alpha=float(user_options.cloud_disc_alpha),
        satellite_opacity=float(user_options.satellite_opacity),
        terrain_horizon_opacity=float(user_options.terrain_horizon_opacity),
        earth_guide_opacity=float(user_options.earth_guide_opacity),
        night_light_opacity=float(user_options.night_light_opacity),
        diffuse_sky_source=str(user_options.diffuse_sky_source),
        road_night_lights_opacity=float(user_options.road_light_opacity),
        precipitation_opacity=float(user_options.precipitation_opacity),
        akari_ir_bands_opacity=float(user_options.akari_ir_bands_opacity),
        urban_outline_opacity=float(user_options.urban_outline_opacity),
        show_urban_outline_layer=float(user_options.urban_outline_opacity) > 0.0,
        aircraft_opacity=float(user_options.aircraft_opacity),
        star_render_expected_width=int(runtime_options.star_render_expected_width),
    )
