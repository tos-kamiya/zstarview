from __future__ import annotations

import math
from typing import Any

import numpy as np
from PySide6.QtCore import QPoint, QPointF, QRect, QRectF, Qt
from PySide6.QtGui import QColor, QFontMetrics, QImage, QPainter, QPainterPath, QTransform

from ..gui.composite import SkyCompositorCache
from ..moon_hover import MoonHoverImage
from ..night_lights import (
    diffuse_sky_artificial_light_attenuation_factor,
    diffuse_sky_sun_altitude_factor,
    night_light_strength_factor,
    post_solar_midnight_activity_factor,
)
from ..paths import ThemeStyle
from ..satellites.types import SatelliteOverlayPoint
from ..search.models import SearchJumpTarget
from ..simplified_view import resolve_simplified_view_mode
from ..solar_hover import SolarHoverImage
from ..tropical_cyclones.models import TropicalCycloneSnapshot
from ..types import (
    CelestialData,
    CelestialObject,
    ScreenGeometry,
    ViewerData,
)
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
from . import tropical_cyclones as render_tropical_cyclones
from .render_types import (
    FrameContext,
    RenderHudState,
    RenderSceneData,
    RenderStyle,
)

URBAN_OUTLINE_FILL_FACTOR_FLOOR = 0.05
STAR_DISC_CLIP_OVERSCAN_DEG = 0.75


def _sun_alt_deg(celestial_data: CelestialData) -> float | None:
    sun_altaz = _sun_altaz(celestial_data)
    return None if sun_altaz is None else sun_altaz[0]


def _sun_altaz(celestial_data: CelestialData) -> tuple[float, float] | None:
    for body in celestial_data.planets:
        if body.name == "sun":
            return float(body.alt), float(body.az)
    return None


def scene_diffuse_sky_opacity_factor(
    scene: RenderSceneData,
    viewer: ViewerData,
    *,
    base_opacity: float = 0.10,
    artificial_lights_enabled: bool = True,
) -> float:
    """Return diffuse-sky opacity for the scene's current Sun altitude."""
    sun_alt_deg = _sun_alt_deg(scene.celestial_data)
    if sun_alt_deg is None:
        return max(0.0, float(base_opacity))
    opacity = max(0.0, float(base_opacity)) * diffuse_sky_sun_altitude_factor(
        sun_alt_deg
    )
    if not artificial_lights_enabled:
        return opacity
    return opacity * diffuse_sky_artificial_light_attenuation_factor(
        scene_post_solar_midnight_activity_factor(scene, viewer)
    )


def scene_post_solar_midnight_activity_factor(
    scene: RenderSceneData, viewer: ViewerData
) -> float:
    """Return the solar-geometry human-activity factor for the current scene."""
    sun_altaz = _sun_altaz(scene.celestial_data)
    if sun_altaz is None:
        return 1.0
    return post_solar_midnight_activity_factor(
        sun_altaz[0],
        sun_altaz[1],
        float(viewer.location[0]),
    )


def scene_urban_outline_fill_factor(
    scene: RenderSceneData, viewer: ViewerData
) -> float:
    """Return the solar-altitude factor for illuminated building roofs."""
    sun_alt_deg = _sun_alt_deg(scene.celestial_data)
    if sun_alt_deg is None:
        return URBAN_OUTLINE_FILL_FACTOR_FLOOR
    solar_activity_factor = night_light_strength_factor(
        sun_alt_deg
    ) * scene_post_solar_midnight_activity_factor(scene, viewer)
    return max(URBAN_OUTLINE_FILL_FACTOR_FLOOR, solar_activity_factor)


def scene_night_light_opacity_factor(
    scene: RenderSceneData, viewer: ViewerData
) -> float:
    """Return roof-fill activity with its five-percent floor removed."""
    return max(
        0.0,
        scene_urban_outline_fill_factor(scene, viewer) - URBAN_OUTLINE_FILL_FACTOR_FLOOR,
    )


def _ground_reset_rgba_for_theme(theme: ThemeStyle) -> tuple[int, int, int, int]:
    """Return the below-horizon reset fill tuned for the active theme."""
    red, green, blue = (18, 18, 18)
    alpha = int(round(255.0 * float(theme.sky_disc.opacity)))
    return (
        int(red),
        int(green),
        int(blue),
        int(np.clip(alpha, 0, 255)),
    )


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


def _should_draw_water_overlay(scene: RenderSceneData, style: RenderStyle) -> bool:
    return float(style.terrain_horizon_opacity) > 0.0 and scene.terrain_horizon_profile is not None


def _simplified_view_active(hud: RenderHudState) -> bool:
    return _effective_simplified_view_mode(hud) != "normal"


def _effective_simplified_view_mode(hud: RenderHudState) -> str:
    return resolve_simplified_view_mode(
        base_enabled=bool(hud.simplified_view_enabled),
        labels_enabled=bool(hud.simplified_view_labels_enabled),
    )


def _simplified_view_labels_visible(hud: RenderHudState) -> bool:
    return _effective_simplified_view_mode(hud) == "labels"


def _is_instrument_presentation(style: RenderStyle) -> bool:
    return str(style.presentation_id).strip().lower() == "instrument"


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
    if _is_instrument_presentation(style):
        from .atlas_pipeline import InstrumentSkyPresentation

        InstrumentSkyPresentation().render_base_scene_into_painter(
            painter,
            frame=frame,
            scene=scene,
            style=style,
            hud=hud,
            compositor=compositor,
            draw_fast_overlays=draw_fast_overlays,
            label_candidates=label_candidates,
            draw_labels=draw_labels,
                draw_direction_labels=draw_direction_labels,
                draw_stars=draw_stars,
                draw_planets=draw_planets,
            )
        return
    from . import zstarview_pipeline

    zstarview_pipeline.render_base_scene_into_painter(
        painter,
        frame=frame,
        scene=scene,
        style=style,
        hud=hud,
        compositor=compositor,
        draw_fast_overlays=draw_fast_overlays,
        label_candidates=label_candidates,
        draw_labels=draw_labels,
        draw_direction_labels=draw_direction_labels,
        draw_stars=draw_stars,
        draw_planets=draw_planets,
        draw_asterisms=draw_asterisms,
    )


def render_fast_overlay_layers_into_painter(
    painter: QPainter,
    *,
    frame: FrameContext,
    scene: RenderSceneData,
    style: RenderStyle,
    highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None = None,
    highlighted_tropical_cyclone: TropicalCycloneSnapshot | None = None,
    label_candidates: list[dict[str, Any]] | None = None,
    draw_labels: bool = True,
    draw_simplified_satellite_labels: bool = False,
    fast_mode: bool = False,
) -> None:
    """Draw dynamic satellite, aircraft, meteor, and cyclone overlays."""
    meteor_opacity = float(style.meteor_opacity)
    if (
        style.satellite_opacity <= 0.0
        and style.aircraft_opacity <= 0.0
        and (fast_mode or meteor_opacity <= 0.0)
        and float(style.tropical_cyclone_opacity) <= 0.0
    ):
        return
    local_label_candidates = label_candidates if label_candidates is not None else []
    _draw_satellite_layer(
        painter,
        geometry=frame.geometry,
        scene=scene,
        viewer=frame.viewer,
        style=style,
        highlighted_satellite=highlighted_satellite,
        draw_simplified_labels=draw_simplified_satellite_labels,
        time_obj=frame.time_obj,
    )
    _draw_aircraft_layer(
        painter,
        geometry=frame.geometry,
        scene=scene,
        viewer=frame.viewer,
        style=style,
        label_candidates=local_label_candidates,
        time_obj=frame.time_obj,
    )
    if meteor_opacity > 0.0 and not fast_mode:
        from . import meteors as render_meteors

        meteor_core_color = (
            render_meteors.METEOR_ATLAS_CORE_COLOR
            if _is_instrument_presentation(style)
            else render_meteors.METEOR_CORE_COLOR
        )
        render_meteors.draw_meteor_trails(
            painter,
            frame.geometry,
            viewer_data=frame.viewer,
            trails=scene.meteor_trails,
            time_obj=frame.time_obj,
            opacity=meteor_opacity,
            core_color=meteor_core_color,
        )
    if draw_labels:
        render_text._draw_label_candidates(painter, local_label_candidates, style.text_font)
    for snapshot in scene.tropical_cyclone_snapshots or ():
        render_tropical_cyclones.draw_tropical_cyclone_overlay(
            painter,
            geometry=frame.geometry,
            viewer=frame.viewer,
            snapshot=snapshot,
            when_utc=frame.time_obj.to_datetime() if frame.time_obj is not None else None,
            theme=style.theme,
            opacity=float(style.tropical_cyclone_opacity),
            highlighted=bool(
                highlighted_tropical_cyclone is not None
                and snapshot == highlighted_tropical_cyclone
            ),
            enabled=bool(style.show_tropical_cyclone_overlay and style.tropical_cyclone_opacity > 0.0),
        )


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
    highlighted_tropical_cyclone: tuple[TropicalCycloneSnapshot, QPointF] | None = None,
    external_moon_image: MoonHoverImage | None = None,
    external_solar_image: SolarHoverImage | None = None,
    label_candidates: list[dict[str, Any]] | None = None,
    search_overlay_target: SearchJumpTarget | None = None,
) -> None:
    if hud.viewport_interaction_mode:
        if hud.status_message:
            render_text._draw_status_line_text(
                painter=painter,
                message=hud.status_message,
                status_line_font=style.status_line_font,
                viewport_rect=frame.viewport_rect,
                theme=style.theme,
            )
        if hud.mode_status_message:
            render_text._draw_mode_status_line_text(
                painter=painter,
                message=hud.mode_status_message,
                status_line_font=style.status_line_font,
                viewport_rect=frame.viewport_rect,
                theme=style.theme,
                below_message=hud.status_message,
            )
        return

    simplified_view_active = _simplified_view_active(hud)
    simplified_view_labels_visible = _simplified_view_labels_visible(hud)
    render_satellites.draw_satellite_highlight_overlay(
        painter,
        highlighted_satellite,
        opacity=float(style.satellite_opacity),
        marker_scale=compute_star_render_upscale_factor(
            frame.geometry.radius * 2,
            style.star_render_expected_width,
        ),
        theme=style.theme,
    )
    if highlighted_tropical_cyclone is not None:
        render_tropical_cyclones.draw_tropical_cyclone_overlay(
            painter,
            geometry=frame.geometry,
            viewer=frame.viewer,
            snapshot=highlighted_tropical_cyclone[0],
            when_utc=frame.time_obj.to_datetime() if frame.time_obj is not None else None,
            theme=style.theme,
            opacity=float(style.tropical_cyclone_opacity),
            highlighted=True,
            enabled=bool(style.show_tropical_cyclone_overlay and style.tropical_cyclone_opacity > 0.0),
        )
    _draw_hover_overlay_layer(
        painter,
        geometry=frame.geometry,
        viewport_rect=frame.viewport_rect,
        scene=scene,
        viewer=frame.viewer,
        style=style,
        mouse_pos=hud.mouse_pos,
        highlighted_object=highlighted_object,
        highlighted_dso=highlighted_dso,
        highlighted_satellite=highlighted_satellite,
        external_moon_image=external_moon_image,
        external_solar_image=external_solar_image,
        label_candidates=label_candidates,
        draw_simplified_satellite_labels=simplified_view_labels_visible,
        time_obj=frame.time_obj,
    )
    if search_overlay_target is not None:
        render_search_overlay.draw_search_target_overlay(
            painter,
            frame.geometry,
            search_overlay_target,
            viewer_data=frame.viewer,
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
    hide_simplified_instrument_labels = (
        _is_instrument_presentation(style)
        and simplified_view_active
        and not simplified_view_labels_visible
    )
    if label_candidates and not hide_simplified_instrument_labels:
        render_text._draw_label_candidates(painter, label_candidates, style.text_font)
    if (
        simplified_view_active
        and simplified_view_labels_visible
        and not _is_instrument_presentation(style)
    ):
        _draw_simplified_named_star_labels(
            painter,
            geometry=frame.geometry,
            viewport_rect=frame.viewport_rect,
            scene=scene,
            viewer=frame.viewer,
            style=style,
            highlighted_object=highlighted_object,
        )
    if not simplified_view_active:
        _draw_static_observation_overlay(
            painter,
            geometry=frame.geometry,
            viewport_rect=frame.viewport_rect,
            scene=scene,
            viewer=frame.viewer,
            style=style,
            overlay_info_bottom_left=hud.overlay_info_bottom_left,
            highlighted_object=None,
            highlighted_dso=None,
            label_reservations=[],
            label_candidates=label_candidates,
            status_message=hud.status_message,
            mode_status_message=hud.mode_status_message,
        )
    if hud.status_message:
        render_text._draw_status_line_text(
            painter=painter,
            message=hud.status_message,
            status_line_font=style.status_line_font,
            viewport_rect=frame.viewport_rect,
            theme=style.theme,
        )
    if hud.mode_status_message:
        render_text._draw_mode_status_line_text(
            painter=painter,
            message=hud.mode_status_message,
            status_line_font=style.status_line_font,
            viewport_rect=frame.viewport_rect,
            theme=style.theme,
            below_message=hud.status_message,
        )


def _draw_guide_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewport_rect: QRect,
    viewer: ViewerData,
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
            viewer,
            (int(viewport_rect.width()), int(viewport_rect.height())),
        )
    if draw_direction_labels:
        render_guides.draw_direction_labels(
            painter,
            geometry,
            viewer,
            style.text_font,
            None,
            theme=style.theme,
        )
    render_guides.draw_zenith_marker(
        painter,
        geometry,
        viewer,
        theme=style.theme,
    )
    render_guides.draw_celestial_pole_markers(
        painter,
        geometry,
        viewer,
        theme=style.theme,
    )


def _draw_main_terrain_profile_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    viewer: ViewerData,
    style: RenderStyle,
    line_width_scale: float,
    fast_mode: bool,
) -> None:
    if style.terrain_horizon_opacity <= 0.0:
        return
    render_terrain._draw_terrain_profile_layer(
        painter,
        geometry,
        viewer,
        scene.terrain_horizon_profile,
        scene.terrain_horizon_profile_distances_m,
        spec=render_terrain.TerrainHorizonRenderSpec(
            opacity=style.terrain_horizon_opacity,
            base_width=render_terrain.TERRAIN_HORIZON_FAST_WIDTH,
            far_base_width=render_terrain.TERRAIN_HORIZON_FAR_BASE_WIDTH,
            fg_alpha=render_terrain.terrain_horizon_line_alpha(
                style.terrain_horizon_opacity
                * style.theme.overlays.terrain_horizon.alpha_scale
            ),
            line_width_scale=(
                line_width_scale * style.theme.overlays.terrain_horizon.width_scale
            ),
            color_rgb=style.theme.overlays.terrain_horizon.rgb,
            fast_mode=fast_mode,
            distance_widths=True,
        ),
        is_in_fov_func=render_terrain.is_in_fov,
        altaz_to_normalized_xy_func=render_terrain.altaz_to_normalized_xy,
        normalized_to_screen_xy_func=render_terrain.normalized_to_screen_xy,
        split_by_gaps_func=render_terrain.split_by_gaps,
    )


def _draw_dso_hover_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    style: RenderStyle,
    highlighted_dso: tuple[CelestialObject, QPointF] | None,
) -> None:
    if not style.show_dso:
        return
    render_deep_sky_objects.draw_dso_hover_info(
        painter,
        geometry,
        viewer,
        highlighted_dso,
        style.text_font,
        theme=style.theme,
    )


def _draw_urban_outline_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    viewer: ViewerData,
    style: RenderStyle,
) -> None:
    if not style.show_urban_outline_layer:
        return
    render_terrain.draw_urban_outlines(
        painter,
        geometry,
        viewer,
        scene.urban_outlines,
        opacity=style.urban_outline_opacity,
        fill_opacity_factor=scene_urban_outline_fill_factor(scene, viewer),
        line_width_scale=1.0,
        layer_style=style.theme.overlays.urban_outline,
        inverted_city=bool(style.inverted_city_enabled),
    )


def _draw_aircraft_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    viewer: ViewerData,
    style: RenderStyle,
    label_candidates: list[dict[str, Any]],
    time_obj: Any | None,
) -> None:
    line_width_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    render_aircraft.draw_aircraft_overlay(
        painter,
        geometry,
        viewer_data=viewer,
        aircraft_snapshots=scene.aircraft_snapshots,
        time_obj=time_obj,
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
    viewer: ViewerData,
    style: RenderStyle,
    star_render_surface_size: tuple[int, int] | None = None,
    draw_vmag_limit: float | None = None,
    draw_vmag_min_exclusive: float | None = None,
    fast_mode: bool = False,
    separate_bright_stars: bool = False,
    bright_stars_only: bool = False,
    twinkle_targets: tuple[tuple[int, float], ...] = (),
    clip_to_disc: bool = True,
    render_cache: render_stars.StarRenderCache | None = None,
) -> None:
    if fast_mode:
        twinkle_targets = ()
    draw_data = scene.celestial_data
    win_w, win_h = int(viewport_rect.width()), int(viewport_rect.height())
    outline_render_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    outline_bright_bodies = str(style.bright_bodies_mode) == "outline"
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
    content_fov_deg = float(viewer.content_fov_deg)
    split_bright_stars = (
        separate_bright_stars
        and
        draw_vmag_limit is None
        and float(style.sky_disc_alpha) > 0.0
        and not style.light_background_star_outline
        and hasattr(painter, "save")
    )

    def draw_star_pass(
        target: QPainter,
        pass_geometry: ScreenGeometry,
        pass_size: tuple[int, int],
        *,
        draw_vmag_min_exclusive: float | None = None,
        draw_vmag_limit_override: float | None = None,
        screen_positions: np.ndarray | None = None,
        apply_disc_clip: bool = clip_to_disc,
    ) -> None:
        draw_stars = (
            render_stars.draw_stars_fast if fast_mode else render_stars.draw_stars_normal
        )
        target_save = getattr(target, "save", None)
        target_restore = getattr(target, "restore", None)
        if apply_disc_clip and callable(target_save) and callable(target_restore):
            target_save()
            _set_star_disc_clip(
                target,
                pass_geometry,
                edge_fov_deg=float(viewer.edge_fov_deg),
                content_fov_deg=content_fov_deg,
            )
        try:
            draw_stars(
                target,
                pass_geometry,
                draw_data,
                viewer,
                style.star_base_radius,
                visibility_boost=style.star_visibility_boost,
                outline_bright_bodies=outline_bright_bodies,
                outline_render_scale=outline_render_scale,
                light_background_outline=style.light_background_star_outline,
                draw_vmag_limit=(
                    draw_vmag_limit_override
                    if draw_vmag_limit_override is not None
                    else (draw_vmag_limit if draw_vmag_limit is not None else style.vmag_limit)
                ),
                draw_vmag_min_exclusive=draw_vmag_min_exclusive,
                viewport_size=pass_size,
                twinkle_targets=twinkle_targets,
                screen_positions=screen_positions,
                render_cache=render_cache,
            )
        finally:
            if apply_disc_clip and callable(target_save) and callable(target_restore):
                target_restore()

    def draw_bright_star_pass(target: QPainter) -> None:
        if clip_to_disc:
            target.save()
            _set_star_disc_clip(
                target,
                geometry,
                edge_fov_deg=float(viewer.edge_fov_deg),
                content_fov_deg=content_fov_deg,
            )
        try:
            _draw_bright_star_pass(target)
        finally:
            if clip_to_disc:
                target.restore()

    def _draw_bright_star_pass(target: QPainter) -> None:
        render_stars.draw_bright_star_underlay(
            target, geometry, draw_data, viewer, style.star_base_radius,
            outline_bright_bodies=outline_bright_bodies,
            outline_render_scale=outline_render_scale,
        )
        draw_star_pass(target, geometry, (win_w, win_h), draw_vmag_limit_override=4.0)

    if bright_stars_only:
        draw_bright_star_pass(painter)
        return

    if low_w == win_w and low_h == win_h:
        if split_bright_stars:
            draw_star_pass(
                painter,
                geometry,
                (win_w, win_h),
                draw_vmag_min_exclusive=4.0,
            )
            draw_bright_star_pass(painter)
        else:
            draw_star_pass(
                painter,
                geometry,
                (win_w, win_h),
                draw_vmag_min_exclusive=draw_vmag_min_exclusive,
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
    if split_bright_stars:
        # Only faint stars belong to the transformed/downsampled surface.
        draw_star_pass(
            low_painter,
            low_geometry,
            (low_w, low_h),
            draw_vmag_min_exclusive=4.0,
            apply_disc_clip=False,
        )
    else:
        draw_star_pass(
            low_painter,
            low_geometry,
            (low_w, low_h),
            draw_vmag_min_exclusive=draw_vmag_min_exclusive,
            apply_disc_clip=False,
        )
    low_painter.end()

    painter.save()
    if clip_to_disc:
        _set_star_disc_clip(
            painter,
            geometry,
            edge_fov_deg=float(viewer.edge_fov_deg),
            content_fov_deg=content_fov_deg,
        )
    painter.setRenderHint(QPainter.SmoothPixmapTransform, False)
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
    painter.drawImage(viewport_rect, low_img)
    painter.restore()

    if split_bright_stars:
        # Bright stars and their contrast underlay are composited at viewport
        # resolution. The underlay must remain SourceOver while the star pass
        # itself uses the renderer's additive composition mode.
        draw_bright_star_pass(painter)


def _set_star_disc_clip(
    painter: QPainter,
    geometry: ScreenGeometry,
    *,
    edge_fov_deg: float,
    content_fov_deg: float,
) -> None:
    radius = float(geometry.radius) * max(
        0.0,
        (float(content_fov_deg) + STAR_DISC_CLIP_OVERSCAN_DEG)
        / max(1.0e-6, float(edge_fov_deg)),
    )
    path = QPainterPath()
    path.addEllipse(QPointF(float(geometry.center[0]), float(geometry.center[1])), radius, radius)
    painter.setClipPath(path, Qt.ClipOperation.IntersectClip)


def _set_painter_homography(painter: QPainter, matrix: np.ndarray) -> None:
    """Set a row-vector NumPy homography on a Qt painter."""
    h = np.asarray(matrix, dtype=float)
    transform = QTransform()
    transform.setMatrix(
        float(h[0, 0]),
        float(h[1, 0]),
        float(h[2, 0]),
        float(h[0, 1]),
        float(h[1, 1]),
        float(h[2, 1]),
        float(h[0, 2]),
        float(h[1, 2]),
        float(h[2, 2]),
    )
    painter.setWorldTransform(transform, True)


def _draw_transformed_star_surface(
    painter: QPainter,
    image: QImage,
    *,
    geometry: ScreenGeometry,
    edge_fov_deg: float,
    content_fov_deg: float,
    viewport_rect: QRect,
) -> None:
    """Composite a cached faint-star image at its interpolated position."""
    painter.save()
    _set_star_disc_clip(
        painter,
        geometry,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
    )
    painter.setRenderHint(QPainter.SmoothPixmapTransform, False)
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
    painter.drawImage(viewport_rect, image)
    painter.restore()


def _draw_twinkle_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    viewer: ViewerData,
    style: RenderStyle,
    twinkle_targets: tuple[tuple[int, float], ...],
    fast_mode: bool = False,
) -> None:
    """Draw transient twinkle masks without invalidating the cached star surface."""
    if fast_mode or not twinkle_targets:
        return
    render_stars.draw_twinkle_overlay(
        painter,
        geometry,
        scene.celestial_data,
        viewer,
        style.star_base_radius,
        twinkle_targets=twinkle_targets,
    )


def _draw_planet_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    viewer: ViewerData,
    style: RenderStyle,
    outline_bright_bodies: bool = False,
    dark_contrast_enabled: bool = False,
    label_candidates: list[dict[str, Any]],
    draw_labels: bool = True,
    draw_markers: bool = True,
    suppress_moon_marker: bool = False,
    external_moon_image: MoonHoverImage | None = None,
) -> None:
    marker_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    def draw_bodies() -> None:
        render_solar_system.draw_solar_system_bodies(
            painter,
            geometry,
            viewer,
            scene.celestial_data,
            outline_bright_bodies=outline_bright_bodies,
            text_font=style.text_font,
            label_candidates=label_candidates,
            draw_labels=draw_labels,
            theme=style.theme,
            edge_fov_deg=float(viewer.edge_fov_deg),
            content_fov_deg=float(viewer.content_fov_deg),
            marker_scale=marker_scale,
            instrument_presentation=_is_instrument_presentation(style),
            dark_contrast_enabled=dark_contrast_enabled,
            planet_bodies=scene.dynamic_planets,
            draw_markers=draw_markers,
            suppress_moon_marker=suppress_moon_marker,
            moon_style=str(style.moon_style),
            moon_scale=int(style.moon_scale),
            external_moon_image=external_moon_image,
        )

    # Solar-system positions are supplied for the display tick when
    # available.  Never apply the star snapshot transform to them: their
    # position is independent of the star-surface interpolation path.
    draw_bodies()


def _draw_static_observation_overlay(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    viewer: ViewerData,
    style: RenderStyle,
    overlay_info_bottom_left: bool,
    highlighted_object: tuple[CelestialObject, QPointF] | None,
    highlighted_dso: tuple[CelestialObject, QPointF] | None,
    label_reservations: list[QRectF],
    label_candidates: list[dict[str, Any]],
    status_message: str | None = None,
    mode_status_message: str | None = None,
) -> None:
    if not style.show_observation_info:
        return
    status_line_count = (
        len(render_text.status_line_text_lines(status_message))
        if status_message
        else 0
    )
    if mode_status_message:
        status_line_count += len(render_text.status_line_text_lines(mode_status_message))
    bottom_reserved_height = float(
        QFontMetrics(style.status_line_font).lineSpacing() * status_line_count
    )
    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        scene.celestial_data,
        viewer,
        vmag_limit=style.vmag_limit,
        highlighted_dso=highlighted_dso,
        highlighted_object=highlighted_object,
        text_font=style.text_font,
        label_candidates=label_candidates,
        label_reservations=label_reservations,
        viewport_rect=viewport_rect,
        bottom_left=overlay_info_bottom_left,
        bottom_reserved_height=(
            bottom_reserved_height if overlay_info_bottom_left else 0.0
        ),
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
    viewer: ViewerData,
    style: RenderStyle,
    highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None,
    draw_simplified_labels: bool = False,
    time_obj: Any | None = None,
) -> None:
    render_satellites.draw_satellite_overlay(
        painter,
        geometry,
        viewer_data=viewer,
        satellite_records_by_group=scene.satellite_records_by_group,
        time_obj=time_obj,
        opacity=style.satellite_opacity,
        highlighted_satellite=(
            highlighted_satellite[0] if highlighted_satellite is not None else None
        ),
        marker_scale=compute_star_render_upscale_factor(
            geometry.radius * 2,
            style.star_render_expected_width,
        ),
        draw_simplified_labels=draw_simplified_labels,
        text_font=style.text_font,
        theme=style.theme,
    )


def _draw_hover_overlay_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    viewer: ViewerData,
    style: RenderStyle,
    mouse_pos: QPoint | None = None,
    highlighted_object: tuple[CelestialObject, QPointF] | None,
    highlighted_dso: tuple[CelestialObject, QPointF] | None,
    highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None = None,
    label_candidates: list[dict[str, Any]] | None = None,
    external_moon_image: MoonHoverImage | None = None,
    external_solar_image: SolarHoverImage | None = None,
    draw_simplified_satellite_labels: bool = False,
    time_obj: Any | None = None,
) -> None:
    line_width_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    if (
        style.show_asterisms
        and highlighted_object is not None
    ):
        render_asterisms.draw_asterisms(
            painter,
            geometry,
            viewer,
            scene.celestial_data,
            highlighted_object,
            style.text_font,
            theme=style.theme,
            line_width_scale=line_width_scale,
            content_fov_deg=float(viewer.content_fov_deg),
            draw_base=False,
            draw_highlight=True,
            base_line_alpha_scale=float(style.asterism_visibility_boost),
            opacity=style.asterism_opacity,
            label_candidates=label_candidates,
        )
    if str(style.moon_style) != "image":
        render_solar_system.draw_hovered_moon_overlay(
            painter,
            geometry,
            viewer,
            scene.celestial_data,
            highlighted_object,
            marker_scale=line_width_scale,
            theme=style.theme,
            external_moon_image=external_moon_image,
        )
    render_solar_system.draw_hovered_sun_overlay(
        painter,
        geometry,
        viewer,
        scene.celestial_data,
        highlighted_object,
        time_obj=time_obj,
        marker_scale=line_width_scale,
        text_font=style.text_font,
        theme=style.theme,
        external_solar_image=external_solar_image,
    )
    _draw_dso_hover_layer(
        painter,
        geometry=geometry,
        viewer=viewer,
        style=style,
        highlighted_dso=highlighted_dso,
    )
    direction_hover = None
    if style.show_guidelines and mouse_pos is not None:
        direction_hover = render_guides.resolve_direction_marker_hover(
            geometry,
            viewer,
            mouse_pos,
        )
    if direction_hover is not None:
        if style.sky_disc_altaz_rings_hover == "dimalt":
            dimalt_sample_color = render_background.sample_background_disc_edge_color(
                QRectF(viewport_rect),
                geometry,
                theme=style.theme,
                edge_fov_deg=float(viewer.edge_fov_deg),
                content_fov_deg=float(viewer.content_fov_deg),
                opaque=not style.show_custom_window_frame,
            )
            render_background.draw_altitude_ring_overlay(
                painter,
                geometry,
                view_center=viewer.view_center,
                edge_fov_deg=float(viewer.edge_fov_deg),
                content_fov_deg=float(viewer.content_fov_deg),
                ring_color=render_background.dimalt_ring_pen_color_from_color(
                    dimalt_sample_color
                ),
            )
        elif style.sky_disc_altaz_rings_hover == "altaz":
            render_guides.draw_direction_grid_overlay(
                painter,
                geometry,
                viewer,
                (int(viewport_rect.width()), int(viewport_rect.height())),
            )
    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        scene.celestial_data,
        viewer,
        vmag_limit=style.vmag_limit,
        highlighted_dso=highlighted_dso,
        highlighted_object=highlighted_object,
        text_font=style.text_font,
        highlighted_satellite=highlighted_satellite,
        label_candidates=label_candidates,
        theme=style.theme,
        draw_static_info=False,
        draw_hover_info=True,
        draw_simplified_satellite_labels=draw_simplified_satellite_labels,
        draw_outlined_text_func=render_text.draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
    )

def _draw_simplified_named_star_labels(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewport_rect: QRect,
    scene: RenderSceneData,
    viewer: ViewerData,
    style: RenderStyle,
    highlighted_object: tuple[CelestialObject, QPointF] | None,
) -> None:
    if scene.celestial_data is None:
        return
    star_positions = render_stars.collect_visible_named_star_labels(
        scene.celestial_data,
        viewer,
        geometry,
        style.star_base_radius,
        outline_bright_bodies=str(style.bright_bodies_mode) == "outline",
        outline_render_scale=compute_star_render_upscale_factor(
            geometry.radius * 2,
            style.star_render_expected_width,
        ),
        draw_vmag_limit=style.vmag_limit,
        viewport_size=(int(viewport_rect.width()), int(viewport_rect.height())),
    )
    if not star_positions:
        return

    highlighted_pos = highlighted_object[1] if highlighted_object is not None else None
    for star_name, star_pos, star_rgb in star_positions:
        if highlighted_pos is not None:
            if abs(float(star_pos.x()) - float(highlighted_pos.x())) < 1e-6 and abs(
                float(star_pos.y()) - float(highlighted_pos.y())
            ) < 1e-6:
                continue
        draw_pos = star_pos
        if _is_instrument_presentation(style):
            label_color = QColor(*style.theme.text.foreground_rgb[:3], 255)
        else:
            label_color = render_text.blend_color_toward_white(
                QColor(*star_rgb),
                amount=render_text.LABEL_COLOR_WHITE_BLEND_AMOUNT,
            )
            label_color.setAlpha(round(255.0 * 0.4))
        label_font = style.text_font
        text_bounds = render_text._text_bounds_at_baseline(star_name, label_font, QPointF(0.0, 0.0))
        label_pos = QPointF(
            float(draw_pos.x()) - float(text_bounds.left()),
            float(draw_pos.y()) - float(text_bounds.bottom()),
        )
        label_style = render_text.ResolvedTextStyle(
            font=label_font,
            text_color=label_color,
            outline_color=QColor(0, 0, 0, 0),
            outline_width=0.0,
        )
        render_text.draw_outlined_text(
            painter,
            star_name,
            label_pos,
            style=label_style,
        )
