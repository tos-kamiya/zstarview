from __future__ import annotations

import math
from typing import Any

import numpy as np
from PySide6.QtCore import QPoint, QPointF, QRect, QRectF, Qt
from PySide6.QtGui import QColor, QImage, QPainter, QTransform

from ..gui.composite import SkyCompositorCache
from ..night_lights import akari_midnight_opacity, night_activity_factor
from ..paths import ThemeStyle
from ..satellites.types import SatelliteOverlayPoint
from ..search.models import SearchJumpTarget
from ..simplified_view import resolve_simplified_view_mode
from ..tropical_cyclones.models import TropicalCycloneSnapshot
from ..types import (
    CelestialData,
    CelestialObject,
    ScreenGeometry,
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


def _sun_alt_deg(celestial_data: CelestialData) -> float | None:
    sun_altaz = _sun_altaz(celestial_data)
    return None if sun_altaz is None else sun_altaz[0]


def _sun_altaz(celestial_data: CelestialData) -> tuple[float, float] | None:
    for body in celestial_data.planets:
        if body.name == "sun":
            return float(body.alt), float(body.az)
    return None


def scene_night_activity_factor(
    scene: RenderSceneData,
    *,
    time_obj: Any | None = None,
) -> float:
    """Return the shared local-clock factor for artificial-light layers."""
    current_time = time_obj if time_obj is not None else scene.time_obj
    return night_activity_factor(
        current_time,
        scene.viewer.timezone_name,
        sun_alt_deg=_sun_alt_deg(scene.celestial_data),
    )


def scene_akari_opacity_factor(
    scene: RenderSceneData,
    *,
    time_obj: Any | None = None,
    base_opacity: float = 0.10,
) -> float:
    """Return the AKARI opacity for the current scene time."""
    return akari_midnight_opacity(
        base_opacity,
        scene_night_activity_factor(scene, time_obj=time_obj),
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
        labels_enabled=bool(getattr(hud, "simplified_view_labels_enabled", True)),
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
            draw_asterisms=draw_asterisms,
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
) -> None:
    """Draw dynamic satellite, aircraft, and cyclone overlays and their labels."""
    if (
        style.satellite_opacity <= 0.0
        and style.aircraft_opacity <= 0.0
        and float(style.tropical_cyclone_opacity) <= 0.0
    ):
        return
    local_label_candidates = label_candidates if label_candidates is not None else []
    _draw_satellite_layer(
        painter,
        geometry=frame.geometry,
        scene=scene,
        style=style,
        highlighted_satellite=highlighted_satellite,
        draw_simplified_labels=draw_simplified_satellite_labels,
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
    for snapshot in scene.tropical_cyclone_snapshots or ():
        render_tropical_cyclones.draw_tropical_cyclone_overlay(
            painter,
            geometry=frame.geometry,
            viewer=scene.viewer,
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
            viewer=scene.viewer,
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
        style=style,
        mouse_pos=hud.mouse_pos,
        highlighted_object=highlighted_object,
        highlighted_dso=highlighted_dso,
        highlighted_satellite=highlighted_satellite,
        label_candidates=label_candidates,
        draw_simplified_satellite_labels=simplified_view_labels_visible,
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
        from .zstarview_pipeline import _star_interpolation_matrix

        painter.save()
        label_matrix = _star_interpolation_matrix(frame=frame, scene=scene)
        if label_matrix is not None:
            _set_painter_homography(painter, label_matrix)
        _draw_simplified_named_star_labels(
            painter,
            geometry=frame.geometry,
            viewport_rect=frame.viewport_rect,
            scene=scene,
            style=style,
            highlighted_object=highlighted_object,
        )
        painter.restore()
    if not simplified_view_active:
        _draw_static_observation_overlay(
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
    if hud.status_message:
        render_text._draw_status_line_text(
            painter=painter,
            message=hud.status_message,
            status_line_font=style.status_line_font,
            viewport_rect=frame.viewport_rect,
            theme=style.theme,
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
            scene.viewer,
            (int(viewport_rect.width()), int(viewport_rect.height())),
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
    render_guides.draw_celestial_pole_markers(
        painter,
        geometry,
        scene.viewer,
        theme=style.theme,
    )


def _draw_main_terrain_profile_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    line_width_scale: float,
    fast_mode: bool,
) -> None:
    if style.terrain_horizon_opacity <= 0.0:
        return
    render_terrain._draw_terrain_profile_layer(
        painter,
        geometry,
        scene.viewer,
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
    time_obj: Any | None = None,
) -> None:
    if not style.show_urban_outline_layer:
        return
    render_terrain.draw_urban_outlines(
        painter,
        geometry,
        scene.viewer,
        scene.urban_outlines,
        opacity=style.urban_outline_opacity,
        fill_opacity_factor=scene_night_activity_factor(scene, time_obj=time_obj),
        line_width_scale=1.0,
        layer_style=style.theme.overlays.urban_outline,
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
        viewer_data=scene.viewer,
        aircraft_snapshots=scene.aircraft_snapshots,
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
    draw_vmag_min_exclusive: float | None = None,
    fast_mode: bool = False,
    star_interpolation_matrix: np.ndarray | None = None,
    separate_bright_stars: bool = False,
    bright_stars_only: bool = False,
) -> None:
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
    content_fov_deg = float(scene.viewer.content_fov_deg)
    split_bright_stars = (
        separate_bright_stars
        and
        draw_vmag_limit is None
        and float(getattr(style, "sky_disc_alpha", 0.0)) > 0.0
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
    ) -> None:
        draw_stars = (
            render_stars.draw_stars_fast if fast_mode else render_stars.draw_stars_normal
        )
        draw_stars(
            target,
            pass_geometry,
            draw_data,
            scene.viewer,
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
            content_fov_deg=content_fov_deg,
        )

    def draw_bright_star_pass(target: QPainter) -> None:
        target.save()
        if star_interpolation_matrix is not None:
            _set_painter_homography(target, star_interpolation_matrix)
        try:
            render_stars.draw_bright_star_underlay(
                target,
                geometry,
                draw_data,
                scene.viewer,
                style.star_base_radius,
                outline_bright_bodies=outline_bright_bodies,
                outline_render_scale=outline_render_scale,
                viewport_size=(win_w, win_h),
                content_fov_deg=content_fov_deg,
            )
            draw_star_pass(
                target,
                geometry,
                (win_w, win_h),
                draw_vmag_limit_override=4.0,
            )
        finally:
            target.restore()

    if bright_stars_only:
        draw_bright_star_pass(painter)
        return

    if low_w == win_w and low_h == win_h:
        if split_bright_stars:
            # Keep the faint-star raster separate from the bright-star pass so
            # the latter can remain crisp when a future interpolation transform
            # is applied to the faint-star surface.
            painter.save()
            if star_interpolation_matrix is not None:
                _set_painter_homography(painter, star_interpolation_matrix)
            draw_star_pass(
                painter,
                geometry,
                (win_w, win_h),
                draw_vmag_min_exclusive=4.0,
            )
            painter.restore()
            draw_bright_star_pass(painter)
        else:
            if star_interpolation_matrix is None:
                draw_star_pass(
                    painter,
                    geometry,
                    (win_w, win_h),
                    draw_vmag_min_exclusive=draw_vmag_min_exclusive,
                )
            else:
                painter.save()
                _set_painter_homography(painter, star_interpolation_matrix)
                draw_star_pass(
                    painter,
                    geometry,
                    (win_w, win_h),
                    draw_vmag_min_exclusive=draw_vmag_min_exclusive,
                )
                painter.restore()
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
        )
    else:
        draw_star_pass(
            low_painter,
            low_geometry,
            (low_w, low_h),
            draw_vmag_min_exclusive=draw_vmag_min_exclusive,
        )
    low_painter.end()

    painter.save()
    painter.setRenderHint(QPainter.SmoothPixmapTransform, False)
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
    if star_interpolation_matrix is None:
        painter.drawImage(viewport_rect, low_img)
    else:
        sx = low_w / max(1.0, float(win_w))
        sy = low_h / max(1.0, float(win_h))
        low_matrix = np.asarray(star_interpolation_matrix, dtype=float) @ np.array(
            [[1.0 / sx, 0.0, 0.0], [0.0, 1.0 / sy, 0.0], [0.0, 0.0, 1.0]],
            dtype=float,
        )
        _set_painter_homography(painter, low_matrix)
        painter.drawImage(0, 0, low_img)
    painter.restore()

    if split_bright_stars:
        # Bright stars and their contrast underlay are composited at viewport
        # resolution. The underlay must remain SourceOver while the star pass
        # itself uses the renderer's additive composition mode.
        draw_bright_star_pass(painter)


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
    viewport_rect: QRect,
    star_interpolation_matrix: np.ndarray | None,
) -> None:
    """Composite a cached faint-star image at its interpolated position."""
    painter.save()
    painter.setRenderHint(QPainter.SmoothPixmapTransform, False)
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
    if star_interpolation_matrix is None:
        painter.drawImage(viewport_rect, image)
    else:
        _set_painter_homography(painter, star_interpolation_matrix)
        painter.drawImage(0, 0, image)
    painter.restore()


def _draw_planet_layer(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    scene: RenderSceneData,
    style: RenderStyle,
    enlarge_moon: bool,
    outline_bright_bodies: bool = False,
    dark_contrast_enabled: bool = False,
    label_candidates: list[dict[str, Any]],
    draw_labels: bool = True,
    interpolation_matrix: np.ndarray | None = None,
) -> None:
    marker_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    def draw_bodies() -> None:
        render_solar_system.draw_solar_system_bodies(
            painter,
            geometry,
            scene.viewer,
            scene.celestial_data,
            enlarge_moon,
            outline_bright_bodies=outline_bright_bodies,
            text_font=style.text_font,
            label_candidates=label_candidates,
            draw_labels=draw_labels,
            theme=style.theme,
            edge_fov_deg=float(scene.viewer.edge_fov_deg),
            content_fov_deg=float(scene.viewer.content_fov_deg),
            marker_scale=marker_scale,
            instrument_presentation=_is_instrument_presentation(style),
            dark_contrast_enabled=dark_contrast_enabled,
            planet_bodies=scene.dynamic_planets,
        )

    if interpolation_matrix is None:
        draw_bodies()
        return
    painter.save()
    _set_painter_homography(painter, interpolation_matrix)
    try:
        draw_bodies()
    finally:
        painter.restore()


def _draw_static_observation_overlay(
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
    draw_simplified_labels: bool = False,
) -> None:
    render_satellites.draw_satellite_overlay(
        painter,
        geometry,
        viewer_data=scene.viewer,
        satellite_records_by_group=scene.satellite_records_by_group,
        time_obj=scene.time_obj,
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
    style: RenderStyle,
    mouse_pos: QPoint | None = None,
    highlighted_object: tuple[CelestialObject, QPointF] | None,
    highlighted_dso: tuple[CelestialObject, QPointF] | None,
    highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None = None,
    label_candidates: list[dict[str, Any]] | None = None,
    draw_simplified_satellite_labels: bool = False,
) -> None:
    line_width_scale = compute_star_render_upscale_factor(
        geometry.radius * 2,
        style.star_render_expected_width,
    )
    if style.show_asterisms and highlighted_object is not None:
        render_asterisms.draw_asterisms(
            painter,
            geometry,
            scene.viewer,
            scene.celestial_data,
            highlighted_object,
            style.text_font,
            theme=style.theme,
            line_width_scale=line_width_scale,
            content_fov_deg=float(scene.viewer.content_fov_deg),
            draw_base=False,
            draw_highlight=True,
            label_candidates=label_candidates,
        )
    render_solar_system.draw_hovered_moon_overlay(
        painter,
        geometry,
        scene.viewer,
        scene.celestial_data,
        highlighted_object,
        marker_scale=line_width_scale,
        outline_bright_bodies=str(style.bright_bodies_mode) == "outline",
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
                content_fov_deg=float(scene.viewer.content_fov_deg),
                opaque=not style.show_custom_window_frame,
            )
            render_background.draw_altitude_ring_overlay(
                painter,
                QRectF(viewport_rect),
                geometry,
                view_center=scene.viewer.view_center,
                theme=style.theme,
                edge_fov_deg=float(scene.viewer.edge_fov_deg),
                content_fov_deg=float(scene.viewer.content_fov_deg),
                ring_color=render_background.dimalt_ring_pen_color_from_color(
                    dimalt_sample_color
                ),
            )
        elif style.sky_disc_altaz_rings_hover == "altaz":
            render_guides.draw_direction_grid_overlay(
                painter,
                geometry,
                scene.viewer,
                (int(viewport_rect.width()), int(viewport_rect.height())),
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
    style: RenderStyle,
    highlighted_object: tuple[CelestialObject, QPointF] | None,
) -> None:
    if scene.celestial_data is None:
        return
    star_positions = render_stars.collect_visible_named_star_labels(
        scene.celestial_data,
        scene.viewer,
        geometry,
        style.star_base_radius,
        outline_bright_bodies=str(style.bright_bodies_mode) == "outline",
        outline_render_scale=compute_star_render_upscale_factor(
            geometry.radius * 2,
            style.star_render_expected_width,
        ),
        draw_vmag_limit=style.vmag_limit,
        content_fov_deg=float(scene.viewer.content_fov_deg),
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
            float(star_pos.x()) - float(text_bounds.left()),
            float(star_pos.y()) - float(text_bounds.bottom()),
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
