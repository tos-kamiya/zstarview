from __future__ import annotations

import logging
import time
from typing import Callable, cast

from PySide6.QtCore import QPoint, QPointF, Qt
from PySide6.QtGui import QFont, QImage, QPainter, QPaintEvent

from ..astro import altaz_to_normalized_xy, resolve_star_names
from ..render import deep_sky_objects as render_deep_sky_objects
from ..render import geometry as render_geometry
from ..render import satellites as render_satellites
from ..render import stars as render_stars
from ..render import text as render_text
from ..render.search_overlay import draw_search_target_overlay
from ..render.pipeline import (
    RenderSceneData,
    RenderHudState,
    RenderStyle,
    compute_star_render_surface_size,
    compute_star_render_upscale_factor,
    render_base_scene_into_painter,
    render_fast_overlay_layers_into_painter,
    render_hud_overlay_into_painter,
    render_status_line_into_painter,
)
from ..satellites.types import SatelliteOverlayPoint
from ..types import CelestialData, CelestialObject, ScreenGeometry, ViewerData

logger = logging.getLogger(__name__)


def _resolve_hover_targets(
    *,
    celestial_data: CelestialData,
    render_viewer: ViewerData,
    mouse_pos: QPoint | None,
    geometry: ScreenGeometry,
    satellite_overlay_points: list[object] | None,
    show_dso: bool,
) -> tuple[
    tuple[CelestialObject, QPointF] | None,
    tuple[CelestialObject, QPointF] | None,
    tuple[SatelliteOverlayPoint, QPointF] | None,
]:
    highlighted_object = None
    highlighted_dso = None
    highlighted_satellite = None
    if mouse_pos is None:
        return highlighted_object, highlighted_dso, highlighted_satellite

    highlighted_object = render_stars.find_highlighted_object(
        celestial_data,
        render_viewer,
        mouse_pos,
        geometry,
    )
    if show_dso:
        highlighted_dso = render_deep_sky_objects.find_highlighted_dso(
            celestial_data,
            render_viewer,
            mouse_pos,
            geometry,
        )
    highlighted_satellite = render_satellites.find_highlighted_satellite(
        satellite_overlay_points,
        mouse_pos,
        geometry,
        render_viewer.view_center,
        edge_fov_deg=float(render_viewer.edge_fov_deg),
        content_fov_deg=render_viewer.content_fov_deg,
    )
    return highlighted_object, highlighted_dso, highlighted_satellite


class SkyWindowRenderMixin:
    def _render_cache_stamp(self, value: object) -> object:
        if value is None:
            return None
        if hasattr(value, "cacheKey"):
            return int(value.cacheKey())
        return id(value)

    def _render_frame_cache_key(
        self,
        *,
        geometry: ScreenGeometry,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
        include_fast_overlays: bool = True,
    ) -> tuple[object, ...]:
        key_parts: list[object] = [
            int(self.client_width()),
            int(self.client_height()),
            tuple(geometry.center),
            int(geometry.radius),
            self.visual_preset,
            bool(self.state.viewport_interaction_mode),
            tuple(float(v) for v in self.state.render_view_center),
            tuple(float(v) for v in render_viewer.location),
            render_viewer.city_name,
            render_viewer.timezone_name,
            float(render_viewer.observer_height_m),
            getattr(render_viewer, "ground_elevation_m", None),
            float(render_viewer.content_fov_deg),
            render_viewer.location_height_label,
            render_viewer.location_height_m,
            bool(render_viewer.show_observer_height),
            bool(self.show_dso),
            bool(self.show_asterisms),
            bool(self.show_guidelines),
            bool(self.show_observation_info),
            bool(self.enlarge_moon),
            getattr(self, "bright_bodies_mode", "outline"),
            round(float(self.vmag_limit), 3),
            round(float(self.sky_disc_alpha), 3),
            getattr(self, "sky_disc_style", "smooth"),
            round(float(self.cloud_disc_alpha), 3),
            round(float(self.terrain_horizon_opacity), 3),
            round(float(self.earth_guide_opacity), 3),
            round(float(getattr(self, "night_light_opacity", 0.0)), 3),
            round(float(self.urban_outline_opacity), 3),
            bool(self.show_urban_outline_layer),
            self._render_cache_stamp(celestial_data),
            self._render_cache_stamp(self.state.sky_disc_image),
            self._render_cache_stamp(self.state.night_light_glow_profile),
            self._render_cache_stamp(self.cloud_state.image),
            self._render_cache_stamp(self.cloud_state.missing_mask),
            None
            if self.cloud_state.cloud_amount_field is None
            else int(self.cloud_state.cloud_amount_field.source_cache_key),
            self._render_cache_stamp(self.state.terrain_horizon_profile),
            self._render_cache_stamp(self.state.terrain_horizon_profile_distances_m),
            self._render_cache_stamp(self.state.terrain_horizon_secondary_profile_altaz_layers),
            self._render_cache_stamp(self.state.terrain_horizon_secondary_profile_distances_m_layers),
            self._render_cache_stamp(self.state.urban_outlines),
        ]
        if include_fast_overlays:
            key_parts.extend(
                [
                    round(float(self.satellite_opacity), 3),
                    round(float(self.aircraft_opacity), 3),
                    None
                    if self.state.satellite_overlay_points is None
                    else tuple(self.state.satellite_overlay_points),
                    None
                    if self.state.aircraft_overlay_points is None
                    else tuple(self.state.aircraft_overlay_points),
                ]
            )
        return tuple(key_parts)

    def _draw_cached_frame(
        self,
        painter: QPainter,
        frame_key: tuple[object, ...],
        render_fn: Callable[[QPainter], None],
    ) -> None:
        frame_cache_image = SkyWindowRenderMixin._render_cached_frame_image(
            self,
            frame_key=frame_key,
            render_fn=render_fn,
            cache_key_attr="_frame_cache_key",
            cache_image_attr="_frame_cache_image",
        )
        painter.drawImage(0, 0, frame_cache_image)

    def _render_cached_frame_image(
        self,
        *,
        frame_key: tuple[object, ...],
        render_fn: Callable[[QPainter], None],
        cache_key_attr: str,
        cache_image_attr: str,
    ) -> QImage:
        frame_cache_key = getattr(self, cache_key_attr, None)
        frame_cache_image = cast(QImage | None, getattr(self, cache_image_attr, None))
        if frame_cache_key != frame_key or frame_cache_image is None:
            frame = QImage(
                self.client_size(),
                QImage.Format.Format_ARGB32_Premultiplied,
            )
            frame.fill(Qt.GlobalColor.transparent)
            frame_painter = QPainter(frame)
            frame_painter.setRenderHint(QPainter.Antialiasing)
            frame_painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
            try:
                render_fn(frame_painter)
            finally:
                frame_painter.end()
            setattr(self, cache_image_attr, frame)
            setattr(self, cache_key_attr, frame_key)
            return frame
        return cast(QImage, frame_cache_image)

    def _present_frame_cache_key(
        self,
        *,
        base_frame_key: tuple[object, ...],
        hud: RenderHudState,
    ) -> tuple[object, ...]:
        mouse_pos = hud.mouse_pos
        mouse_key = None
        if mouse_pos is not None:
            mouse_key = (int(mouse_pos.x()), int(mouse_pos.y()))
        jump_key = (
            self.state.jump_highlight_name,
            self.state.jump_highlight_altaz,
            round(float(self.state.jump_highlight_until_ms), 3),
        )
        search_target = getattr(self.state, "persistent_search_target", None)
        search_key = None
        if search_target is not None:
            search_key = (
                getattr(search_target, "label", None),
                getattr(search_target, "alt_deg", None),
                getattr(search_target, "az_deg", None),
                bool(getattr(search_target, "persistent_keep_marker", False)),
            )
        return (
            "present-frame",
            base_frame_key,
            str(getattr(self, "sky_disc_altaz_rings", "dimalt")),
            str(getattr(self, "sky_disc_altaz_rings_hover", "altaz")),
            round(float(self.satellite_opacity), 3),
            round(float(self.aircraft_opacity), 3),
            self._render_cache_stamp(self.state.satellite_overlay_points),
            self._render_cache_stamp(self.state.aircraft_overlay_points),
            mouse_key,
            bool(hud.overlay_info_bottom_left),
            bool(hud.viewport_interaction_mode),
            jump_key,
            search_key,
        )

    def _render_present_frame_image(
        self,
        *,
        base_frame_key: tuple[object, ...],
        geometry: ScreenGeometry,
        scene: RenderSceneData,
        style: RenderStyle,
        hud: RenderHudState,
        highlighted_object: tuple[CelestialObject, QPointF] | None,
        highlighted_dso: tuple[CelestialObject, QPointF] | None,
        highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None,
    ) -> QImage:
        base_label_candidates: list[dict[str, object]] = []
        base_frame_image = SkyWindowRenderMixin._render_cached_frame_image(
            self,
            frame_key=base_frame_key,
            render_fn=lambda frame_painter: (
                render_base_scene_into_painter(
                    frame_painter,
                    geometry=geometry,
                    viewport_rect=self.client_rect(),
                    scene=scene,
                    style=style,
                    hud=hud,
                    compositor=self._compositor,
                    draw_fast_overlays=False,
                    label_candidates=base_label_candidates,
                    draw_labels=False,
                ),
                setattr(
                    self,
                    "_cached_base_label_candidates",
                    list(base_label_candidates),
                ),
            ),
            cache_key_attr="_frame_cache_key",
            cache_image_attr="_frame_cache_image",
        )
        cached_base_label_candidates = getattr(self, "_cached_base_label_candidates", [])
        present_frame_key = SkyWindowRenderMixin._present_frame_cache_key(
            self,
            base_frame_key=base_frame_key,
            hud=hud,
        )
        return SkyWindowRenderMixin._render_cached_frame_image(
            self,
            frame_key=present_frame_key,
            render_fn=lambda frame_painter: (
                SkyWindowRenderMixin._draw_present_frame_layers(
                    self,
                    frame_painter=frame_painter,
                    base_frame_image=base_frame_image,
                    base_label_candidates=cached_base_label_candidates,
                    geometry=geometry,
                    scene=scene,
                    style=style,
                    hud=hud,
                    highlighted_object=highlighted_object,
                    highlighted_dso=highlighted_dso,
                    highlighted_satellite=highlighted_satellite,
                )
            ),
            cache_key_attr="_present_frame_cache_key",
            cache_image_attr="_present_frame_cache_image",
        )

    def _render_fast_frame_image(
        self,
        *,
        base_frame_key: tuple[object, ...],
        geometry: ScreenGeometry,
        scene: RenderSceneData,
        style: RenderStyle,
        hud: RenderHudState,
        highlighted_object: tuple[CelestialObject, QPointF] | None,
        highlighted_dso: tuple[CelestialObject, QPointF] | None,
        highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None,
    ) -> QImage:
        return self._render_present_frame_image(
            base_frame_key=base_frame_key,
            geometry=geometry,
            scene=scene,
            style=style,
            hud=hud,
            highlighted_object=highlighted_object,
            highlighted_dso=highlighted_dso,
            highlighted_satellite=highlighted_satellite,
        )

    def _render_normal_frame_image(
        self,
        *,
        base_frame_key: tuple[object, ...],
        geometry: ScreenGeometry,
        scene: RenderSceneData,
        style: RenderStyle,
        hud: RenderHudState,
        highlighted_object: tuple[CelestialObject, QPointF] | None,
        highlighted_dso: tuple[CelestialObject, QPointF] | None,
        highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None,
    ) -> QImage:
        return self._render_present_frame_image(
            base_frame_key=base_frame_key,
            geometry=geometry,
            scene=scene,
            style=style,
            hud=hud,
            highlighted_object=highlighted_object,
            highlighted_dso=highlighted_dso,
            highlighted_satellite=highlighted_satellite,
        )

    def _draw_present_frame_layers(
        self,
        *,
        frame_painter: QPainter,
        base_frame_image: QImage,
        base_label_candidates: list[dict[str, object]] | tuple[dict[str, object], ...] | None,
        geometry: ScreenGeometry,
        scene: RenderSceneData,
        style: RenderStyle,
        hud: RenderHudState,
        highlighted_object: tuple[CelestialObject, QPointF] | None,
        highlighted_dso: tuple[CelestialObject, QPointF] | None,
        highlighted_satellite: tuple[SatelliteOverlayPoint, QPointF] | None,
    ) -> None:
        frame_painter.drawImage(0, 0, base_frame_image)
        if hud.viewport_interaction_mode:
            render_status_line_into_painter(
                frame_painter,
                viewport_rect=self.client_rect(),
                style=style,
                hud=hud,
            )
            return

        label_candidates: list[dict[str, object]] = list(base_label_candidates or [])
        render_fast_overlay_layers_into_painter(
            frame_painter,
            geometry=geometry,
            scene=scene,
            style=style,
            highlighted_satellite=highlighted_satellite,
            label_candidates=label_candidates,
            draw_labels=False,
        )
        render_hud_overlay_into_painter(
            frame_painter,
            geometry=geometry,
            viewport_rect=self.client_rect(),
            scene=scene,
            style=style,
            hud=hud,
            highlighted_object=highlighted_object,
            highlighted_dso=highlighted_dso,
            highlighted_satellite=highlighted_satellite,
            label_candidates=label_candidates,
            search_overlay_target=getattr(self.state, "persistent_search_target", None),
        )

    def _viewer_data_for_render(self) -> ViewerData:
        return ViewerData(
            location=self.viewer_data.location,
            timezone_name=self.viewer_data.timezone_name,
            city_name=self.viewer_data.city_name,
            view_center=self.state.render_view_center,
            edge_fov_deg=self.viewer_data.edge_fov_deg,
            content_fov_deg=self.viewer_data.content_fov_deg,
            observer_height_m=self.viewer_data.observer_height_m,
            ground_elevation_m=self.viewer_data.ground_elevation_m,
            location_height_label=self.viewer_data.location_height_label,
            location_height_m=self.viewer_data.location_height_m,
            show_observer_height=self.viewer_data.show_observer_height,
        )

    def _render_inputs(
        self,
        *,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
    ) -> tuple[RenderSceneData, RenderStyle, RenderHudState]:
        return (
            SkyWindowRenderMixin._render_scene_data(
                self,
                celestial_data=celestial_data,
                render_viewer=render_viewer,
            ),
            SkyWindowRenderMixin._render_style(self),
            SkyWindowRenderMixin._render_hud_state(self),
        )

    def _render_scene_data(
        self,
        *,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
    ) -> RenderSceneData:
        state = self.state
        cloud_state = self.cloud_state
        return RenderSceneData(
            viewer=render_viewer,
            celestial_data=celestial_data,
            sky_disc_image=getattr(state, "sky_disc_image", None),
            cloud_image=getattr(cloud_state, "image", None),
            cloud_missing_mask=getattr(cloud_state, "missing_mask", None),
            cloud_amount_field=getattr(cloud_state, "cloud_amount_field", None),
            terrain_horizon_profile=getattr(state, "terrain_horizon_profile", None),
            terrain_horizon_profile_distances_m=getattr(
                state,
                "terrain_horizon_profile_distances_m",
                None,
            ),
            terrain_horizon_secondary_profile_altaz_layers=getattr(
                state,
                "terrain_horizon_secondary_profile_altaz_layers",
                None,
            ),
            terrain_horizon_secondary_profile_distances_m_layers=getattr(
                state,
                "terrain_horizon_secondary_profile_distances_m_layers",
                None,
            ),
            urban_outlines=getattr(state, "urban_outlines", None),
            water_overlay_points=getattr(state, "water_overlay_points", None),
            satellite_overlay_points=getattr(state, "satellite_overlay_points", None),
            aircraft_overlay_points=getattr(state, "aircraft_overlay_points", None),
            night_light_glow_profile=getattr(state, "night_light_glow_profile", None),
        )

    def _render_style(self) -> RenderStyle:
        status_line_font = self.status_line_font
        return RenderStyle(
            theme=self.theme,
            visual_preset=self.visual_preset,
            text_font=self.text_font,
            status_line_font=cast(QFont, status_line_font),
            show_background_gradient=True,
            show_custom_window_frame=bool(self._frameless_window),
            show_observation_info=bool(self.show_observation_info),
            show_dso=bool(self.show_dso),
            show_asterisms=bool(self.show_asterisms),
            show_guidelines=bool(self.show_guidelines),
            enlarge_moon=bool(self.enlarge_moon),
            bright_bodies_mode=str(getattr(self, "bright_bodies_mode", "outline")),
            star_base_radius=float(self.star_base_radius),
            star_visibility_boost=float(self.star_visibility_boost),
            asterism_visibility_boost=float(self.asterism_visibility_boost),
            earth_guide_visibility_boost=float(self.earth_guide_visibility_boost),
            vmag_limit=float(self.vmag_limit),
            sky_disc_altaz_rings=str(getattr(self, "sky_disc_altaz_rings", "dimalt")),
            sky_disc_altaz_rings_hover=str(getattr(self, "sky_disc_altaz_rings_hover", "altaz")),
            cloud_disc_alpha=float(self.cloud_disc_alpha),
            satellite_opacity=float(self.satellite_opacity),
            terrain_horizon_opacity=float(self.terrain_horizon_opacity),
            earth_guide_opacity=float(self.earth_guide_opacity),
            night_light_opacity=float(getattr(self, "night_light_opacity", 0.0)),
            urban_outline_opacity=float(self.urban_outline_opacity),
            show_urban_outline_layer=bool(self.show_urban_outline_layer),
            water_overlay_opacity=float(getattr(self, "water_overlay_opacity", 0.12)),
            aircraft_opacity=float(self.aircraft_opacity),
            star_render_expected_width=int(self._star_render_expected_width),
        )

    def _render_hud_state(self) -> RenderHudState:
        status_message = None
        if hasattr(self, "_status_line_message"):
            status_message = self._status_line_message()
        mouse_pos = self.state.mouse_pos
        if getattr(self, "_startup_input_blocked", lambda: False)():
            mouse_pos = None
        overlay_info_bottom_left = bool(
            getattr(self.state, "overlay_info_bottom_left", False)
        )
        if mouse_pos is not None:
            # Respect pinned CLI modes: when pinned, do not move the overlay based on mouse.
            if not self.observation_info_pinned:
                window_height = max(1, int(self.client_height()))
                upper_threshold = float(window_height) / 3.0
                lower_threshold = 2.0 * float(window_height) / 3.0
                mouse_y = float(mouse_pos.y())
                if mouse_y <= upper_threshold:
                    overlay_info_bottom_left = True
                elif mouse_y >= lower_threshold:
                    overlay_info_bottom_left = False
                self.state.overlay_info_bottom_left = overlay_info_bottom_left
            else:
                # Ensure the HUD uses the pinned position from state; do not override.
                overlay_info_bottom_left = bool(
                    getattr(self.state, "overlay_info_bottom_left", False)
                )
        return RenderHudState(
            mouse_pos=mouse_pos,
            overlay_info_bottom_left=overlay_info_bottom_left,
            viewport_interaction_mode=bool(self.state.viewport_interaction_mode),
            viewport_interaction_stars=self.state.viewport_interaction_stars,
            status_message=status_message,
        )

    def _update_star_render_stats(self, geometry: ScreenGeometry) -> None:
        win_w = int(self.client_width())
        win_h = int(self.client_height())
        low_w, low_h = compute_star_render_surface_size(
            win_w,
            win_h,
            geometry.radius * 2,
            int(self._star_render_expected_width),
        )
        stats = (win_w, win_h, low_w, low_h)
        if stats != self.state.last_star_render_stats:
            logger.info(
                "Star render resolution: window=%dx%d draw=%dx%d",
                win_w,
                win_h,
                low_w,
                low_h,
            )
            self.state.last_star_render_stats = stats

    def _active_jump_highlight_object(self, geometry):
        jump_highlight_name = self.state.jump_highlight_name
        if not jump_highlight_name:
            return None
        jump_highlight_until_ms = self.state.jump_highlight_until_ms
        if time.monotonic() * 1000.0 > jump_highlight_until_ms:
            self.state.jump_highlight_name = None
            self.state.jump_highlight_altaz = None
            self.state.jump_highlight_until_ms = 0.0
            return None
        celestial_data = self.state.celestial_data
        if not celestial_data:
            return None

        target_name = jump_highlight_name
        stars = celestial_data.stars
        star_names = resolve_star_names(stars, celestial_data.star_catalog_meta)
        best_idx = None
        best_vmag = float("inf")
        for idx, raw_name in enumerate(star_names):
            name = str(raw_name).strip()
            if name != target_name:
                continue
            vmag = float(stars["vmag"][idx])
            if vmag < best_vmag:
                best_vmag = vmag
                best_idx = idx

        if best_idx is not None:
            alt = float(stars["alt"][best_idx])
            az = float(stars["az"][best_idx])
        else:
            jump_highlight_altaz = self.state.jump_highlight_altaz
            if jump_highlight_altaz is not None:
                alt, az = jump_highlight_altaz
            else:
                return None

        nx, ny = altaz_to_normalized_xy(
            alt,
            az,
            self.state.render_view_center,
            edge_fov_deg=float(self.viewer_data.edge_fov_deg),
        )
        px, py = render_geometry.normalized_to_screen_xy(nx, ny, geometry)
        return ({"name": target_name}, QPointF(px, py))

    def _draw_persistent_search_overlay(
        self,
        painter: QPainter,
        geometry,
    ) -> None:
        target = getattr(self.state, "persistent_search_target", None)
        if target is None:
            return
        if not bool(getattr(target, "persistent_keep_marker", False)):
            return
        draw_search_target_overlay(
            painter,
            geometry,
            target,
            view_center=self.state.render_view_center,
            edge_fov_deg=float(self.viewer_data.edge_fov_deg),
            content_fov_deg=float(self.viewer_data.content_fov_deg),
            theme=self.theme,
            text_font=self.text_font,
            draw_marker=True,
            draw_label=True,
            marker_scale=compute_star_render_upscale_factor(
                geometry.radius * 2,
                int(self._star_render_expected_width),
            ),
        )

    def render_current_image(self, *, include_hud: bool = False) -> QImage:
        """Render the current window state into an off-screen image."""
        image = QImage(
            self.client_size(),
            QImage.Format.Format_ARGB32_Premultiplied,
        )
        image.fill(Qt.GlobalColor.transparent)
        painter = QPainter(image)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        try:
            celestial_data = self.state.celestial_data
            if celestial_data is None:
                loading_color, _ = render_text.get_text_style(self.theme)
                painter.setPen(loading_color)
                painter.setFont(self.text_font)
                painter.drawText(
                    self.client_rect(),
                    Qt.AlignmentFlag.AlignCenter,
                    "Loading celestial data...",
                )
                return image

            render_viewer = self._viewer_data_for_render()
            geometry = render_geometry.get_screen_geometry(
                int(self.client_width()),
                int(self.client_height()),
                render_viewer.view_center[0],
            )
            scene, style, hud = self._render_inputs(
                celestial_data=celestial_data,
                render_viewer=render_viewer,
            )
            if include_hud:
                label_candidates: list[dict[str, object]] = []
                render_base_scene_into_painter(
                    painter,
                    geometry=geometry,
                    viewport_rect=self.client_rect(),
                    scene=scene,
                    style=style,
                    hud=hud,
                    compositor=self._compositor,
                    label_candidates=label_candidates,
                    draw_labels=False,
                )
                highlighted_object = None
                highlighted_dso = None
                jump_highlight = self._active_jump_highlight_object(geometry)
                if jump_highlight is not None:
                    highlighted_object = jump_highlight
                render_hud_overlay_into_painter(
                    painter,
                    geometry=geometry,
                    viewport_rect=self.client_rect(),
                    scene=scene,
                    style=style,
                    hud=hud,
                    highlighted_object=highlighted_object,
                    highlighted_dso=highlighted_dso,
                    label_candidates=label_candidates,
                    search_overlay_target=getattr(self.state, "persistent_search_target", None),
                )
            else:
                render_base_scene_into_painter(
                    painter,
                    geometry=geometry,
                    viewport_rect=self.client_rect(),
                    scene=scene,
                    style=style,
                    hud=hud,
                    compositor=self._compositor,
                )
            return image
        finally:
            painter.end()

    def paintEvent(self, event: QPaintEvent) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)

        celestial_data = self.state.celestial_data
        if celestial_data is None:
            loading_color, _ = render_text.get_text_style(self.theme)
            painter.setPen(loading_color)
            painter.setFont(self.text_font)
            painter.drawText(
                self.client_rect(),
                Qt.AlignmentFlag.AlignCenter,
                "Loading celestial data...",
            )
            return

        render_viewer = self._viewer_data_for_render()
        alt = render_viewer.view_center[0]
        geometry = render_geometry.get_screen_geometry(
            int(self.client_width()),
            int(self.client_height()),
            alt,
        )
        frame_key = self._render_frame_cache_key(
            geometry=geometry,
            celestial_data=celestial_data,
            render_viewer=render_viewer,
            include_fast_overlays=False,
        )
        self._update_star_render_stats(geometry)
        scene, style, hud = self._render_inputs(
            celestial_data=celestial_data,
            render_viewer=render_viewer,
        )
        if self.state.viewport_interaction_mode:
            present_frame = self._render_fast_frame_image(
                base_frame_key=frame_key,
                geometry=geometry,
                scene=scene,
                style=style,
                hud=hud,
                highlighted_object=None,
                highlighted_dso=None,
                highlighted_satellite=None,
            )
        else:
            mouse_pos = self.state.mouse_pos
            if getattr(self, "_startup_input_blocked", lambda: False)():
                mouse_pos = None

            highlighted_object, highlighted_dso, highlighted_satellite = (
                _resolve_hover_targets(
                    celestial_data=celestial_data,
                    render_viewer=render_viewer,
                    mouse_pos=mouse_pos,
                    geometry=geometry,
                    satellite_overlay_points=self.state.satellite_overlay_points,
                    show_dso=bool(self.show_dso),
                )
            )
            jump_highlight = self._active_jump_highlight_object(geometry)
            if jump_highlight is not None:
                highlighted_object = jump_highlight
            present_frame = self._render_normal_frame_image(
                base_frame_key=frame_key,
                geometry=geometry,
                scene=scene,
                style=style,
                hud=hud,
                highlighted_object=highlighted_object,
                highlighted_dso=highlighted_dso,
                highlighted_satellite=highlighted_satellite,
            )
        painter.drawImage(0, 0, present_frame)
