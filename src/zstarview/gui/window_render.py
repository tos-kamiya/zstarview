from __future__ import annotations

import logging
import time
from typing import Any

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QImage, QPainter, QPaintEvent

from ..astro import resolve_star_names
from ..render import draw as render_draw
from ..render.pipeline import (
    RenderSceneData,
    RenderHudState,
    RenderStyle,
    compute_star_render_surface_size,
    render_base_scene_into_painter,
    render_hud_overlay_into_painter,
)
from ..types import CelestialData, ViewerData

logger = logging.getLogger(__name__)

class SkyWindowRenderMixin:
    def _render_cache_stamp(self, value: Any) -> Any:
        if value is None:
            return None
        if hasattr(value, "cacheKey"):
            return int(value.cacheKey())
        return id(value)

    def _render_frame_cache_key(
        self,
        *,
        geometry: render_draw.ScreenGeometry,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
    ) -> tuple[Any, ...]:
        return (
            int(self.width()),
            int(self.height()),
            tuple(geometry.center),
            int(geometry.radius),
            self.visual_preset,
            bool(self.state.viewport_interaction_mode),
            tuple(float(v) for v in self.state.render_view_center),
            tuple(float(v) for v in render_viewer.location),
            render_viewer.city_name,
            render_viewer.timezone_name,
            float(render_viewer.observer_height_m),
            float(render_viewer.content_fov_deg),
            render_viewer.location_height_label,
            render_viewer.location_height_m,
            bool(render_viewer.show_observer_height),
            bool(self.show_dso),
            bool(self.show_asterisms),
            bool(getattr(self, "show_guidelines", True)),
            bool(self.enlarge_moon),
            round(float(self.vmag_limit), 3),
            round(float(self.sky_disc_alpha), 3),
            round(float(self.cloud_disc_alpha), 3),
            round(float(getattr(self, "satellite_opacity", 0.0)), 3),
            round(float(getattr(self, "aircraft_opacity", 0.0)), 3),
            round(float(self.terrain_horizon_opacity), 3),
            round(float(self.urban_outline_opacity), 3),
            bool(getattr(self, "show_urban_outline_layer", True)),
            self._render_cache_stamp(celestial_data),
            self._render_cache_stamp(self.state.sky_disc_image),
            self._render_cache_stamp(self.cloud_state.image),
            self._render_cache_stamp(self.cloud_state.missing_mask),
            None if self.cloud_state.stripe_density is None else int(self.cloud_state.stripe_density.source_cache_key),
            self._render_cache_stamp(self.state.terrain_horizon_profile),
            self._render_cache_stamp(self.state.urban_outlines),
            self._render_cache_stamp(getattr(self.state, "satellite_overlay_points", None)),
            self._render_cache_stamp(self.state.aircraft_overlay_points),
        )

    def _draw_cached_frame(
        self,
        painter: QPainter,
        frame_key: tuple[Any, ...],
        render_fn: Any,
    ) -> None:
        if self._frame_cache_key != frame_key or self._frame_cache_image is None:
            frame = QImage(self.size(), QImage.Format.Format_ARGB32_Premultiplied)
            frame.fill(Qt.GlobalColor.transparent)
            frame_painter = QPainter(frame)
            frame_painter.setRenderHint(QPainter.Antialiasing)
            frame_painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
            try:
                render_fn(frame_painter)
            finally:
                frame_painter.end()
            self._frame_cache_image = frame
            self._frame_cache_key = frame_key
        painter.drawImage(0, 0, self._frame_cache_image)

    def _viewer_data_for_render(self) -> ViewerData:
        return ViewerData(
            location=self.viewer_data.location,
            timezone_name=self.viewer_data.timezone_name,
            city_name=self.viewer_data.city_name,
            view_center=self.state.render_view_center,
            content_fov_deg=self.viewer_data.content_fov_deg,
            observer_height_m=self.viewer_data.observer_height_m,
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
        cloud_state = getattr(self, "cloud_state", None)
        return RenderSceneData(
            viewer=render_viewer,
            celestial_data=celestial_data,
            sky_disc_image=getattr(state, "sky_disc_image", None),
            cloud_image=getattr(cloud_state, "image", None),
            cloud_missing_mask=getattr(cloud_state, "missing_mask", None),
            cloud_stripe_density=getattr(cloud_state, "stripe_density", None),
            terrain_horizon_profile=getattr(state, "terrain_horizon_profile", None),
            urban_outlines=getattr(state, "urban_outlines", None),
            satellite_overlay_points=getattr(state, "satellite_overlay_points", None),
            aircraft_overlay_points=getattr(state, "aircraft_overlay_points", None),
        )

    def _render_style(self) -> RenderStyle:
        return RenderStyle(
            visual_preset=self.visual_preset,
            text_font=self.text_font,
            status_line_font=getattr(self, "status_line_font", self.text_font),
            show_background_gradient=True,
            show_overlay_info=True,
            show_dso=bool(getattr(self, "show_dso", False)),
            show_asterisms=bool(getattr(self, "show_asterisms", False)),
            show_guidelines=bool(getattr(self, "show_guidelines", True)),
            enlarge_moon=bool(getattr(self, "enlarge_moon", False)),
            star_base_radius=float(getattr(self, "star_base_radius", 1.0)),
            star_visibility_boost=float(getattr(self, "star_visibility_boost", 1.0)),
            vmag_limit=float(getattr(self, "vmag_limit", 6.0)),
            cloud_disc_alpha=float(getattr(self, "cloud_disc_alpha", 0.0)),
            satellite_opacity=float(getattr(self, "satellite_opacity", 0.0)),
            terrain_horizon_opacity=float(getattr(self, "terrain_horizon_opacity", 0.0)),
            urban_outline_opacity=float(getattr(self, "urban_outline_opacity", 0.2)),
            show_urban_outline_layer=bool(getattr(self, "show_urban_outline_layer", True)),
            aircraft_opacity=float(getattr(self, "aircraft_opacity", 1.0)),
            star_render_expected_width=int(self._star_render_expected_width),
        )

    def _render_hud_state(self) -> RenderHudState:
        status_message = None
        if hasattr(self, "_status_line_message"):
            status_message = self._status_line_message()
        return RenderHudState(
            mouse_pos=self.state.mouse_pos,
            viewport_interaction_mode=bool(self.state.viewport_interaction_mode),
            viewport_interaction_stars=self.state.viewport_interaction_stars,
            status_message=status_message,
        )

    def _update_star_render_stats(self, geometry: render_draw.ScreenGeometry) -> None:
        win_w = int(self.width())
        win_h = int(self.height())
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

        nx, ny = render_draw.altaz_to_normalized_xy(
            alt,
            az,
            self.state.render_view_center,
        )
        px, py = render_draw.normalized_to_screen_xy(nx, ny, geometry)
        return ({"name": target_name}, QPointF(px, py))

    def render_current_image(self, *, include_hud: bool = False) -> QImage:
        """Render the current window state into an off-screen image."""
        image = QImage(self.size(), QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(Qt.GlobalColor.transparent)
        painter = QPainter(image)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        try:
            celestial_data = self.state.celestial_data
            if celestial_data is None:
                loading_color, _ = render_draw.get_text_style(self.visual_preset)
                painter.setPen(loading_color)
                painter.setFont(self.text_font)
                painter.drawText(
                    self.rect(),
                    Qt.AlignmentFlag.AlignCenter,
                    "Loading celestial data...",
                )
                return image

            render_viewer = self._viewer_data_for_render()
            geometry = render_draw.get_screen_geometry(
                self.width(),
                self.height(),
                render_viewer.view_center[0],
            )
            scene, style, hud = self._render_inputs(
                celestial_data=celestial_data,
                render_viewer=render_viewer,
            )
            render_base_scene_into_painter(
                painter,
                geometry=geometry,
                viewport_rect=self.rect(),
                scene=scene,
                style=style,
                hud=hud,
                compositor=self._compositor,
            )
            if include_hud:
                highlighted_object = None
                highlighted_dso = None
                jump_highlight = self._active_jump_highlight_object(geometry)
                if jump_highlight is not None:
                    highlighted_object = jump_highlight
                render_hud_overlay_into_painter(
                    painter,
                    geometry=geometry,
                    viewport_rect=self.rect(),
                    scene=scene,
                    style=style,
                    hud=hud,
                    highlighted_object=highlighted_object,
                    highlighted_dso=highlighted_dso,
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
            loading_color, _ = render_draw.get_text_style(self.visual_preset)
            painter.setPen(loading_color)
            painter.setFont(self.text_font)
            painter.drawText(
                self.rect(), Qt.AlignmentFlag.AlignCenter, "Loading celestial data..."
            )
            return

        render_viewer = self._viewer_data_for_render()
        alt = render_viewer.view_center[0]
        geometry = render_draw.get_screen_geometry(self.width(), self.height(), alt)

        highlighted_object = None
        highlighted_dso = None
        mouse_pos = self.state.mouse_pos
        if mouse_pos is not None:
            highlighted_object = render_draw.find_highlighted_object(
                celestial_data,
                render_viewer,
                mouse_pos,
                geometry,
            )
            if self.show_dso:
                highlighted_dso = render_draw.find_highlighted_dso(
                    celestial_data,
                    render_viewer,
                    mouse_pos,
                    geometry,
                )
        jump_highlight = self._active_jump_highlight_object(geometry)
        if jump_highlight is not None:
            highlighted_object = jump_highlight
        frame_key = self._render_frame_cache_key(
            geometry=geometry,
            celestial_data=celestial_data,
            render_viewer=render_viewer,
        )
        self._update_star_render_stats(geometry)
        scene, style, hud = self._render_inputs(
            celestial_data=celestial_data,
            render_viewer=render_viewer,
        )
        self._draw_cached_frame(
            painter,
            frame_key,
            lambda frame_painter: render_base_scene_into_painter(
                frame_painter,
                geometry=geometry,
                viewport_rect=self.rect(),
                scene=scene,
                style=style,
                hud=hud,
                compositor=self._compositor,
            ),
        )
        render_hud_overlay_into_painter(
            painter,
            geometry=geometry,
            viewport_rect=self.rect(),
            scene=scene,
            style=style,
            hud=hud,
            highlighted_object=highlighted_object,
            highlighted_dso=highlighted_dso,
        )
