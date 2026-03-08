from __future__ import annotations

import logging
import time
from typing import Any

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QImage, QPainter, QPaintEvent

from ..render import draw as render_draw
from ..types import CelestialData, ViewerData

logger = logging.getLogger(__name__)


class SkyWindowRenderMixin:
    def _viewer_data_for_render(self) -> ViewerData:
        return ViewerData(
            location=self.viewer_data.location,
            timezone_name=self.viewer_data.timezone_name,
            city_name=self.viewer_data.city_name,
            view_center=self.state.render_view_center,
            observer_height_m=self.viewer_data.observer_height_m,
        )

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
        best_idx = None
        best_vmag = float("inf")
        for idx, raw_name in enumerate(stars["name"]):
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

        self._clear_background_layer(painter)
        self._draw_background_layer(painter, geometry)
        self._draw_sky_cloud_layers(painter, geometry)
        label_reservations: list[QRectF] = []
        label_candidates: list[dict[str, Any]] = []
        self._draw_terrain_layers(
            painter,
            geometry,
            celestial_data,
            render_viewer,
            highlighted_object,
            highlighted_dso,
            label_reservations,
            label_candidates,
        )
        self._draw_star_layer(painter, geometry, celestial_data, render_viewer)

        enlarge_moon = self.enlarge_moon
        if highlighted_object is not None:
            obj = highlighted_object[0]
            name = getattr(obj, "name", "") if hasattr(obj, "name") else obj.get("name", "")
            enlarge_moon = enlarge_moon or name == "moon"

        self._draw_planet_layer(
            painter,
            geometry,
            celestial_data,
            render_viewer,
            enlarge_moon,
            label_candidates,
        )
        self._draw_overlay_layer(
            painter,
            geometry,
            celestial_data,
            render_viewer,
            highlighted_object,
            highlighted_dso,
            enlarge_moon,
            label_reservations,
            label_candidates,
        )
        self._draw_label_layer(painter, label_candidates)
        self._draw_status_line(painter)

    def _clear_background_layer(self, painter: QPainter) -> None:
        painter.save()
        painter.setCompositionMode(QPainter.CompositionMode_Clear)
        painter.fillRect(self.rect(), Qt.transparent)
        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)
        painter.restore()

    def _draw_background_layer(self, painter: QPainter, geometry: render_draw.ScreenGeometry) -> None:
        render_draw.draw_radial_background(
            painter, QRectF(self.rect()), geometry, preset=self.visual_preset
        )

    def _draw_sky_cloud_layers(self, painter: QPainter, geometry: render_draw.ScreenGeometry) -> None:
        self._compositor.draw(
            painter,
            geometry,
            self.state.sky_disc_image,
            self.cloud_state.image,
            cloud_alpha=self.cloud_disc_alpha,
            stripe_density=self.cloud_state.stripe_density,
            missing_mask=self.cloud_state.missing_mask,
        )

    def _draw_terrain_layers(
        self,
        painter: QPainter,
        geometry: render_draw.ScreenGeometry,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
        highlighted_object: Any | None,
        highlighted_dso: Any | None,
        label_reservations: list[QRectF],
        label_candidates: list[dict[str, Any]],
    ) -> None:
        if self.show_dso:
            render_draw.draw_deep_sky_shapes(
                painter,
                geometry,
                celestial_data,
                render_viewer,
                preset=self.visual_preset,
            )
            render_draw.draw_dso_hover_info(
                painter,
                geometry,
                render_viewer,
                highlighted_dso,
                self.text_font,
                preset=self.visual_preset,
            )
        if self.show_asterisms:
            render_draw.draw_asterisms(
                painter,
                geometry,
                celestial_data,
                render_viewer,
                highlighted_object,
                self.text_font,
                label_reservations,
                label_candidates=label_candidates,
                preset=self.visual_preset,
            )
        render_draw.draw_sky_reference_lines(painter, geometry, celestial_data)
        render_draw.draw_terrain_horizon_line(
            painter,
            geometry,
            self.state.terrain_horizon_profile,
            render_viewer.view_center,
            opacity=self.terrain_horizon_opacity,
        )
        render_draw.draw_direction_labels(
            painter,
            geometry,
            render_viewer.view_center,
            self.text_font,
            self.state.mouse_pos,
            preset=self.visual_preset,
        )
        render_draw.draw_zenith_marker(painter, geometry, render_viewer.view_center)

    def _draw_star_layer(
        self,
        painter: QPainter,
        geometry: render_draw.ScreenGeometry,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
    ) -> None:
        win_w = self.width()
        win_h = self.height()
        low_w, low_h = self.compute_star_render_surface_size(
            win_w,
            win_h,
            geometry.radius * 2,
            self._star_render_expected_width,
        )
        stats = (win_w, win_h, low_w, low_h)
        last_star_render_stats = self.state.last_star_render_stats
        if stats != last_star_render_stats:
            logger.info(
                "Star render resolution: window=%dx%d draw=%dx%d",
                win_w,
                win_h,
                low_w,
                low_h,
            )
            self.state.last_star_render_stats = stats

        if low_w == win_w and low_h == win_h:
            render_draw.draw_stars(
                painter,
                geometry,
                celestial_data,
                render_viewer,
                self.star_base_radius,
                visibility_boost=self.star_visibility_boost,
                draw_vmag_limit=self.vmag_limit,
                viewport_size=(win_w, win_h),
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
            celestial_data,
            render_viewer,
            self.star_base_radius,
            visibility_boost=self.star_visibility_boost,
            draw_vmag_limit=self.vmag_limit,
            viewport_size=(low_w, low_h),
        )
        low_painter.end()

        painter.save()
        painter.setRenderHint(QPainter.SmoothPixmapTransform, False)
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
        painter.drawImage(self.rect(), low_img)
        painter.restore()

    def _draw_planet_layer(
        self,
        painter: QPainter,
        geometry: render_draw.ScreenGeometry,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
        enlarge_moon: bool,
        label_candidates: list[dict[str, Any]],
    ) -> None:
        render_draw.draw_solar_system_bodies(
            painter,
            geometry,
            celestial_data,
            render_viewer,
            enlarge_moon,
            text_font=self.text_font,
            label_candidates=label_candidates,
            draw_labels=True,
            preset=self.visual_preset,
        )

    def _draw_overlay_layer(
        self,
        painter: QPainter,
        geometry: render_draw.ScreenGeometry,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
        highlighted_object: Any | None,
        highlighted_dso: Any | None,
        enlarge_moon: bool,
        label_reservations: list[QRectF],
        label_candidates: list[dict[str, Any]],
    ) -> None:
        render_draw.draw_overlay_info(
            painter,
            geometry,
            celestial_data,
            render_viewer,
            self.vmag_limit,
            enlarge_moon,
            highlighted_dso,
            highlighted_object,
            self.text_font,
            label_candidates=label_candidates,
            label_reservations=label_reservations,
            preset=self.visual_preset,
        )

    def _draw_label_layer(self, painter: QPainter, label_candidates: list[dict[str, Any]]) -> None:
        render_draw.draw_label_candidates(painter, label_candidates, self.text_font)

    def _draw_status_line(self, painter: QPainter) -> None:
        message = self._status_line_message()
        if message:
            render_draw.draw_status_line_text(
                painter=painter,
                message=message,
                status_line_font=self.status_line_font,
                viewport_rect=self.rect(),
                preset=self.visual_preset,
            )
