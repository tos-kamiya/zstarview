from __future__ import annotations

from collections.abc import Callable
from typing import cast

from PySide6.QtCore import QRect, QSize, Qt
from PySide6.QtGui import QImage, QPainter

from ..render.pipeline import FrameContext, RenderHudState
from ..types import CelestialData, ScreenGeometry, ViewerData
from .display_tone_curve import (
    DISPLAY_TONE_CURVE_LUT_VERSION,
    apply_display_tone_curve,
)

_FAST_FRAME_MAX_EDGE_PX = 600


class SkyWindowRenderCacheMixin:
    def _display_frame_image(self, present_frame: QImage) -> QImage:
        curve = getattr(self, "display_tone_curve", None)
        if curve is None:
            self._display_frame_cache_key = None
            self._display_frame_cache_image = None
            return present_frame
        key = (
            "display-frame",
            int(present_frame.cacheKey()),
            tuple(curve),
            DISPLAY_TONE_CURVE_LUT_VERSION,
        )
        cached = cast(
            QImage | None, getattr(self, "_display_frame_cache_image", None)
        )
        if getattr(self, "_display_frame_cache_key", None) != key or cached is None:
            cached = apply_display_tone_curve(present_frame, curve)
            self._display_frame_cache_key = key
            self._display_frame_cache_image = cached
        return cached

    def _render_cloud_state(self):
        return self._active_cloud_state()

    def _render_cache_stamp(self, value: object) -> object:
        if value is None:
            return None
        if hasattr(value, "cacheKey"):
            return int(value.cacheKey())
        return id(value)

    def _tropical_cyclone_snapshot_cache_value(self) -> object:
        state = self.tropical_cyclone_state
        if state is None:
            return None
        snapshots = state.snapshots
        if snapshots:
            return snapshots
        snapshot_collection = state.snapshot_collection
        if snapshot_collection is not None:
            return snapshot_collection
        return None

    def _render_frame_cache_key(
        self,
        *,
        frame: FrameContext | None = None,
        geometry: ScreenGeometry | None = None,
        celestial_data: CelestialData,
        render_viewer: ViewerData | None = None,
        include_fast_overlays: bool = True,
    ) -> tuple[object, ...]:
        if frame is None:
            if geometry is None or render_viewer is None:
                raise TypeError("frame or geometry/render_viewer must be provided")
            viewport_rect = QRect(
                0, 0, int(self.client_width()), int(self.client_height())
            )
            frame = FrameContext(
                viewer=render_viewer,
                time_obj=celestial_data.time,
                geometry=geometry,
                viewport_rect=viewport_rect,
                sky_update_interval=int(self.sky_update_interval),
            )
        cyclone_time_bucket = None
        try:
            current_time_obj = frame.time_obj
            if current_time_obj is None:
                current_time_obj = self._current_time_obj()
            cyclone_time_bucket = int(float(current_time_obj.unix) // 2.0)
        except Exception:
            cyclone_time_bucket = None
        key_parts: list[object] = [
            int(self.client_width()),
            int(self.client_height()),
            tuple(frame.geometry.center),
            int(frame.geometry.radius),
            self.visual_preset,
            self.presentation_id,
            bool(self.state.viewport_interaction_mode),
            tuple(float(v) for v in self.state.render_view_center),
            tuple(float(v) for v in frame.viewer.location),
            frame.viewer.city_name,
            frame.viewer.timezone_name,
            float(frame.viewer.observer_height_m),
            float(frame.viewer.height_add_m),
            frame.viewer.ground_elevation_m,
            float(frame.viewer.content_fov_deg),
            frame.viewer.location_height_label,
            frame.viewer.location_height_m,
            bool(self.show_dso),
            bool(self.show_asterisms),
            bool(self.show_guidelines),
            bool(self.show_observation_info),
            bool(self.state.simplified_view_enabled),
            bool(self.state.simplified_view_labels_enabled),
            bool(self.show_tropical_cyclone_overlay),
            round(float(self.tropical_cyclone_opacity), 3),
            self.bright_bodies_mode,
            str(self.moon_style),
            int(self.moon_scale),
            round(float(self.vmag_limit), 3),
            round(float(self.sky_disc_alpha), 3),
            self.sky_disc_style,
            round(float(self.cloud_disc_alpha), 3),
            round(float(self.terrain_horizon_opacity), 3),
            round(float(self.earth_guide_opacity), 3),
            round(float(self.night_light_opacity), 3),
            round(float(self.akari_ir_bands_opacity), 3),
            round(float(self.ridge_glow_opacity), 3),
            round(float(self.urban_outline_opacity), 3),
            bool(self.show_urban_outline_layer),
            getattr(
                self.state,
                "current_display_mode",
                "inverted-city"
                if bool(self.inverted_city_enabled)
                else "normal",
            ),
            bool(self.inverted_city_enabled),
            self._render_cache_stamp(celestial_data),
            self._render_cache_stamp(self.state.sky_disc_image),
            self._render_cache_stamp(self.state.night_light_glow_profile),
            self._render_cache_stamp(self._render_cloud_state().image),
            self._render_cache_stamp(self._render_cloud_state().missing_mask),
            None
            if self._render_cloud_state().cloud_amount_field is None
            else int(self._render_cloud_state().cloud_amount_field.source_cache_key),
            self._render_cache_stamp(self.state.terrain_horizon_profile),
            self._render_cache_stamp(self.state.terrain_horizon_profile_distances_m),
            self._render_cache_stamp(self.state.terrain_secondary_ridges_altaz_layers),
            self._render_cache_stamp(
                self.state.terrain_secondary_ridges_distances_m_layers
            ),
            self._render_cache_stamp(self.state.urban_outlines),
            self._render_cache_stamp(self.state.water_overlay_dots),
        ]
        if include_fast_overlays:
            cyclone_state = self.tropical_cyclone_state
            overlay_time_bucket = cyclone_time_bucket
            satellite_overlay_source = self.satellite_state.records_by_group
            aircraft_overlay_source = self.aircraft_state.snapshots
            key_parts.extend(
                [
                    round(float(self.satellite_opacity), 3),
                    round(float(self.aircraft_opacity), 3),
                    overlay_time_bucket,
                    self._render_cache_stamp(satellite_overlay_source),
                    self._render_cache_stamp(aircraft_overlay_source),
                    self._render_cache_stamp(
                        SkyWindowRenderCacheMixin._tropical_cyclone_snapshot_cache_value(
                            self
                        )
                    ),
                    cyclone_state.banner_text if cyclone_state is not None else None,
                ]
            )
        return tuple(key_parts)

    def _render_cached_frame_image(
        self,
        *,
        frame_key: tuple[object, ...],
        render_fn: Callable[[QPainter], None],
        cache_key_attr: str,
        cache_image_attr: str,
    ) -> QImage:
        return SkyWindowRenderCacheMixin._render_cached_image(
            self,
            image_size=self.client_size(),
            frame_key=frame_key,
            render_fn=render_fn,
            cache_key_attr=cache_key_attr,
            cache_image_attr=cache_image_attr,
        )

    def _render_cached_image(
        self,
        *,
        image_size: QSize,
        frame_key: tuple[object, ...],
        render_fn: Callable[[QPainter], None],
        cache_key_attr: str,
        cache_image_attr: str,
    ) -> QImage:
        frame_cache_key = self.__dict__.get(cache_key_attr)
        frame_cache_image = cast(QImage | None, self.__dict__.get(cache_image_attr))
        if frame_cache_key != frame_key or frame_cache_image is None:
            frame = QImage(
                image_size,
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
            self.__dict__[cache_image_attr] = frame
            self.__dict__[cache_key_attr] = frame_key
            return frame
        return cast(QImage, frame_cache_image)

    def _fast_frame_render_size(self, frame: FrameContext) -> QSize:
        win_w = max(1, int(frame.viewport_rect.width()))
        win_h = max(1, int(frame.viewport_rect.height()))
        max_edge = max(1, int(_FAST_FRAME_MAX_EDGE_PX))
        scale = min(1.0, float(max_edge) / float(max(win_w, win_h)))
        return QSize(
            max(1, int(round(win_w * scale))),
            max(1, int(round(win_h * scale))),
        )

    def _present_frame_cache_key(
        self,
        *,
        base_frame_key: tuple[object, ...],
        hud: RenderHudState,
    ) -> tuple[object, ...]:
        overlay_time_bucket = None
        star_interpolation_bucket = None
        try:
            current_time_obj = self._current_time_obj()
            overlay_time_bucket = int(float(current_time_obj.unix) // 2.0)
            celestial_data = self.state.celestial_data
            if (
                str(self.presentation_id).strip().lower()
                == "scenic"
                and celestial_data is not None
                and (celestial_data.star_time or celestial_data.time) is not None
            ):
                star_time = celestial_data.star_time or celestial_data.time
                elapsed_seconds = float(current_time_obj.unix - star_time.unix)
                interval_seconds = max(
                    1.0, float(self.sky_update_interval)
                )
                half_interval_seconds = interval_seconds / 2.0
                bucket_width_seconds = interval_seconds / 6.0
                star_interpolation_bucket = int(
                    max(
                        0.0,
                        min(interval_seconds, elapsed_seconds + half_interval_seconds),
                    )
                    // bucket_width_seconds
                )
        except Exception:
            overlay_time_bucket = None
            star_interpolation_bucket = None
        return (
            "present-frame",
            base_frame_key,
            star_interpolation_bucket,
            self._render_cache_stamp(self.state.dynamic_planets),
            str(self.sky_disc_altaz_rings),
            str(self.sky_disc_altaz_rings_hover),
            round(float(self.satellite_opacity), 3),
            round(float(self.aircraft_opacity), 3),
            bool(self.show_tropical_cyclone_overlay),
            round(float(self.tropical_cyclone_opacity), 3),
            overlay_time_bucket,
            self._render_cache_stamp(self.satellite_state.records_by_group),
            self._render_cache_stamp(self.aircraft_state.snapshots),
            self._render_cache_stamp(
                SkyWindowRenderCacheMixin._tropical_cyclone_snapshot_cache_value(self)
            ),
            self.tropical_cyclone_state.banner_text,
        )
