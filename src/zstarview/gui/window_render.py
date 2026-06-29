from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from typing import Callable, cast

import astropy.time
from PySide6.QtCore import QPoint, QPointF, QRect, QSize, Qt
from PySide6.QtGui import QFont, QImage, QPainter, QPaintEvent

from ..astro import altaz_to_normalized_xy, resolve_star_names
from ..render import deep_sky_objects as render_deep_sky_objects
from ..render import geometry as render_geometry
from ..render import guides as render_guides
from ..render import satellites as render_satellites
from ..render import tropical_cyclones as render_tropical_cyclones
from ..render import stars as render_stars
from ..render import text as render_text
from ..render.pipeline import (
    FrameContext,
    RenderHudState,
    RenderSceneData,
    RenderStyle,
    _simplified_view_labels_visible,
    compute_star_render_surface_size,
    render_base_scene_into_painter,
    render_fast_overlay_layers_into_painter,
    render_hud_overlay_into_painter,
)
from ..satellites.types import SatelliteOverlayPoint
from ..tropical_cyclones.models import TropicalCycloneSnapshot
from ..types import CelestialData, CelestialObject, ScreenGeometry, ViewerData

logger = logging.getLogger(__name__)
_FAST_FRAME_MAX_EDGE_PX = 600


@dataclass(frozen=True, slots=True)
class HoverTargets:
    object: tuple[CelestialObject, QPointF] | None = None
    dso: tuple[CelestialObject, QPointF] | None = None
    satellite: tuple[SatelliteOverlayPoint, QPointF] | None = None
    tropical_cyclone: tuple[TropicalCycloneSnapshot, QPointF] | None = None


@dataclass(frozen=True, slots=True)
class RenderInputs:
    scene: RenderSceneData
    style: RenderStyle
    hud: RenderHudState


def _resolve_hover_targets(
    *,
    celestial_data: CelestialData,
    render_viewer: ViewerData,
    mouse_pos: QPoint | None,
    geometry: ScreenGeometry,
    satellite_records_by_group: object | None = None,
    tropical_cyclone_snapshots: object | None = None,
    time_obj: object | None = None,
    show_dso: bool = False,
) -> HoverTargets:
    highlighted_object = None
    highlighted_dso = None
    highlighted_satellite = None
    highlighted_tropical_cyclone = None
    if mouse_pos is None:
        return HoverTargets()

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
        satellite_records_by_group,
        mouse_pos,
        geometry,
        viewer_data=render_viewer,
        time_obj=time_obj,
    )
    highlighted_tropical_cyclone = (
        render_tropical_cyclones.find_highlighted_tropical_cyclone(
            tropical_cyclone_snapshots,
            mouse_pos,
            geometry,
            viewer_data=render_viewer,
            time_obj=time_obj.to_datetime()
            if hasattr(time_obj, "to_datetime")
            else None,
        )
    )
    return HoverTargets(
        object=highlighted_object,
        dso=highlighted_dso,
        satellite=highlighted_satellite,
        tropical_cyclone=highlighted_tropical_cyclone,
    )


class SkyWindowRenderMixin:
    def _startup_splash_visible(self) -> bool:
        overlay = self._owner._startup_log_overlay
        return bool(overlay is not None and overlay.isVisible())

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
            bool(getattr(self.state, "simplified_view_labels_enabled", True)),
            bool(self.show_tropical_cyclone_overlay),
            round(float(self.tropical_cyclone_opacity), 3),
            bool(self.enlarge_moon),
            self.bright_bodies_mode,
            round(float(self.vmag_limit), 3),
            round(float(self.sky_disc_alpha), 3),
            self.sky_disc_style,
            round(float(self.cloud_disc_alpha), 3),
            round(float(self.terrain_horizon_opacity), 3),
            round(float(self.earth_guide_opacity), 3),
            round(float(self.night_light_opacity), 3),
            round(float(self.ridge_glow_opacity), 3),
            round(float(self.urban_outline_opacity), 3),
            bool(self.show_urban_outline_layer),
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
                        SkyWindowRenderMixin._tropical_cyclone_snapshot_cache_value(
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
        return SkyWindowRenderMixin._render_cached_image(
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
        try:
            current_time_obj = self._current_time_obj()
            overlay_time_bucket = int(float(current_time_obj.unix) // 2.0)
        except Exception:
            overlay_time_bucket = None
        return (
            "present-frame",
            base_frame_key,
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
                SkyWindowRenderMixin._tropical_cyclone_snapshot_cache_value(self)
            ),
            self.tropical_cyclone_state.banner_text,
        )

    def _render_present_frame_image(
        self,
        *,
        base_frame_key: tuple[object, ...],
        frame: FrameContext,
        render_inputs: RenderInputs,
    ) -> QImage:
        base_label_candidates: list[dict[str, object]] = []
        base_frame_image = SkyWindowRenderMixin._render_cached_frame_image(
            self,
            frame_key=base_frame_key,
            render_fn=lambda frame_painter: (
                render_base_scene_into_painter(
                    frame_painter,
                    frame=frame,
                    scene=render_inputs.scene,
                    style=render_inputs.style,
                    hud=render_inputs.hud,
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
        cached_base_label_candidates = self._cached_base_label_candidates
        present_label_candidates: list[dict[str, object]] = []
        present_frame_key = SkyWindowRenderMixin._present_frame_cache_key(
            self,
            base_frame_key=base_frame_key,
            hud=render_inputs.hud,
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
                    present_label_candidates=present_label_candidates,
                    frame=frame,
                    render_inputs=render_inputs,
                ),
                setattr(
                    self,
                    "_cached_present_label_candidates",
                    list(present_label_candidates),
                ),
            ),
            cache_key_attr="_present_frame_cache_key",
            cache_image_attr="_present_frame_cache_image",
        )

    def _render_fast_frame_image(
        self,
        *,
        base_frame_key: tuple[object, ...],
        frame: FrameContext,
        render_inputs: RenderInputs,
    ) -> QImage:
        # Fast mode renders the heavy scene into a capped-size buffer and then
        # scales it up into the final window-sized frame.
        fast_frame_size = SkyWindowRenderMixin._fast_frame_render_size(self, frame)
        fast_viewport_rect = QRect(
            0,
            0,
            fast_frame_size.width(),
            fast_frame_size.height(),
        )
        fast_geometry = render_geometry.get_screen_geometry(
            fast_frame_size.width(),
            fast_frame_size.height(),
            frame.viewer.view_center[0],
        )
        fast_frame = FrameContext(
            viewer=frame.viewer,
            time_obj=frame.time_obj,
            geometry=fast_geometry,
            viewport_rect=fast_viewport_rect,
        )
        fast_base_frame_key = (
            base_frame_key,
            int(fast_frame_size.width()),
            int(fast_frame_size.height()),
            "fast-base",
        )
        overlay_time_bucket = None
        try:
            current_time_obj = frame.time_obj
            if current_time_obj is None:
                current_time_obj = self._current_time_obj()
            overlay_time_bucket = int(float(current_time_obj.unix) // 2.0)
        except Exception:
            overlay_time_bucket = None
        cache_stamp = self._render_cache_stamp
        fast_base_frame_image = SkyWindowRenderMixin._render_cached_image(
            self,
            image_size=fast_frame_size,
            frame_key=fast_base_frame_key,
            render_fn=lambda frame_painter: render_base_scene_into_painter(
                frame_painter,
                frame=fast_frame,
                scene=render_inputs.scene,
                style=render_inputs.style,
                hud=render_inputs.hud,
                compositor=self._compositor,
                draw_fast_overlays=False,
                label_candidates=[],
                draw_labels=False,
                draw_direction_labels=False,
            ),
            cache_key_attr="_fast_frame_base_cache_key",
            cache_image_attr="_fast_frame_base_cache_image",
        )
        return SkyWindowRenderMixin._render_cached_frame_image(
            self,
            frame_key=(
                "fast-present",
                base_frame_key,
                int(fast_frame_size.width()),
                int(fast_frame_size.height()),
                round(float(self.satellite_opacity), 3),
                round(float(self.aircraft_opacity), 3),
                overlay_time_bucket,
                cache_stamp(self.satellite_state.records_by_group),
                cache_stamp(self.aircraft_state.snapshots),
                cache_stamp(
                    SkyWindowRenderMixin._tropical_cyclone_snapshot_cache_value(self)
                ),
                self.tropical_cyclone_state.banner_text,
            ),
            render_fn=lambda frame_painter: (
                frame_painter.drawImage(frame.viewport_rect, fast_base_frame_image),
                render_fast_overlay_layers_into_painter(
                    frame_painter,
                    frame=frame,
                    scene=render_inputs.scene,
                    style=render_inputs.style,
                    draw_labels=True,
                ),
                render_guides.draw_direction_labels(
                    frame_painter,
                    frame.geometry,
                    frame.viewer,
                    render_inputs.style.text_font,
                    None,
                    theme=render_inputs.style.theme,
                ),
            ),
            cache_key_attr="_fast_frame_cache_key",
            cache_image_attr="_fast_frame_cache_image",
        )

    def _render_normal_frame_image(
        self,
        *,
        base_frame_key: tuple[object, ...],
        frame: FrameContext,
        render_inputs: RenderInputs,
        hover_targets: HoverTargets,
    ) -> QImage:
        return self._render_present_frame_image(
            base_frame_key=base_frame_key,
            frame=frame,
            render_inputs=render_inputs,
        )

    def _draw_present_frame_layers(
        self,
        *,
        frame_painter: QPainter,
        base_frame_image: QImage,
        base_label_candidates: list[dict[str, object]]
        | tuple[dict[str, object], ...]
        | None,
        present_label_candidates: list[dict[str, object]],
        frame: FrameContext,
        render_inputs: RenderInputs,
    ) -> None:
        frame_painter.drawImage(0, 0, base_frame_image)
        label_candidates: list[dict[str, object]] = list(base_label_candidates or [])
        render_fast_overlay_layers_into_painter(
            frame_painter,
            frame=frame,
            scene=render_inputs.scene,
            style=render_inputs.style,
            highlighted_satellite=None,
            highlighted_tropical_cyclone=None,
            label_candidates=label_candidates,
            draw_labels=False,
            draw_simplified_satellite_labels=_simplified_view_labels_visible(
                render_inputs.hud
            ),
        )
        present_label_candidates[:] = label_candidates

    def _draw_volatile_overlay_layers(
        self,
        painter: QPainter,
        *,
        frame: FrameContext,
        render_inputs: RenderInputs,
        hover_targets: HoverTargets,
    ) -> None:
        label_candidates = list(
            getattr(self, "_cached_present_label_candidates", ())
            or getattr(self, "_cached_base_label_candidates", ())
            or ()
        )
        render_hud_overlay_into_painter(
            painter,
            frame=frame,
            scene=render_inputs.scene,
            style=render_inputs.style,
            hud=render_inputs.hud,
            highlighted_object=hover_targets.object,
            highlighted_dso=hover_targets.dso,
            highlighted_satellite=hover_targets.satellite,
            highlighted_tropical_cyclone=hover_targets.tropical_cyclone,
            label_candidates=label_candidates,
            search_overlay_target=self.state.persistent_search_target,
        )

    def _compose_aircraft_debug_snapshot_image(
        self,
        present_frame: QImage,
        *,
        hover_targets: HoverTargets | None,
        frame: FrameContext,
        render_inputs: RenderInputs,
    ) -> QImage:
        debug_snapshot_frame = QImage(present_frame)
        if hover_targets is None:
            return debug_snapshot_frame
        debug_painter = QPainter(debug_snapshot_frame)
        debug_painter.setRenderHint(QPainter.Antialiasing)
        debug_painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        try:
            self._draw_volatile_overlay_layers(
                debug_painter,
                frame=frame,
                render_inputs=render_inputs,
                hover_targets=hover_targets,
            )
        finally:
            debug_painter.end()
        return debug_snapshot_frame

    def _frame_context_for_render(
        self, *, viewport_rect: QRect | None = None
    ) -> FrameContext:
        viewer = ViewerData(
            location=self.viewer_data.location,
            timezone_name=self.viewer_data.timezone_name,
            city_name=self.viewer_data.city_name,
            view_center=self.state.render_view_center,
            edge_fov_deg=self.viewer_data.edge_fov_deg,
            content_fov_deg=self.viewer_data.content_fov_deg,
            observer_height_m=self.viewer_data.observer_height_m,
            height_add_m=self.viewer_data.height_add_m,
            ground_elevation_m=self.viewer_data.ground_elevation_m,
            location_height_label=self.viewer_data.location_height_label,
            location_height_m=self.viewer_data.location_height_m,
        )
        if viewport_rect is None:
            viewport_rect = self.client_rect()
        geometry = render_geometry.get_screen_geometry(
            int(self.client_width()),
            int(self.client_height()),
            viewer.view_center[0],
        )
        time_obj = self._current_time_obj()
        return FrameContext(
            viewer=viewer,
            time_obj=time_obj,
            geometry=geometry,
            viewport_rect=viewport_rect,
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
            height_add_m=self.viewer_data.height_add_m,
            ground_elevation_m=self.viewer_data.ground_elevation_m,
            location_height_label=self.viewer_data.location_height_label,
            location_height_m=self.viewer_data.location_height_m,
        )

    def _render_inputs(
        self,
        *,
        celestial_data: CelestialData,
        frame: FrameContext,
    ) -> RenderInputs:
        return RenderInputs(
            scene=SkyWindowRenderMixin._render_scene_data(
                self,
                celestial_data=celestial_data,
                render_viewer=frame.viewer,
                time_obj=frame.time_obj,
            ),
            style=SkyWindowRenderMixin._render_style(self),
            hud=SkyWindowRenderMixin._render_hud_state(self),
        )

    def _render_scene_data(
        self,
        *,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
        time_obj: astropy.time.Time | None,
    ) -> RenderSceneData:
        state = self.state
        cloud_state = self._render_cloud_state()
        tropical_cyclone_state = self.tropical_cyclone_state
        tropical_cyclone_snapshots = tropical_cyclone_state.snapshots
        return RenderSceneData(
            viewer=render_viewer,
            celestial_data=celestial_data,
            sky_disc_image=state.sky_disc_image,
            cloud_missing_mask=cloud_state.missing_mask,
            cloud_altaz_grid=cloud_state.altaz_grid,
            terrain_horizon_profile=state.terrain_horizon_profile,
            terrain_horizon_profile_distances_m=state.terrain_horizon_profile_distances_m,
            terrain_secondary_ridges_altaz_layers=state.terrain_secondary_ridges_altaz_layers,
            terrain_secondary_ridges_distances_m_layers=state.terrain_secondary_ridges_distances_m_layers,
            urban_outlines=state.urban_outlines,
            water_overlay_dots=state.water_overlay_dots,
            tropical_cyclone_snapshots=tropical_cyclone_snapshots,
            satellite_element_epoch_utc=self.satellite_state.element_epoch_utc,
            satellite_records_by_group=self.satellite_state.records_by_group,
            aircraft_snapshots=self.aircraft_state.snapshots,
            time_obj=time_obj,
            night_light_glow_profile=state.night_light_glow_profile,
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
            bright_bodies_mode=str(self.bright_bodies_mode),
            star_base_radius=float(self.star_base_radius),
            star_visibility_boost=float(self.star_visibility_boost),
            asterism_visibility_boost=float(self.asterism_visibility_boost),
            earth_guide_visibility_boost=float(self.earth_guide_visibility_boost),
            vmag_limit=float(self.vmag_limit),
            sky_disc_altaz_rings=str(self.sky_disc_altaz_rings),
            sky_disc_altaz_rings_hover=str(self.sky_disc_altaz_rings_hover),
            cloud_disc_alpha=float(self.cloud_disc_alpha),
            satellite_opacity=float(self.satellite_opacity),
            terrain_horizon_opacity=float(self.terrain_horizon_opacity),
            earth_guide_opacity=float(self.earth_guide_opacity),
            night_light_opacity=float(self.night_light_opacity),
            ridge_glow_opacity=float(self.ridge_glow_opacity),
            urban_outline_opacity=float(self.urban_outline_opacity),
            show_urban_outline_layer=bool(self.show_urban_outline_layer),
            water_overlay_opacity=float(self.water_overlay_opacity),
            aircraft_opacity=float(self.aircraft_opacity),
            tropical_cyclone_opacity=float(self.tropical_cyclone_opacity),
            show_tropical_cyclone_overlay=bool(self.show_tropical_cyclone_overlay),
            star_render_expected_width=int(self._star_render_expected_width),
        )

    def _resolve_overlay_info_bottom_left(self, mouse_pos: QPoint | None) -> bool:
        overlay_info_bottom_left = bool(self.state.overlay_info_bottom_left)
        if mouse_pos is None:
            return overlay_info_bottom_left
        if self.observation_info_pinned:
            # Pinned CLI modes own the placement; mouse movement must not override it.
            return overlay_info_bottom_left
        window_height = max(1, int(self.client_height()))
        upper_threshold = float(window_height) / 3.0
        lower_threshold = 2.0 * float(window_height) / 3.0
        mouse_y = float(mouse_pos.y())
        if mouse_y <= upper_threshold:
            return True
        if mouse_y >= lower_threshold:
            return False
        return overlay_info_bottom_left

    def _render_hud_state(self) -> RenderHudState:
        status_message = self._status_line_message()
        mouse_pos = self.state.mouse_pos
        if self._startup_input_blocked():
            mouse_pos = None
        # Keep this state write centralized until overlay placement moves to input handlers.
        overlay_info_bottom_left = (
            SkyWindowRenderMixin._resolve_overlay_info_bottom_left(
                self,
                mouse_pos,
            )
        )
        self.state.overlay_info_bottom_left = overlay_info_bottom_left
        return RenderHudState(
            mouse_pos=mouse_pos,
            overlay_info_bottom_left=overlay_info_bottom_left,
            viewport_interaction_mode=bool(self.state.viewport_interaction_mode),
            viewport_interaction_stars=self.state.viewport_interaction_stars,
            simplified_view_enabled=bool(self._simplified_view_enabled()),
            simplified_view_labels_enabled=bool(self._simplified_view_labels_enabled()),
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

    def _draw_current_image_with_hud(
        self,
        painter: QPainter,
        *,
        frame: FrameContext,
        render_inputs: RenderInputs,
    ) -> None:
        label_candidates: list[dict[str, object]] = []
        render_base_scene_into_painter(
            painter,
            frame=frame,
            scene=render_inputs.scene,
            style=render_inputs.style,
            hud=render_inputs.hud,
            compositor=self._compositor,
            label_candidates=label_candidates,
            draw_labels=False,
        )
        render_fast_overlay_layers_into_painter(
            painter,
            frame=frame,
            scene=render_inputs.scene,
            style=render_inputs.style,
            highlighted_tropical_cyclone=None,
            label_candidates=label_candidates,
            draw_labels=True,
            draw_simplified_satellite_labels=_simplified_view_labels_visible(
                render_inputs.hud
            ),
        )
        highlighted_object = None
        highlighted_dso = None
        jump_highlight = self._active_jump_highlight_object(frame.geometry)
        if jump_highlight is not None:
            highlighted_object = jump_highlight
        render_hud_overlay_into_painter(
            painter,
            frame=frame,
            scene=render_inputs.scene,
            style=render_inputs.style,
            hud=render_inputs.hud,
            highlighted_object=highlighted_object,
            highlighted_dso=highlighted_dso,
            label_candidates=label_candidates,
            search_overlay_target=self.state.persistent_search_target,
        )

    def _draw_current_image_without_hud(
        self,
        painter: QPainter,
        *,
        frame: FrameContext,
        render_inputs: RenderInputs,
    ) -> None:
        render_base_scene_into_painter(
            painter,
            frame=frame,
            scene=render_inputs.scene,
            style=render_inputs.style,
            hud=render_inputs.hud,
            compositor=self._compositor,
        )
        render_fast_overlay_layers_into_painter(
            painter,
            frame=frame,
            scene=render_inputs.scene,
            style=render_inputs.style,
            highlighted_tropical_cyclone=None,
            draw_labels=True,
            draw_simplified_satellite_labels=_simplified_view_labels_visible(
                render_inputs.hud
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

            frame = self._frame_context_for_render(viewport_rect=self.client_rect())
            render_inputs = self._render_inputs(
                celestial_data=celestial_data,
                frame=frame,
            )
            if include_hud:
                self._draw_current_image_with_hud(
                    painter,
                    frame=frame,
                    render_inputs=render_inputs,
                )
            else:
                self._draw_current_image_without_hud(
                    painter,
                    frame=frame,
                    render_inputs=render_inputs,
                )
            return image
        finally:
            painter.end()

    def paintEvent(self, _event: QPaintEvent) -> None:
        if self._startup_splash_visible():
            return
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

        frame = self._frame_context_for_render(viewport_rect=self.client_rect())
        geometry = frame.geometry
        frame_key = self._render_frame_cache_key(
            frame=frame,
            celestial_data=celestial_data,
            include_fast_overlays=False,
        )
        self._update_star_render_stats(geometry)
        render_inputs = self._render_inputs(
            celestial_data=celestial_data,
            frame=frame,
        )
        hover_targets: HoverTargets | None = None
        if self.state.viewport_interaction_mode:
            present_frame = self._render_fast_frame_image(
                base_frame_key=frame_key,
                frame=frame,
                render_inputs=render_inputs,
            )
            hover_targets = HoverTargets()
        else:
            mouse_pos = self.state.mouse_pos
            if self._startup_input_blocked():
                mouse_pos = None

            hover_targets = _resolve_hover_targets(
                celestial_data=celestial_data,
                render_viewer=frame.viewer,
                mouse_pos=mouse_pos,
                geometry=geometry,
                satellite_records_by_group=self.satellite_state.records_by_group,
                tropical_cyclone_snapshots=self.tropical_cyclone_state.snapshots,
                time_obj=render_inputs.scene.time_obj,
                show_dso=bool(self.show_dso),
            )
            jump_highlight = self._active_jump_highlight_object(geometry)
            if jump_highlight is not None:
                hover_targets = HoverTargets(
                    object=jump_highlight,
                    dso=hover_targets.dso,
                    satellite=hover_targets.satellite,
                    tropical_cyclone=hover_targets.tropical_cyclone,
                )
            present_frame = self._render_normal_frame_image(
                base_frame_key=frame_key,
                frame=frame,
                render_inputs=render_inputs,
                hover_targets=hover_targets,
            )
        painter.drawImage(0, 0, present_frame)
        if hover_targets is not None:
            self._draw_volatile_overlay_layers(
                painter,
                frame=frame,
                render_inputs=render_inputs,
                hover_targets=hover_targets,
            )
        if self._pending_aircraft_debug_snapshot_path is not None:
            debug_snapshot_frame = SkyWindowRenderMixin._compose_aircraft_debug_snapshot_image(
                self,
                present_frame,
                hover_targets=hover_targets,
                frame=frame,
                render_inputs=render_inputs,
            )
            self._flush_aircraft_debug_snapshot_save(debug_snapshot_frame)
