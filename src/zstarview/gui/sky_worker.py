"""
Background worker for celestial/sky-disc calculations.

This module extracts heavy sky computations from the main window so UI code can
focus on orchestration and rendering.
"""
from __future__ import annotations

import logging
import math
import threading
import time
from collections.abc import Callable
from concurrent.futures import Future
from datetime import datetime, timedelta, timezone

import astropy
import astropy.time
import numpy as np
import polars as pl
from PySide6.QtCore import QObject, Signal
from PySide6.QtGui import QImage

from ..astro import (
    DeepSkyCatalogArrays,
    StarCatalogArrays,
    calculate_celestial_equator_points,
    calculate_ecliptic_points,
    calculate_horizon_points,
    calculate_planets,
    calculate_visible_deep_sky_objects,
    calculate_visible_stars,
    eclipse_factor_from_info,
)
from ..night_lights import compute_night_light_glow_profile, is_night_light_enabled
from ..paths import ThemeStyle
from ..render import sky_disc
from ..render.aerosol_profile import bundled_aod550_or_default
from ..types import CelestialData, ScreenGeometry, StarCatalogMeta, ViewerData
from .application_services import ApplicationServices
from .worker_pool import wait_for_gui_futures

logger = logging.getLogger(__name__)


def _sky_disc_render_surface(
    geometry: ScreenGeometry,
    image_size: tuple[int, int] | None,
    render_scale: float = sky_disc.SKY_DISC_RENDER_SCALE,
) -> tuple[ScreenGeometry, tuple[int, int] | None]:
    """Return the reduced sky-disc surface geometry and image size."""
    if image_size is None:
        return geometry, None
    render_scale = max(0.01, min(1.0, float(render_scale)))
    disc_width_px = max(2, int(geometry.radius) * 2)
    if disc_width_px > sky_disc.SKY_DISC_REFERENCE_WIDTH_PX:
        render_scale *= math.sqrt(
            sky_disc.SKY_DISC_REFERENCE_WIDTH_PX / float(disc_width_px)
        )
        render_scale = max(0.01, min(1.0, render_scale))
    render_size = (
        max(2, int(math.ceil(max(2, int(image_size[0])) * render_scale))),
        max(2, int(math.ceil(max(2, int(image_size[1])) * render_scale))),
    )
    render_geometry = ScreenGeometry(
        center=(
            int(round(geometry.center[0] * render_scale)),
            int(round(geometry.center[1] * render_scale)),
        ),
        radius=max(1, int(math.ceil(geometry.radius * render_scale))),
    )
    return render_geometry, render_size


def compute_sky_snapshot(
    *,
    ephemeris: object,
    viewer_data: ViewerData,
    geometry: ScreenGeometry,
    star_catalog: pl.DataFrame | StarCatalogArrays,
    dso_catalog: DeepSkyCatalogArrays | None,
    star_vmag_limit: float | None,
    star_subset_indices: np.ndarray | None,
    star_data_policy: str = "scenic_view_scoped",
    delta_t: timedelta,
    sky_update_interval: float = 60.0,
    sky_disc_alpha: float,
    theme: ThemeStyle,
    star_catalog_meta: StarCatalogMeta | None = None,
    image_size: tuple[int, int] | None = None,
    sky_disc_render_scale: float = 1.0,
    terrain_horizon_profile_altaz: list[tuple[float, float]] | None = None,
    terrain_horizon_profile_distances_m: list[float] | None = None,
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
    terrain_secondary_ridges_distances_m_layers: list[list[float]] | None = None,
    terrain_sample_distances_m: np.ndarray | None = None,
    terrain_sample_terrain_elevation_m: np.ndarray | None = None,
    night_light_glow_profile: object | None = None,
    night_light_opacity: float = 0.0,
    render_generation: int = 0,
) -> dict[str, object]:
    """Compute celestial data and sky-disc image synchronously."""
    now = datetime.now(timezone.utc) + delta_t
    time_obj = astropy.time.Time(now)
    lat, lon = viewer_data.location
    view_center = viewer_data.view_center
    observer_height_m = float(viewer_data.observer_height_m)
    observer_elevation_m = max(
        0.0,
        float(viewer_data.ground_elevation_m) + observer_height_m,
    )
    edge_fov_deg = float(viewer_data.edge_fov_deg)
    content_fov_deg = float(viewer_data.content_fov_deg)
    star_time_obj = time_obj + astropy.time.TimeDelta(
        max(0.0, float(sky_update_interval)) / 2.0,
        format="sec",
    )

    stars, loc = calculate_visible_stars(
        star_catalog,
        lat,
        lon,
        observer_height_m,
        star_time_obj,
        view_center,
        content_fov_deg=content_fov_deg,
        max_vmag=star_vmag_limit,
        subset_indices=star_subset_indices,
        star_data_policy=star_data_policy,
    )
    if dso_catalog is None:
        empty_obj = np.array([], dtype=object)
        empty_float = np.array([], dtype=float)
        deep_sky_objects = {
            "id": empty_obj,
            "name": empty_obj,
            "type": empty_obj,
            "alt": empty_float,
            "az": empty_float,
            "vmag": empty_float,
            "major_arcmin": empty_float,
            "minor_arcmin": empty_float,
            "pa_deg": empty_float,
        }
    else:
        deep_sky_objects = calculate_visible_deep_sky_objects(
            dso_catalog,
            lat,
            lon,
            observer_height_m,
            time_obj,
            view_center,
            content_fov_deg=content_fov_deg,
        )
    planets = calculate_planets(
        lat,
        lon,
        observer_height_m,
        time_obj,
        view_center,
        ephemeris,
        content_fov_deg=content_fov_deg,
    )
    celestial_equator_points = calculate_celestial_equator_points(loc, time_obj)
    ecliptic_points = calculate_ecliptic_points(loc, time_obj)
    horizon_points = calculate_horizon_points()
    celestial_data = CelestialData(
        time=time_obj,
        planets=planets,
        stars=stars,
        deep_sky_objects=deep_sky_objects,
        celestial_equator_points=celestial_equator_points,
        ecliptic_points=ecliptic_points,
        horizon_points=horizon_points,
        star_catalog_meta=star_catalog_meta,
        star_time=star_time_obj,
    )

    sun_altaz = None
    solar_eclipse_info = None
    for body in planets:
        if body.name == "sun":
            sun_altaz = (body.alt, body.az)
            solar_eclipse_info = body.solar_eclipse_info
            break

    sky_disc_img: QImage | None = None
    cached_night_light_glow_profile = night_light_glow_profile
    if sun_altaz is not None:
        terrain_ready_for_night_light = bool(
            terrain_horizon_profile_altaz
            and terrain_horizon_profile_distances_m is not None
        )
        sky_disc_geometry, render_image_size = _sky_disc_render_surface(
            geometry,
            image_size,
            sky_disc_render_scale,
        )
        ef = eclipse_factor_from_info(solar_eclipse_info)
        disc_opacity = float(theme.sky_disc.opacity)
        aerosol_optical_depth = bundled_aod550_or_default(
            float(lat),
            float(lon),
            int(time_obj.datetime.month),
        )
        if sky_disc_alpha > 0.0:
            sky_disc_img = sky_disc.draw_sky_color_disc(
                sky_disc_geometry,
                view_center,
                edge_fov_deg=edge_fov_deg,
                content_fov_deg=content_fov_deg + sky_disc.SKY_DISC_OVERSCAN_DEG,
                sun_altaz=sun_altaz,
                alpha=sky_disc_alpha,
                disc_opacity=disc_opacity,
                eclipse_factor=ef,
                observer_height_m=observer_elevation_m,
                time_obj=time_obj,
                timezone_name=viewer_data.timezone_name,
                image_size=render_image_size,
                aerosol_optical_depth=aerosol_optical_depth,
            )
        else:
            sky_disc_img = sky_disc.draw_uniform_sky_color_disc(
                sky_disc_geometry,
                view_center,
                edge_fov_deg=edge_fov_deg,
                content_fov_deg=content_fov_deg + sky_disc.SKY_DISC_OVERSCAN_DEG,
                disc_opacity=disc_opacity,
                image_size=render_image_size,
            )
        if (
            cached_night_light_glow_profile is None
            and is_night_light_enabled(float(sun_altaz[0]))
            and terrain_ready_for_night_light
        ):
            try:
                logger.info("Calculating night light alpha grid...")
                cached_night_light_glow_profile = compute_night_light_glow_profile(
                    observer_lat_deg=float(lat),
                    observer_lon_deg=float(lon),
                    observer_elevation_m=observer_elevation_m,
                    sun_alt_deg=float(sun_altaz[0]),
                    terrain_profile_altaz=terrain_horizon_profile_altaz,
                    terrain_profile_distances_m=terrain_horizon_profile_distances_m,
                    terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
                    terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
                    terrain_sample_distances_m=terrain_sample_distances_m,
                    terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
                    include_night_light_tiles=float(night_light_opacity) > 0.0,
                )
                if cached_night_light_glow_profile is not None:
                    logger.info("Night light alpha grid computed")
            except Exception as exc:
                logger.warning("Night light overlay unavailable: %s", exc)

    payload: dict[str, object] = {
        "celestial": celestial_data,
        "sky_disc": sky_disc_img,
        "night_light_glow_profile": cached_night_light_glow_profile,
    }
    payload["view_center"] = (float(view_center[0]), float(view_center[1]))
    payload["geometry"] = geometry
    payload["render_generation"] = int(render_generation)
    return payload


def compute_sky_disc_snapshot(
    *,
    ephemeris: object,
    viewer_data: ViewerData,
    geometry: ScreenGeometry,
    delta_t: timedelta,
    sky_disc_alpha: float,
    theme: ThemeStyle,
    image_size: tuple[int, int] | None = None,
    sky_disc_render_scale: float = sky_disc.SKY_DISC_RENDER_SCALE,
    render_generation: int = 0,
) -> dict[str, object]:
    """Compute only the low-resolution sky-colour disc."""
    now = datetime.now(timezone.utc) + delta_t
    time_obj = astropy.time.Time(now)
    lat, lon = viewer_data.location
    observer_height_m = float(viewer_data.observer_height_m)
    observer_elevation_m = max(
        0.0,
        float(viewer_data.ground_elevation_m) + observer_height_m,
    )
    planets = calculate_planets(
        lat,
        lon,
        observer_height_m,
        time_obj,
        viewer_data.view_center,
        ephemeris,
        content_fov_deg=float(viewer_data.content_fov_deg),
    )
    sun_altaz = None
    solar_eclipse_info = None
    for body in planets:
        if body.name == "sun":
            sun_altaz = (body.alt, body.az)
            solar_eclipse_info = body.solar_eclipse_info
            break

    sky_disc_img: QImage | None = None
    if sun_altaz is not None:
        sky_disc_geometry, render_image_size = _sky_disc_render_surface(
            geometry,
            image_size,
            sky_disc_render_scale,
        )
        aerosol_optical_depth = bundled_aod550_or_default(
            float(lat),
            float(lon),
            int(time_obj.datetime.month),
        )
        eclipse_factor = eclipse_factor_from_info(solar_eclipse_info)
        if sky_disc_alpha > 0.0:
            sky_disc_img = sky_disc.draw_sky_color_disc(
                sky_disc_geometry,
                viewer_data.view_center,
                edge_fov_deg=float(viewer_data.edge_fov_deg),
                content_fov_deg=float(viewer_data.content_fov_deg)
                + sky_disc.SKY_DISC_OVERSCAN_DEG,
                sun_altaz=sun_altaz,
                alpha=float(sky_disc_alpha),
                disc_opacity=float(theme.sky_disc.opacity),
                eclipse_factor=eclipse_factor,
                observer_height_m=observer_elevation_m,
                time_obj=time_obj,
                timezone_name=viewer_data.timezone_name,
                image_size=render_image_size,
                aerosol_optical_depth=aerosol_optical_depth,
            )
        else:
            sky_disc_img = sky_disc.draw_uniform_sky_color_disc(
                sky_disc_geometry,
                viewer_data.view_center,
                edge_fov_deg=float(viewer_data.edge_fov_deg),
                content_fov_deg=float(viewer_data.content_fov_deg)
                + sky_disc.SKY_DISC_OVERSCAN_DEG,
                disc_opacity=float(theme.sky_disc.opacity),
                image_size=render_image_size,
            )

    return {
        "sky_disc": sky_disc_img,
        "sun_alt_deg": None if sun_altaz is None else float(sun_altaz[0]),
        "view_center": tuple(viewer_data.view_center),
        "geometry": geometry,
        "render_generation": int(render_generation),
    }


class SkyDataWorker(QObject):
    """Compute sky data in a Python background thread and emit results."""

    data_ready = Signal(object)
    sky_disc_ready = Signal(object)
    planet_data_ready = Signal(object)

    def __init__(
        self,
        services: ApplicationServices | None = None,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._owns_services = services is None
        self._services = services or ApplicationServices()
        self._lock = threading.Lock()
        self._running = False
        self._sky_disc_running = False
        self._planet_running = False
        self._stopping = False
        self._active_workers: set[Future[None]] = set()

    def shutdown(self, *, wait_timeout_s: float | None = None) -> None:
        """Stop accepting/emitting updates during application shutdown."""
        with self._lock:
            self._stopping = True
        self._wait_for_workers(wait_timeout_s)
        if self._owns_services:
            self._services.shutdown(wait=True)

    def has_in_flight_update(self) -> bool:
        with self._lock:
            return bool(
                self._running
                or self._sky_disc_running
                or self._planet_running
                or self._active_workers
            )

    def update_sky_disc(
        self,
        *,
        ephemeris: object,
        viewer_data: ViewerData,
        geometry: ScreenGeometry,
        delta_t: timedelta,
        sky_disc_alpha: float,
        theme: ThemeStyle,
        image_size: tuple[int, int] | None = None,
        sky_disc_render_scale: float = sky_disc.SKY_DISC_RENDER_SCALE,
        render_generation: int = 0,
    ) -> bool:
        """Start a background update for only the sky-colour disc."""
        with self._lock:
            if self._stopping or self._running or self._sky_disc_running:
                return False
            self._sky_disc_running = True
        self._spawn_worker(
            target=self._run_sky_disc_update,
            kwargs={
                "ephemeris": ephemeris,
                "viewer_data": viewer_data,
                "geometry": geometry,
                "delta_t": delta_t,
                "sky_disc_alpha": sky_disc_alpha,
                "theme": theme,
                "image_size": image_size,
                "sky_disc_render_scale": sky_disc_render_scale,
                "render_generation": render_generation,
            },
        )
        return True

    def update_planets(
        self,
        *,
        ephemeris: object,
        viewer_data: ViewerData,
        time_obj: astropy.time.Time,
        render_generation: int,
    ) -> bool:
        """Calculate only solar-system positions for a display refresh."""
        with self._lock:
            if self._stopping or self._planet_running:
                return False
            self._planet_running = True
        self._spawn_worker(
            target=self._run_planet_update,
            kwargs={
                "ephemeris": ephemeris,
                "viewer_data": viewer_data,
                "time_obj": time_obj,
                "render_generation": render_generation,
            },
        )
        return True

    def update(
        self,
        *,
        ephemeris: object,
        viewer_data: ViewerData,
        geometry: ScreenGeometry,
        star_catalog: pl.DataFrame | StarCatalogArrays,
        dso_catalog: DeepSkyCatalogArrays | None = None,
        star_vmag_limit: float | None = None,
        star_subset_indices: np.ndarray | None = None,
        star_data_policy: str = "scenic_view_scoped",
        delta_t: timedelta,
        sky_update_interval: float = 60.0,
        sky_disc_alpha: float,
        theme: ThemeStyle,
        star_catalog_meta: StarCatalogMeta | None = None,
        image_size: tuple[int, int] | None = None,
        sky_disc_render_scale: float = sky_disc.SKY_DISC_RENDER_SCALE,
        terrain_horizon_profile_altaz: list[tuple[float, float]] | None = None,
        terrain_horizon_profile_distances_m: list[float] | None = None,
        terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
        terrain_secondary_ridges_distances_m_layers: list[list[float]] | None = None,
        terrain_sample_distances_m: np.ndarray | None = None,
        terrain_sample_terrain_elevation_m: np.ndarray | None = None,
        night_light_glow_profile: object | None = None,
        night_light_opacity: float = 0.0,
        render_generation: int = 0,
    ) -> bool:
        """Start background computation if idle; return False when already running."""
        with self._lock:
            if self._stopping or self._running:
                return False
            self._running = True

        self._spawn_worker(
            target=self._run_update,
            kwargs={
                "ephemeris": ephemeris,
                "viewer_data": viewer_data,
                "geometry": geometry,
                "star_catalog": star_catalog,
                "dso_catalog": dso_catalog,
                "star_vmag_limit": star_vmag_limit,
                "star_subset_indices": star_subset_indices,
                "star_data_policy": star_data_policy,
                "delta_t": delta_t,
                "sky_update_interval": sky_update_interval,
                "sky_disc_alpha": sky_disc_alpha,
                "theme": theme,
                "star_catalog_meta": star_catalog_meta,
                "image_size": image_size,
                "sky_disc_render_scale": sky_disc_render_scale,
                "terrain_horizon_profile_altaz": terrain_horizon_profile_altaz,
                "terrain_horizon_profile_distances_m": terrain_horizon_profile_distances_m,
                "terrain_secondary_ridges_altaz_layers": terrain_secondary_ridges_altaz_layers,
                "terrain_secondary_ridges_distances_m_layers": terrain_secondary_ridges_distances_m_layers,
                "terrain_sample_distances_m": terrain_sample_distances_m,
                "terrain_sample_terrain_elevation_m": terrain_sample_terrain_elevation_m,
                "night_light_glow_profile": night_light_glow_profile,
                "night_light_opacity": night_light_opacity,
                "render_generation": render_generation,
            },
        )
        return True

    def _spawn_worker(
        self,
        *,
        target: Callable[..., None],
        kwargs: dict[str, object],
    ) -> None:
        future = self._services.submit(target, **kwargs)
        with self._lock:
            if self._stopping:
                return
            self._active_workers.add(future)
            if future.done():
                self._active_workers.discard(future)
                return
        future.add_done_callback(self._unregister_worker)

    def _unregister_worker(self, worker: Future[None]) -> None:
        with self._lock:
            self._active_workers.discard(worker)

    def _wait_for_workers(self, wait_timeout_s: float | None) -> None:
        deadline = None if wait_timeout_s is None else time.monotonic() + max(0.0, float(wait_timeout_s))
        while True:
            with self._lock:
                workers = tuple(self._active_workers)
            if not workers:
                return
            if deadline is None:
                wait_for_gui_futures(workers, None)
                continue
            remaining = deadline - time.monotonic()
            if remaining <= 0.0:
                logger.warning(
                    "Timed out waiting for %d sky worker task(s) to finish during shutdown",
                    len(workers),
                )
                return
            wait_for_gui_futures(workers, remaining)

    def _run_planet_update(
        self,
        *,
        ephemeris: object,
        viewer_data: ViewerData,
        time_obj: astropy.time.Time,
        render_generation: int,
    ) -> None:
        try:
            with self._services.native_work_lock:
                planets = calculate_planets(
                    viewer_data.location[0],
                    viewer_data.location[1],
                    float(viewer_data.observer_height_m),
                    time_obj,
                    viewer_data.view_center,
                    ephemeris,
                    content_fov_deg=float(viewer_data.content_fov_deg),
                )
            with self._lock:
                if self._stopping:
                    return
            self.planet_data_ready.emit(
                {
                    "planets": planets,
                    "time_unix": float(time_obj.unix),
                    "render_generation": render_generation,
                }
            )
        except Exception:
            logger.exception("Error in planet position update")
        finally:
            with self._lock:
                self._planet_running = False

    def _run_sky_disc_update(
        self,
        *,
        ephemeris: object,
        viewer_data: ViewerData,
        geometry: ScreenGeometry,
        delta_t: timedelta,
        sky_disc_alpha: float,
        theme: ThemeStyle,
        image_size: tuple[int, int] | None,
        sky_disc_render_scale: float,
        render_generation: int,
    ) -> None:
        try:
            with self._services.native_work_lock:
                payload = compute_sky_disc_snapshot(
                    ephemeris=ephemeris,
                    viewer_data=viewer_data,
                    geometry=geometry,
                    delta_t=delta_t,
                    sky_disc_alpha=sky_disc_alpha,
                    theme=theme,
                    image_size=image_size,
                    sky_disc_render_scale=sky_disc_render_scale,
                    render_generation=render_generation,
                )
            with self._lock:
                if self._stopping:
                    return
            self.sky_disc_ready.emit(payload)
        except Exception:
            logger.exception("Error in background sky-disc update thread")
        finally:
            with self._lock:
                self._sky_disc_running = False

    def _run_update(
        self,
        *,
        ephemeris: object,
        viewer_data: ViewerData,
        geometry: ScreenGeometry,
        star_catalog: pl.DataFrame | StarCatalogArrays,
        dso_catalog: DeepSkyCatalogArrays | None,
        star_vmag_limit: float | None,
        star_subset_indices: np.ndarray | None,
        star_data_policy: str,
        delta_t: timedelta,
        sky_update_interval: float,
        sky_disc_alpha: float,
        theme: ThemeStyle,
        star_catalog_meta: StarCatalogMeta | None,
        image_size: tuple[int, int] | None,
        sky_disc_render_scale: float,
        terrain_horizon_profile_altaz: list[tuple[float, float]] | None,
        terrain_horizon_profile_distances_m: list[float] | None,
        terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None,
        terrain_secondary_ridges_distances_m_layers: list[list[float]] | None,
        terrain_sample_distances_m: np.ndarray | None,
        terrain_sample_terrain_elevation_m: np.ndarray | None,
        night_light_glow_profile: object | None,
        night_light_opacity: float,
        render_generation: int,
    ) -> None:
        try:
            with self._services.native_work_lock:
                payload = compute_sky_snapshot(
                    ephemeris=ephemeris,
                    viewer_data=viewer_data,
                    geometry=geometry,
                    star_catalog=star_catalog,
                    dso_catalog=dso_catalog,
                    star_vmag_limit=star_vmag_limit,
                    star_subset_indices=star_subset_indices,
                    star_data_policy=star_data_policy,
                    delta_t=delta_t,
                    sky_update_interval=sky_update_interval,
                    sky_disc_alpha=sky_disc_alpha,
                    theme=theme,
                    star_catalog_meta=star_catalog_meta,
                    image_size=image_size,
                    sky_disc_render_scale=sky_disc_render_scale,
                    terrain_horizon_profile_altaz=terrain_horizon_profile_altaz,
                    terrain_horizon_profile_distances_m=terrain_horizon_profile_distances_m,
                    terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
                    terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
                    terrain_sample_distances_m=terrain_sample_distances_m,
                    terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
                    night_light_glow_profile=night_light_glow_profile,
                    night_light_opacity=float(night_light_opacity),
                    render_generation=render_generation,
                )
            with self._lock:
                if self._stopping:
                    return
            try:
                self.data_ready.emit(payload)
            except RuntimeError:
                # QObject can be deleted while background thread is still unwinding.
                logger.debug("Skip sky data emit during shutdown.")
        except Exception:
            logger.exception("Error in background sky update thread")
        finally:
            with self._lock:
                self._running = False
