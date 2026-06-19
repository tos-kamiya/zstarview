# -*- coding: utf-8 -*-
"""
Background worker for celestial/sky-disc calculations.

This module extracts heavy sky computations from the main window so UI code can
focus on orchestration and rendering.
"""
from __future__ import annotations

import logging
import threading
import time
from concurrent.futures import Future
from datetime import datetime, timedelta, timezone
from typing import Callable, Dict

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
from ..night_lights import compute_night_light_glow_profile
from ..paths import ThemeStyle
from ..render import sky_disc
from ..types import CelestialData, ScreenGeometry, StarCatalogMeta, ViewerData
from .native_work_lock import HEAVY_NATIVE_WORK_LOCK
from .worker_pool import submit_gui_work, wait_for_gui_futures

logger = logging.getLogger(__name__)


def compute_sky_snapshot(
    *,
    ephemeris: object,
    viewer_data: ViewerData,
    geometry: ScreenGeometry,
    star_catalog: pl.DataFrame | StarCatalogArrays,
    dso_catalog: DeepSkyCatalogArrays | None,
    star_vmag_limit: float | None,
    star_subset_indices: np.ndarray | None,
    delta_t: timedelta,
    sky_disc_alpha: float,
    theme: ThemeStyle,
    star_catalog_meta: StarCatalogMeta | None = None,
    image_size: tuple[int, int] | None = None,
    terrain_horizon_profile_altaz: list[tuple[float, float]] | None = None,
    terrain_horizon_profile_distances_m: list[float] | None = None,
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
    terrain_secondary_ridges_distances_m_layers: list[list[float]] | None = None,
    terrain_sample_distances_m: np.ndarray | None = None,
    terrain_sample_terrain_elevation_m: np.ndarray | None = None,
    night_light_glow_profile: object | None = None,
    render_generation: int = 0,
) -> Dict[str, object]:
    """Compute celestial data and sky-disc image synchronously."""
    now = datetime.now(timezone.utc) + delta_t
    time_obj = astropy.time.Time(now)
    lat, lon = viewer_data.location
    view_center = viewer_data.view_center
    observer_height_m = float(viewer_data.observer_height_m)
    edge_fov_deg = float(viewer_data.edge_fov_deg)
    content_fov_deg = float(viewer_data.content_fov_deg)

    stars, loc = calculate_visible_stars(
        star_catalog,
        lat,
        lon,
        observer_height_m,
        time_obj,
        view_center,
        content_fov_deg=content_fov_deg,
        max_vmag=star_vmag_limit,
        subset_indices=star_subset_indices,
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
        render_image_size = (
            max(2, int(image_size[0])),
            max(2, int(image_size[1])),
        ) if image_size is not None else None
        ef = eclipse_factor_from_info(solar_eclipse_info)
        disc_opacity = float(theme.sky_disc.opacity)
        if sky_disc_alpha > 0.0:
            sky_disc_img = sky_disc.draw_sky_color_disc(
                geometry,
                view_center,
                sun_altaz,
                alpha=sky_disc_alpha,
                disc_opacity=disc_opacity,
                eclipse_factor=ef,
                edge_fov_deg=edge_fov_deg,
                content_fov_deg=content_fov_deg,
                image_size=render_image_size,
            )
        else:
            sky_disc_img = sky_disc.draw_uniform_sky_color_disc(
                geometry,
                view_center,
                edge_fov_deg=edge_fov_deg,
                content_fov_deg=content_fov_deg,
                image_size=render_image_size,
                disc_opacity=disc_opacity,
            )
        if cached_night_light_glow_profile is None and float(sun_altaz[0]) < 0.0:
            try:
                cached_night_light_glow_profile = compute_night_light_glow_profile(
                    observer_lat_deg=float(lat),
                    observer_lon_deg=float(lon),
                    sun_alt_deg=float(sun_altaz[0]),
                    terrain_profile_altaz=terrain_horizon_profile_altaz,
                    terrain_profile_distances_m=terrain_horizon_profile_distances_m,
                    terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
                    terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
                    terrain_sample_distances_m=terrain_sample_distances_m,
                    terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
                )
                if cached_night_light_glow_profile is not None:
                    logger.info("[INFO] Night light alpha grid computed")
            except Exception as exc:
                logger.warning("Night light overlay unavailable: %s", exc)

    payload: Dict[str, object] = {
        "celestial": celestial_data,
        "sky_disc": sky_disc_img,
        "night_light_glow_profile": cached_night_light_glow_profile,
    }
    payload["view_center"] = (float(view_center[0]), float(view_center[1]))
    payload["geometry"] = geometry
    payload["render_generation"] = int(render_generation)
    return payload


class SkyDataWorker(QObject):
    """Compute sky data in a Python background thread and emit results."""

    data_ready = Signal(object)

    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._lock = threading.Lock()
        self._running = False
        self._stopping = False
        self._active_workers: set[Future[None]] = set()

    def shutdown(self, *, wait_timeout_s: float | None = None) -> None:
        """Stop accepting/emitting updates during application shutdown."""
        with self._lock:
            self._stopping = True
        self._wait_for_workers(wait_timeout_s)

    def has_in_flight_update(self) -> bool:
        with self._lock:
            return bool(self._running or self._active_workers)

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
        delta_t: timedelta,
        sky_disc_alpha: float,
        theme: ThemeStyle,
        star_catalog_meta: StarCatalogMeta | None = None,
        image_size: tuple[int, int] | None = None,
        terrain_horizon_profile_altaz: list[tuple[float, float]] | None = None,
        terrain_horizon_profile_distances_m: list[float] | None = None,
        terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
        terrain_secondary_ridges_distances_m_layers: list[list[float]] | None = None,
        terrain_sample_distances_m: np.ndarray | None = None,
        terrain_sample_terrain_elevation_m: np.ndarray | None = None,
        night_light_glow_profile: object | None = None,
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
                "delta_t": delta_t,
                "sky_disc_alpha": sky_disc_alpha,
                "theme": theme,
                "star_catalog_meta": star_catalog_meta,
                "image_size": image_size,
                "terrain_horizon_profile_altaz": terrain_horizon_profile_altaz,
                "terrain_horizon_profile_distances_m": terrain_horizon_profile_distances_m,
                "terrain_secondary_ridges_altaz_layers": terrain_secondary_ridges_altaz_layers,
                "terrain_secondary_ridges_distances_m_layers": terrain_secondary_ridges_distances_m_layers,
                "terrain_sample_distances_m": terrain_sample_distances_m,
                "terrain_sample_terrain_elevation_m": terrain_sample_terrain_elevation_m,
                "night_light_glow_profile": night_light_glow_profile,
                "render_generation": render_generation,
            },
            label="sky",
        )
        return True

    def _spawn_worker(
        self,
        *,
        target: Callable[..., None],
        kwargs: dict[str, object],
        label: str,
    ) -> None:
        def runner() -> None:
            target(**kwargs)

        future = submit_gui_work(runner)
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
        delta_t: timedelta,
        sky_disc_alpha: float,
        theme: ThemeStyle,
        star_catalog_meta: StarCatalogMeta | None,
        image_size: tuple[int, int] | None,
        terrain_horizon_profile_altaz: list[tuple[float, float]] | None,
        terrain_horizon_profile_distances_m: list[float] | None,
        terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None,
        terrain_secondary_ridges_distances_m_layers: list[list[float]] | None,
        terrain_sample_distances_m: np.ndarray | None,
        terrain_sample_terrain_elevation_m: np.ndarray | None,
        night_light_glow_profile: object | None,
        render_generation: int,
    ) -> None:
        try:
            with HEAVY_NATIVE_WORK_LOCK:
                payload = compute_sky_snapshot(
                    ephemeris=ephemeris,
                    viewer_data=viewer_data,
                    geometry=geometry,
                    star_catalog=star_catalog,
                    dso_catalog=dso_catalog,
                    star_vmag_limit=star_vmag_limit,
                    star_subset_indices=star_subset_indices,
                    delta_t=delta_t,
                    sky_disc_alpha=sky_disc_alpha,
                    theme=theme,
                    star_catalog_meta=star_catalog_meta,
                    image_size=image_size,
                    terrain_horizon_profile_altaz=terrain_horizon_profile_altaz,
                    terrain_horizon_profile_distances_m=terrain_horizon_profile_distances_m,
                    terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
                    terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
                    terrain_sample_distances_m=terrain_sample_distances_m,
                    terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
                    night_light_glow_profile=night_light_glow_profile,
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
        except Exception as e:
            logger.error("Error in background sky update thread: %s", e, exc_info=True)
        finally:
            with self._lock:
                self._running = False
