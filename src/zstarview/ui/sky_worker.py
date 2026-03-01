# -*- coding: utf-8 -*-
"""
Background worker for celestial/sky-disc calculations.

This module extracts heavy sky computations from the main window so UI code can
focus on orchestration and rendering.
"""
from __future__ import annotations

import logging
import threading
from datetime import datetime, timedelta, timezone
from typing import Dict, Tuple

import astropy
import polars as pl
from PySide6.QtCore import QObject, Signal
from PySide6.QtGui import QImage

from ..astro import (
    StarCatalogArrays,
    calculate_celestial_equator_points,
    calculate_ecliptic_points,
    calculate_horizon_points,
    calculate_planets,
    calculate_visible_stars,
    eclipse_factor_from_info,
)
from ..render import draw as render_draw
from ..render import draw_sky_disc
from ..types import CelestialData

logger = logging.getLogger(__name__)


class SkyDataWorker(QObject):
    """Compute sky data in a Python background thread and emit results."""

    data_ready = Signal(object)

    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._lock = threading.Lock()
        self._running = False
        self._stopping = False

    def shutdown(self) -> None:
        """Stop accepting/emitting updates during application shutdown."""
        with self._lock:
            self._stopping = True

    def update(
        self,
        *,
        lat: float,
        lon: float,
        view_center: Tuple[float, float],
        star_catalog: pl.DataFrame | StarCatalogArrays,
        star_vmag_limit: float | None = None,
        delta_t: timedelta,
        sky_disc_alpha: float,
        sky_disc_base_size: int,
        debug_eclipses: bool = False,
    ) -> bool:
        """Start background computation if idle; return False when already running."""
        with self._lock:
            if self._stopping or self._running:
                return False
            self._running = True

        t = threading.Thread(
            target=self._run_update,
            kwargs={
                "lat": lat,
                "lon": lon,
                "view_center": view_center,
                "star_catalog": star_catalog,
                "star_vmag_limit": star_vmag_limit,
                "delta_t": delta_t,
                "sky_disc_alpha": sky_disc_alpha,
                "sky_disc_base_size": sky_disc_base_size,
                "debug_eclipses": debug_eclipses,
            },
            daemon=True,
        )
        t.start()
        return True

    def _run_update(
        self,
        *,
        lat: float,
        lon: float,
        view_center: Tuple[float, float],
        star_catalog: pl.DataFrame | StarCatalogArrays,
        star_vmag_limit: float | None,
        delta_t: timedelta,
        sky_disc_alpha: float,
        sky_disc_base_size: int,
        debug_eclipses: bool,
    ) -> None:
        try:
            now = datetime.now(timezone.utc) + delta_t
            time_obj = astropy.time.Time(now)

            stars, loc = calculate_visible_stars(
                star_catalog,
                lat,
                lon,
                time_obj,
                view_center,
                max_vmag=star_vmag_limit,
            )
            planets = calculate_planets(lat, lon, time_obj, view_center)
            celestial_equator_points = calculate_celestial_equator_points(loc, time_obj, view_center)
            ecliptic_points = calculate_ecliptic_points(loc, time_obj, view_center)
            horizon_points = calculate_horizon_points(view_center)
            celestial_data = CelestialData(
                time=time_obj,
                planets=planets,
                stars=stars,
                celestial_equator_points=celestial_equator_points,
                ecliptic_points=ecliptic_points,
                horizon_points=horizon_points,
            )

            if debug_eclipses:
                for body in planets:
                    if body.name == "sun" and body.solar_eclipse_info.is_eclipse:
                        logger.debug("Solar eclipse detected")
                    if body.name == "moon" and body.lunar_eclipse_info.is_eclipse:
                        logger.debug("Lunar eclipse detected")

            sun_altaz = None
            solar_eclipse_info = None
            for body in planets:
                if body.name == "sun":
                    sun_altaz = (body.alt, body.az)
                    solar_eclipse_info = body.solar_eclipse_info
                    break

            sky_disc_img: QImage | None = None
            if sky_disc_alpha > 0.0 and sun_altaz is not None:
                base = sky_disc_base_size
                fixed_geom = render_draw.get_screen_geometry(base, base, base // 2)
                ef = eclipse_factor_from_info(solar_eclipse_info)
                sky_disc_img = draw_sky_disc.draw_sky_color_disc(
                    fixed_geom,
                    view_center,
                    sun_altaz,
                    observer_lat_deg=lat,
                    alpha=sky_disc_alpha,
                    eclipse_factor=ef,
                    mask_fov_deg=93,
                )

            payload: Dict[str, object] = {"celestial": celestial_data, "sky_disc": sky_disc_img}
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
