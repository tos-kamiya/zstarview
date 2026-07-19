# ruff: noqa: F401

from __future__ import annotations

import math

import astropy.time
import numpy as np
from PySide6.QtCore import QPoint, QPointF, QRectF
from PySide6.QtGui import QColor, QFont, QFontMetrics, QImage, QPainter
from PySide6.QtWidgets import QApplication

from zstarview.aircraft.types import AircraftOverlayPoint
from zstarview.paths import (
    ATLAS_THEME_PRESET,
    PALETTE_AIRCRAFT_AND_SATELLITE_RGB,
    THEME_STYLES_BY_PRESET,
)
from zstarview.render import aircraft as render_aircraft
from zstarview.render import background as render_background
from zstarview.render import guides as render_guides
from zstarview.render import overlay_info as render_overlay_info
from zstarview.render import satellites as render_satellites
from zstarview.render import solar_system as render_solar_system
from zstarview.render import stars as render_stars
from zstarview.render import text as render_text
from zstarview.render.deep_sky_objects import DSO_LABEL_TEXT_RGB
from zstarview.satellite_constants import SATELLITE_HORIZONS_CACHE_KEY
from zstarview.satellites.types import SatelliteOverlayPoint
from zstarview.types import (
    CelestialData,
    PlanetBody,
    ScreenGeometry,
    StarCatalogMeta,
    ViewerData,
)


def _empty_star_catalog_meta() -> StarCatalogMeta:
    return StarCatalogMeta(
        name_indices=np.array([], dtype=np.int32),
        names=np.array([], dtype=object),
        source_id_indices=np.array([], dtype=np.int32),
        source_ids=np.array([], dtype=object),
    )


def _empty_celestial_data(planets: list[PlanetBody]) -> CelestialData:
    return CelestialData(
        time=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
        planets=planets,
        stars={
            "star_index": np.array([], dtype=np.int32),
            "alt": np.array([], dtype=float),
            "az": np.array([], dtype=float),
            "vmag": np.array([], dtype=float),
            "bv": np.array([], dtype=float),
            "size_factor": np.array([], dtype=float),
            "color_factor_base": np.array([], dtype=float),
        },
        deep_sky_objects={
            "id": np.array([], dtype=object),
            "name": np.array([], dtype=object),
            "type": np.array([], dtype=object),
            "alt": np.array([], dtype=float),
            "az": np.array([], dtype=float),
            "vmag": np.array([], dtype=float),
            "major_arcmin": np.array([], dtype=float),
            "minor_arcmin": np.array([], dtype=float),
            "pa_deg": np.array([], dtype=float),
        },
        celestial_equator_points=[],
        ecliptic_points=[],
        horizon_points=[],
        star_catalog_meta=_empty_star_catalog_meta(),
    )


def _celestial_data_with_stars(
    stars: dict[str, np.ndarray], star_catalog_meta: StarCatalogMeta | None = None
) -> CelestialData:
    return CelestialData(
        time=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
        planets=[],
        stars=stars,
        deep_sky_objects={
            "id": np.array([], dtype=object),
            "name": np.array([], dtype=object),
            "type": np.array([], dtype=object),
            "alt": np.array([], dtype=float),
            "az": np.array([], dtype=float),
            "vmag": np.array([], dtype=float),
            "major_arcmin": np.array([], dtype=float),
            "minor_arcmin": np.array([], dtype=float),
            "pa_deg": np.array([], dtype=float),
        },
        celestial_equator_points=[],
        ecliptic_points=[],
        horizon_points=[],
        star_catalog_meta=star_catalog_meta,
    )


def _star_table(
    names: list[str],
    *,
    source_ids: list[str] | None = None,
    alt: float = 45.0,
    az: float = 180.0,
) -> tuple[dict[str, np.ndarray], StarCatalogMeta]:
    count = len(names)
    stars = {
        "star_index": np.arange(count, dtype=np.int32),
        "alt": np.full(count, alt, dtype=float),
        "az": np.full(count, az, dtype=float),
        "vmag": np.linspace(1.0, 5.0, count, dtype=float),
        "bv": np.zeros(count, dtype=float),
        "size_factor": np.ones(count, dtype=float),
        "color_factor_base": np.ones(count, dtype=float),
    }
    source_values = [""] * count if source_ids is None else source_ids
    name_indices = np.array(
        [idx for idx, name in enumerate(names) if str(name).strip()], dtype=np.int32
    )
    source_id_indices = np.array(
        [idx for idx, value in enumerate(source_values) if str(value).strip()],
        dtype=np.int32,
    )
    meta = StarCatalogMeta(
        name_indices=name_indices,
        names=np.array([names[idx] for idx in name_indices], dtype=object),
        source_id_indices=source_id_indices,
        source_ids=np.array(
            [source_values[idx] for idx in source_id_indices], dtype=object
        ),
    )
    return stars, meta


__all__ = [name for name in globals() if not name.startswith("__")]
