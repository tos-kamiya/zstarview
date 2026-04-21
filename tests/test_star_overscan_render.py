from __future__ import annotations

import astropy.time
import numpy as np
from PySide6.QtGui import QImage, QPainter

from zstarview.render import stars as render_stars
from zstarview.types import CelestialData, ScreenGeometry, StarCatalogMeta, ViewerData
from zstarview.render.qt_image import qimage_to_np_rgba


def _single_star_celestial_data(
    *,
    alt: float,
    az: float,
    vmag: float = 5.5,
    bv: float = 0.0,
    size_factor: float = 0.18,
    color_factor_base: float = 0.7,
) -> CelestialData:
    return CelestialData(
        time=astropy.time.Time("2026-03-20T00:00:00", scale="utc"),
        planets=[],
        stars={
            "star_index": np.array([0], dtype=np.int32),
            "alt": np.array([alt], dtype=float),
            "az": np.array([az], dtype=float),
            "vmag": np.array([vmag], dtype=float),
            "bv": np.array([bv], dtype=float),
            "size_factor": np.array([size_factor], dtype=float),
            "color_factor_base": np.array([color_factor_base], dtype=float),
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
        star_catalog_meta=StarCatalogMeta(
            name_indices=np.array([0], dtype=np.int32),
            names=np.array(["Faint"], dtype=object),
            source_id_indices=np.array([0], dtype=np.int32),
            source_ids=np.array(["HIP_FAINT"], dtype=object),
        ),
    )


def test_draw_stars_keeps_faint_overscan_star_outside_90_deg_background() -> None:
    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        geometry = ScreenGeometry(center=(120, 120), radius=80)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            content_fov_deg=110.0,
        )
        celestial_data = _single_star_celestial_data(
            # Angular distance is 95 deg from the view center, so this is outside 90 deg
            # but still within the requested 110 deg content FOV.
            alt=-50.0,
            az=180.0,
        )

        render_stars.draw_stars(
            painter,
            geometry,
            celestial_data,
            viewer,
            star_base_radius=4.0,
            viewport_size=(240, 240),
            content_fov_deg=110.0,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    assert int(arr[199, 120, 3]) > 0


def test_draw_stars_uses_peak_channel_as_alpha_for_faint_pixels() -> None:
    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        geometry = ScreenGeometry(center=(120, 120), radius=80)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(90.0, 180.0),
            content_fov_deg=90.0,
        )
        celestial_data = _single_star_celestial_data(
            alt=90.0,
            az=180.0,
            vmag=5.8,
            bv=0.35,
            size_factor=0.18,
            color_factor_base=0.18,
        )

        render_stars.draw_stars(
            painter,
            geometry,
            celestial_data,
            viewer,
            star_base_radius=4.0,
            viewport_size=(240, 240),
            content_fov_deg=90.0,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    center = arr[120, 120, :]

    assert 0 < int(center[3]) < 255
    assert int(np.max(center[:3])) >= 250
    assert int(np.min(center[:3])) > 0
