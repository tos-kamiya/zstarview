from __future__ import annotations

import astropy.time
import numpy as np
from PySide6.QtGui import QImage, QPainter

from zstarview.render import draw as render_draw
from zstarview.types import CelestialData, ScreenGeometry, ViewerData
from zstarview.utils.qt import qimage_to_np_rgba


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
        celestial_data = CelestialData(
            time=astropy.time.Time("2026-03-20T00:00:00", scale="utc"),
            planets=[],
            stars={
                "name": np.array(["Faint"], dtype=object),
                "source_id": np.array(["HIP_FAINT"], dtype=object),
                # Angular distance is 95 deg from the view center, so this is outside 90 deg
                # but still within the requested 110 deg content FOV.
                "alt": np.array([-50.0], dtype=float),
                "az": np.array([180.0], dtype=float),
                "vmag": np.array([5.5], dtype=float),
                "bv": np.array([0.0], dtype=float),
                "size_factor": np.array([0.18], dtype=float),
                "color_factor_base": np.array([0.7], dtype=float),
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
        )

        render_draw.draw_stars(
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
    assert int(arr[204, 120, 3]) > 0
