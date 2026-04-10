import numpy as np
from PySide6.QtGui import QImage, QPainter

from zstarview.render.earth_guide import (
    _effective_visible_altitude_limit_deg,
    draw_earth_guide,
    load_earth_guide_rings,
)
from zstarview.render.qt_image import qimage_to_np_rgba
from zstarview.types import ScreenGeometry


def test_load_earth_guide_rings_has_expected_runtime_payload() -> None:
    rings = load_earth_guide_rings()

    assert len(rings) >= 10
    assert sum(len(ring.points_lonlat_deg) for ring in rings) > 200
    assert all(ring.points_xyz.shape[1] == 3 for ring in rings)


def test_draw_earth_guide_renders_visible_lines_below_horizon() -> None:
    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    draw_earth_guide(
        painter,
        geometry=ScreenGeometry(center=(120, 120), radius=100),
        view_center=(0.0, 180.0),
        observer_lat_deg=35.68,
        observer_lon_deg=139.76,
        observer_height_m=635.0,
        terrain_profile_altaz=None,
        content_fov_deg=100.0,
    )
    painter.end()

    arr = qimage_to_np_rgba(image)
    alpha = arr[..., 3]

    assert int(np.count_nonzero(alpha)) > 50
    top_half = int(np.count_nonzero(alpha[:120, :]))
    bottom_half = int(np.count_nonzero(alpha[120:, :]))
    assert bottom_half > top_half


def test_earth_guide_uses_minus_five_degree_upper_clip() -> None:
    assert _effective_visible_altitude_limit_deg(180.0, None) == -5.0
    assert _effective_visible_altitude_limit_deg(180.0, [(12.0, 180.0)]) == 12.0
