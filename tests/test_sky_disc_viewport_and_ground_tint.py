import numpy as np
from PySide6.QtCore import QRectF
from PySide6.QtGui import QImage, QPainter

from zstarview.render.draw import draw_radial_background, get_screen_geometry
from zstarview.render.draw_sky_disc import draw_sky_color_disc, draw_uniform_sky_color_disc
from zstarview.types import ScreenGeometry
from zstarview.utils.qt import qimage_to_np_rgba


def test_screen_geometry_wide_mode_top_is_always_tangent() -> None:
    width, height = 1000, 600
    g_zenith = get_screen_geometry(width, height, 90.0)
    g_mid = get_screen_geometry(width, height, 45.0)
    g_low = get_screen_geometry(width, height, 5.0)

    assert g_zenith.center[1] - g_zenith.radius == 10
    assert g_mid.center[1] - g_mid.radius == 10
    assert g_low.center[1] - g_low.radius == 10
    # Lower view altitude should increase radius in wide mode.
    assert g_low.radius > g_mid.radius > g_zenith.radius


def test_screen_geometry_uses_height_radius_for_extra_wide_mode() -> None:
    width, height = 1400, 600  # wider than 2:1
    g = get_screen_geometry(width, height, 45.0)
    expected_radius = int((height - 20) / (1.0 + 45.0 / 90.0))
    assert g.center == (width // 2, 10 + expected_radius)
    assert g.radius == expected_radius


def test_screen_geometry_tall_mode_stays_centered() -> None:
    width, height = 600, 1000
    g = get_screen_geometry(width, height, 30.0)

    assert g.center == (width // 2, height // 2)
    assert g.radius == min((width - 20) // 2, (height - 20) // 2)


def test_sky_disc_raw_image_keeps_below_horizon_untinted() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    img = draw_sky_color_disc(
        geom,
        view_center=(0.0, 0.0),
        sun_altaz=(-90.0, 0.0),  # base sky becomes black in current model
        alpha=1.0,
        eclipse_factor=1.0,
    )
    arr = qimage_to_np_rgba(img)

    # Same x-position, one sample above and one below horizon.
    top_rgb = arr[20, 80, :3].astype(int)
    bottom_rgb = arr[140, 80, :3].astype(int)

    # Ground tint is applied later by the compositor, not in the raw sky disc.
    assert int(top_rgb.max()) <= 1
    assert int(bottom_rgb.max()) <= 1


def test_uniform_sky_disc_uses_single_disc_color() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    img = draw_uniform_sky_color_disc(geom, view_center=(0.0, 0.0))
    arr = qimage_to_np_rgba(img)

    center_rgb = arr[80, 80, :3].astype(int)
    top_rgb = arr[20, 80, :3].astype(int)
    lower_rgb = arr[120, 80, :3].astype(int)

    assert int(arr[80, 80, 3]) == 255
    assert np.array_equal(center_rgb, top_rgb)
    assert np.array_equal(center_rgb, np.array([10, 10, 10]))
    assert np.array_equal(center_rgb, lower_rgb)


def test_radial_background_uses_black_inner_disc_for_all_main_themes() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=60)
    rect = QRectF(0.0, 0.0, 160.0, 160.0)

    for preset in ("white", "day", "night", "black"):
        img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
        img.fill(0)
        painter = QPainter(img)
        draw_radial_background(painter, rect, geom, preset=preset)
        painter.end()

        arr = qimage_to_np_rgba(img)
        center_rgb = arr[80, 80, :3].astype(int)

        assert np.array_equal(center_rgb, np.array([4, 4, 4])), preset
