import numpy as np
from PySide6.QtCore import QRectF
from PySide6.QtGui import QImage, QPainter

from zstarview.paths import PALETTE_NEVER_RISES_RGB
from zstarview.render.background import draw_radial_background, draw_window_border
from zstarview.render.geometry import get_screen_geometry
from zstarview.render.sky_disc import NEVER_RISES_TINT_RGB, draw_sky_color_disc, draw_uniform_sky_color_disc
from zstarview.types import ScreenGeometry
from zstarview.render.qt_image import qimage_to_np_rgba


def test_screen_geometry_wide_mode_top_is_always_tangent() -> None:
    width, height = 1000, 600
    g_zenith = get_screen_geometry(width, height, 90.0)
    g_mid = get_screen_geometry(width, height, 45.0)
    g_low = get_screen_geometry(width, height, 5.0)

    assert g_zenith.center[1] - g_zenith.radius == 0
    assert g_mid.center[1] - g_mid.radius == 0
    assert g_low.center[1] - g_low.radius == 0
    # Lower view altitude should increase radius in wide mode.
    assert g_low.radius > g_mid.radius > g_zenith.radius


def test_screen_geometry_uses_height_radius_for_extra_wide_mode() -> None:
    width, height = 1400, 600  # wider than 2:1
    g = get_screen_geometry(width, height, 45.0)
    expected_radius = int(height / (1.0 + 45.0 / 90.0))
    assert g.center == (width // 2, expected_radius)
    assert g.radius == expected_radius


def test_screen_geometry_tall_mode_stays_centered() -> None:
    width, height = 600, 1000
    g = get_screen_geometry(width, height, 30.0)

    assert g.center == (width // 2, height // 2)
    assert g.radius == min(width // 2, height // 2)


def test_sky_disc_raw_image_keeps_below_horizon_untinted() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    img = draw_sky_color_disc(
        geom,
        view_center=(0.0, 0.0),
        sun_altaz=(-90.0, 0.0),  # base sky becomes black in current model
        alpha=1.0,
        eclipse_factor=1.0,
        content_fov_deg=90.0,
    )
    arr = qimage_to_np_rgba(img)

    # Same x-position, one sample above and one below horizon.
    top_rgb = arr[20, 80, :3].astype(int)
    bottom_rgb = arr[140, 80, :3].astype(int)

    # Ground tint is applied later by the compositor, not in the raw sky disc.
    assert int(top_rgb.max()) <= 1
    assert int(bottom_rgb.max()) <= 1


def test_sky_disc_can_reduce_disc_opacity_without_changing_rgb_samples() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    img = draw_sky_color_disc(
        geom,
        view_center=(0.0, 0.0),
        sun_altaz=(-90.0, 0.0),
        alpha=1.0,
        disc_opacity=0.45,
        eclipse_factor=1.0,
        content_fov_deg=90.0,
    )
    arr = qimage_to_np_rgba(img)

    assert int(arr[20, 80, 3]) == 115
    assert int(arr[20, 80, :3].max()) <= 1


def test_sky_disc_ignores_never_rises_tint_by_observer_latitude() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    north_img = draw_sky_color_disc(
        geom,
        view_center=(10.0, 180.0),
        sun_altaz=(45.0, 180.0),
        observer_lat_deg=35.0,
        alpha=1.0,
        eclipse_factor=1.0,
        content_fov_deg=90.0,
    )
    south_img = draw_sky_color_disc(
        geom,
        view_center=(10.0, 180.0),
        sun_altaz=(45.0, 180.0),
        observer_lat_deg=-35.0,
        alpha=1.0,
        eclipse_factor=1.0,
        content_fov_deg=90.0,
    )

    north_arr = qimage_to_np_rgba(north_img)
    south_arr = qimage_to_np_rgba(south_img)

    # The raw sky disc should no longer inject the never-rises tint; that is
    # handled later by the compositor.
    assert np.array_equal(north_arr[120, 80, :3], south_arr[120, 80, :3])


def test_uniform_sky_disc_uses_single_disc_color() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    img = draw_uniform_sky_color_disc(geom, view_center=(0.0, 0.0), content_fov_deg=90.0)
    arr = qimage_to_np_rgba(img)

    center_rgb = arr[80, 80, :3].astype(int)
    top_rgb = arr[20, 80, :3].astype(int)
    lower_rgb = arr[120, 80, :3].astype(int)

    assert int(arr[80, 80, 3]) == 255
    assert np.array_equal(center_rgb, top_rgb)
    assert np.array_equal(center_rgb, np.array([10, 10, 10]))
    assert np.array_equal(center_rgb, lower_rgb)


def test_uniform_sky_disc_can_reduce_disc_opacity() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    img = draw_uniform_sky_color_disc(
        geom,
        view_center=(0.0, 0.0),
        disc_opacity=0.45,
        content_fov_deg=90.0,
    )
    arr = qimage_to_np_rgba(img)

    assert int(arr[80, 80, 3]) == 115
    assert np.all(np.abs(arr[80, 80, :3].astype(int) - np.array([10, 10, 10])) <= 1)


def test_uniform_sky_disc_content_fov_fills_corner_overscan_area() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)

    default_img = draw_uniform_sky_color_disc(geom, view_center=(90.0, 0.0), content_fov_deg=90.0)
    overscan_img = draw_uniform_sky_color_disc(geom, view_center=(90.0, 0.0), content_fov_deg=110.0)

    default_arr = qimage_to_np_rgba(default_img)
    overscan_arr = qimage_to_np_rgba(overscan_img)

    # This sample lies outside the 90-degree inscribed circle but inside a 110-degree square overscan region.
    assert int(default_arr[20, 20, 3]) == 0
    assert int(overscan_arr[20, 20, 3]) == 255


def test_never_rises_tint_uses_first_palette_swatch() -> None:
    expected = np.array(PALETTE_NEVER_RISES_RGB, dtype=np.float32) / 255.0
    assert np.allclose(NEVER_RISES_TINT_RGB, expected)


def test_radial_background_uses_black_inner_disc_for_all_main_themes() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=60)
    rect = QRectF(0.0, 0.0, 160.0, 160.0)

    for preset in ("white", "day", "night", "black", "transparent"):
        img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
        img.fill(0)
        painter = QPainter(img)
        draw_radial_background(painter, rect, geom, preset=preset)
        painter.end()

        arr = qimage_to_np_rgba(img)
        center_rgb = arr[80, 80, :3].astype(int)

        if preset == "transparent":
            assert int(center_rgb.max()) <= 5, preset
        else:
            assert np.array_equal(center_rgb, np.array([4, 4, 4])), preset


def test_radial_background_fades_between_content_fov_and_window_edge() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=60)
    rect = QRectF(0.0, 0.0, 160.0, 160.0)

    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_radial_background(painter, rect, geom, preset="night", content_fov_deg=100.0)
    painter.end()

    arr = qimage_to_np_rgba(img)

    # Around the content-FOV boundary there should still be visible alpha.
    assert int(arr[13, 80, 3]) > 0
    # Toward the window corner, the fade should become more transparent.
    assert int(arr[10, 10, 3]) < int(arr[13, 80, 3])


def test_radial_background_opaque_mode_keeps_full_alpha_at_edges() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=60)
    rect = QRectF(0.0, 0.0, 160.0, 160.0)

    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_radial_background(
        painter,
        rect,
        geom,
        preset="night",
        content_fov_deg=100.0,
        opaque=True,
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[80, 80, 3]) == 255
    assert int(arr[13, 80, 3]) == 255
    assert int(arr[10, 10, 3]) == 255


def test_window_frame_has_no_outer_border_for_black() -> None:
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        preset="black",
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[1, 80, 3]) == 0
    assert int(arr[80, 1, 3]) == 0
    assert 0 < int(arr[0, 145, 3]) < 255
    assert 0 < int(arr[20, 131, 3]) < 255
    assert int(arr[159, 159, 3]) == 0


def test_window_frame_draws_bottom_right_grip_line_inside_frame() -> None:
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        preset="white",
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[133, 151, 3]) > 0
    assert int(arr[136, 148, 3]) > 0
    assert int(arr[142, 142, 3]) > 0
    assert int(arr[150, 150, 3]) == 0
    assert int(arr[128, 128, 3]) == 0


def test_window_frame_draws_top_right_menu_square_inside_frame() -> None:
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        preset="white",
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert 0 < int(arr[0, 145, 3]) < 255
    assert 0 < int(arr[20, 131, 3]) < 255
    assert int(arr[20, 124, 3]) == 0


def test_window_frame_does_not_double_draw_under_menu_panel() -> None:
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        preset="white",
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[2, 145, 3]) == int(arr[24, 145, 3])


def test_window_frame_draws_hamburger_icon_lines() -> None:
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        preset="white",
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[9, 145, 3]) > 0
    assert int(arr[14, 145, 3]) > 0
    assert int(arr[19, 145, 3]) > 0


def test_window_frame_menu_panel_position_is_consistent_across_presets() -> None:
    night_img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    night_img.fill(0)
    night_painter = QPainter(night_img)
    draw_window_border(
        night_painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        preset="night",
    )
    night_painter.end()

    black_img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    black_img.fill(0)
    black_painter = QPainter(black_img)
    draw_window_border(
        black_painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        preset="black",
    )
    black_painter.end()

    night_arr = qimage_to_np_rgba(night_img)
    black_arr = qimage_to_np_rgba(black_img)

    assert 0 < int(night_arr[0, 145, 3]) < 255
    assert 0 < int(black_arr[0, 145, 3]) < 255
    assert 0 < int(night_arr[20, 131, 3]) < 255
    assert 0 < int(black_arr[20, 131, 3]) < 255
    assert int(night_arr[20, 124, 3]) == 0
    assert int(black_arr[20, 124, 3]) == 0


def test_transparent_window_frame_skips_border_but_keeps_menu_and_grip() -> None:
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        preset="transparent",
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[1, 80, 3]) == 0
    assert int(arr[80, 1, 3]) == 0
    assert 0 < int(arr[0, 145, 3]) < 255
    assert 0 < int(arr[20, 131, 3]) < 255
    assert int(arr[133, 151, 3]) > 0
    assert int(arr[136, 148, 3]) > 0
