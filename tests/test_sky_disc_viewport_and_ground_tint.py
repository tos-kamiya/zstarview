import numpy as np
from astropy.time import Time
from PySide6.QtCore import QRectF
from PySide6.QtGui import QColor, QImage, QPainter

from zstarview.astro import altaz_to_normalized_xy
from zstarview.gui.composite import (
    SkyCompositorCache,
    _apply_ground_reset,
    _clip_sky_image_to_disc,
    _dimalt_ring_color_for_sky_image,
)
from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.render.background import (
    atlas_background_tint_rgba,
    dimalt_ring_pen_color_from_color,
    draw_radial_background,
    draw_window_border,
)
from zstarview.render.geometry import get_screen_geometry
from zstarview.render.instrument_background import (
    draw_instrument_background,
    draw_instrument_time_of_day_marker,
)
from zstarview.render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba
from zstarview.render.sky_disc import (
    _apply_sky_intensity,
    _render_sky_color_disc_cached,
    draw_sky_color_disc,
    draw_uniform_sky_color_disc,
    sky_color_samples,
    sky_disc_update_interval,
)
from zstarview.types import ScreenGeometry


def test_sky_disc_refresh_interval_is_short_only_near_sunrise_or_sunset() -> None:
    assert sky_disc_update_interval(15.0) == 15
    assert sky_disc_update_interval(0.0) == 15
    assert sky_disc_update_interval(-15.0) == 15
    assert sky_disc_update_interval(15.01) == 60
    assert sky_disc_update_interval(-15.01) == 60
    assert sky_disc_update_interval(None) == 60


def test_atlas_background_tint_interpolates_day_to_night() -> None:
    day = atlas_background_tint_rgba(6.0)
    twilight = atlas_background_tint_rgba(0.0)
    night = atlas_background_tint_rgba(-6.0)

    assert day == (150, 200, 235, 255)
    assert twilight == (245, 168, 82, 255)
    assert night == (48, 52, 58, 255)
    assert atlas_background_tint_rgba(None) is None


def test_instrument_background_draws_atlas_time_marker_in_top_left() -> None:
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_instrument_background(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0).toRect(),
        theme=THEME_STYLES_BY_PRESET["atlas-white"],
    )
    draw_instrument_time_of_day_marker(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0).toRect(),
        sun_alt_deg=-20.0,
        bottom_left=False,
    )
    painter.end()

    pixels = qimage_to_np_rgba(img)
    assert np.array_equal(pixels[80, 80, :3], np.array([255, 255, 255]))
    assert np.array_equal(pixels[1, 1, :3], np.array([255, 255, 255]))
    assert np.array_equal(pixels[10, 10, :3], np.array([48, 52, 58]))
    assert np.array_equal(pixels[150, 10, :3], np.array([255, 255, 255]))


def test_sky_color_samples_keep_a_small_blue_night_floor() -> None:
    colors = sky_color_samples(
        np.asarray([45.0], dtype=np.float32),
        np.asarray([0.0], dtype=np.float32),
        (-90.0, 0.0),
    )

    np.testing.assert_allclose(
        colors[0],
        np.asarray([0.5, 1.0, 2.5], dtype=np.float32) / 255.0,
    )


def test_sky_intensity_favors_luminance_at_high_strength() -> None:
    colors = np.asarray([[0.15, 0.45, 0.90]], dtype=np.float32)
    weights = np.asarray([0.2126, 0.7152, 0.0722], dtype=np.float32)

    low = _apply_sky_intensity(colors, 0.2)[0]
    high = _apply_sky_intensity(colors, 1.0)[0]
    low_luminance = float(low @ weights)
    high_luminance = float(high @ weights)
    low_chroma = float(np.ptp(low))
    high_chroma = float(np.ptp(high))

    assert high_luminance > low_luminance
    assert high_chroma > low_chroma
    assert high_chroma / high_luminance < low_chroma / low_luminance
    source_luminance = float(colors[0] @ weights)
    expected_high = colors[0] + 0.18 * (source_luminance - colors[0])
    np.testing.assert_allclose(high, expected_high, rtol=1e-6)


def test_zero_sky_intensity_removes_scattered_color() -> None:
    colors = np.asarray([[0.15, 0.45, 0.90]], dtype=np.float32)

    np.testing.assert_array_equal(
        _apply_sky_intensity(colors, 0.0),
        np.zeros_like(colors),
    )


def test_sky_color_disc_lowers_ambient_floor_with_sun_altitude() -> None:
    geom = ScreenGeometry(center=(20, 20), radius=18)
    twilight = draw_sky_color_disc(
        geom,
        view_center=(45.0, 180.0),
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
        sun_altaz=(-9.0, 0.0),
        time_obj=Time("2026-07-30T11:00:00", format="isot", scale="utc"),
        timezone_name="Asia/Tokyo",
        alpha=1.0,
        disc_opacity=1.0,
        image_size=(40, 40),
    )
    astronomical_night = draw_sky_color_disc(
        geom,
        view_center=(45.0, 180.0),
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
        sun_altaz=(-18.0, 0.0),
        time_obj=Time("2026-07-30T17:00:00", format="isot", scale="utc"),
        timezone_name="Asia/Tokyo",
        alpha=1.0,
        disc_opacity=1.0,
        image_size=(40, 40),
    )

    twilight_pixels = qimage_to_np_rgba(twilight)
    astronomical_night_pixels = qimage_to_np_rgba(astronomical_night)
    assert np.all(twilight_pixels[20, 20, :3] >= astronomical_night_pixels[20, 20, :3])
    assert np.any(twilight_pixels[20, 20, :3] > astronomical_night_pixels[20, 20, :3])


def test_sky_disc_cache_keeps_only_recent_qimages() -> None:
    _render_sky_color_disc_cached.cache_clear()
    geom = ScreenGeometry(center=(20, 20), radius=18)

    for sun_alt in (5.0, 6.0, 7.0, 8.0):
        draw_sky_color_disc(
            geom,
            view_center=(45.0, 180.0),
            edge_fov_deg=90.0,
            content_fov_deg=90.0,
            sun_altaz=(sun_alt, 90.0),
            alpha=1.0,
            disc_opacity=1.0,
            image_size=(40, 40),
        )

    assert _render_sky_color_disc_cached.cache_info().maxsize == 2
    assert _render_sky_color_disc_cached.cache_info().currsize <= 2


def test_dimalt_ring_pen_color_darkens_sample_color() -> None:
    bright_sample = QColor(200, 180, 160, 123)
    bright_result = dimalt_ring_pen_color_from_color(bright_sample)

    dark_sample = QColor(20, 18, 16, 123)
    dark_result = dimalt_ring_pen_color_from_color(dark_sample)

    assert bright_result.alpha() == bright_sample.alpha()
    assert bright_result.red() < bright_sample.red()
    assert bright_result.green() < bright_sample.green()
    assert bright_result.blue() < bright_sample.blue()

    assert dark_result.alpha() == dark_sample.alpha()
    assert dark_result.red() > dark_sample.red()
    assert dark_result.green() > dark_sample.green()
    assert dark_result.blue() > dark_sample.blue()


def test_dimalt_ring_color_uses_alt_specific_samples() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    y, x = np.ogrid[:160, :160]
    dist = np.sqrt((x - 80.0) ** 2 + (y - 80.0) ** 2)
    radial = np.clip((dist / 80.0) * 220.0 + 20.0, 0.0, 255.0).astype(np.uint8)
    sky = np.zeros((160, 160, 4), dtype=np.uint8)
    sky[..., 0] = radial
    sky[..., 1] = radial
    sky[..., 2] = radial
    sky[..., 3] = 255
    img = np_rgba_to_qimage(sky)

    bright_ring = _dimalt_ring_color_for_sky_image(
        img,
        geom,
        (90.0, 180.0),
        alt_deg=30.0,
        edge_fov_deg=90.0,
    )
    dark_ring = _dimalt_ring_color_for_sky_image(
        img,
        geom,
        (90.0, 180.0),
        alt_deg=60.0,
        edge_fov_deg=90.0,
    )

    assert bright_ring is not None
    assert dark_ring is not None
    assert bright_ring.red() > dark_ring.red()


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
    edge_fov = 90.0
    content_fov = 110.0
    g = get_screen_geometry(
        width,
        height,
        30.0,
        edge_fov_deg=edge_fov,
        content_fov_deg=content_fov,
    )
    edge_radius = height / 2.0
    content_fit_radius = edge_radius * edge_fov / content_fov
    blend = ((width / height) - 0.5) / 0.5
    expected_radius = int(
        content_fit_radius + (edge_radius - content_fit_radius) * blend
    )

    assert g.center == (width // 2, height // 2)
    assert g.radius == expected_radius
    assert g.radius > width // 2


def test_screen_geometry_tall_mode_near_square_uses_edge_fov_height_radius() -> None:
    width, height = 999, 1000
    g = get_screen_geometry(
        width,
        height,
        30.0,
        edge_fov_deg=90.0,
        content_fov_deg=110.0,
    )

    assert g.center == (width // 2, height // 2)
    assert g.radius == 499


def test_screen_geometry_tall_mode_fits_content_fov_at_one_by_two() -> None:
    width, height = 500, 1000
    edge_fov = 90.0
    content_fov = 110.0
    g = get_screen_geometry(
        width,
        height,
        30.0,
        edge_fov_deg=edge_fov,
        content_fov_deg=content_fov,
    )

    assert g.center == (width // 2, height // 2)
    assert abs((g.radius * content_fov / edge_fov) - (height / 2.0)) <= 1.0


def test_screen_geometry_tall_mode_linearly_blends_between_one_and_one_by_two() -> None:
    width, height = 750, 1000
    edge_fov = 90.0
    content_fov = 110.0
    g = get_screen_geometry(
        width,
        height,
        30.0,
        edge_fov_deg=edge_fov,
        content_fov_deg=content_fov,
    )
    edge_radius = height / 2.0
    content_fit_radius = edge_radius * edge_fov / content_fov
    expected_radius = int((edge_radius + content_fit_radius) / 2.0)

    assert g.center == (width // 2, height // 2)
    assert g.radius == expected_radius


def test_screen_geometry_tall_mode_keeps_content_fit_below_one_by_two() -> None:
    edge_fov = 90.0
    content_fov = 110.0
    g_one_by_two = get_screen_geometry(
        500,
        1000,
        30.0,
        edge_fov_deg=edge_fov,
        content_fov_deg=content_fov,
    )
    g_taller = get_screen_geometry(
        400,
        1000,
        30.0,
        edge_fov_deg=edge_fov,
        content_fov_deg=content_fov,
    )

    assert g_taller.center == (200, 500)
    assert g_taller.radius == g_one_by_two.radius


def test_screen_geometry_tall_mode_same_edge_and_content_fov_uses_height_radius() -> None:
    g = get_screen_geometry(
        600,
        1000,
        30.0,
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
    )

    assert g.center == (300, 500)
    assert g.radius == 500


def test_sky_disc_raw_image_keeps_below_horizon_untinted() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    img = draw_sky_color_disc(
        geom,
        view_center=(0.0, 0.0),
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
        sun_altaz=(-90.0, 0.0),  # base sky becomes black in current model
        alpha=1.0,
        eclipse_factor=1.0,
    )
    arr = qimage_to_np_rgba(img)

    # Same x-position, one sample above and one below horizon.
    top_rgb = arr[20, 80, :3].astype(int)
    bottom_rgb = arr[140, 80, :3].astype(int)

    # Ground tint is applied later by the compositor, not in the raw sky disc.
    assert int(top_rgb.max()) <= 20
    assert int(bottom_rgb.max()) <= 20
    assert np.array_equal(top_rgb, bottom_rgb)


def test_sky_disc_can_reduce_disc_opacity_without_changing_rgb_samples() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    img = draw_sky_color_disc(
        geom,
        view_center=(0.0, 0.0),
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
        sun_altaz=(-90.0, 0.0),
        alpha=1.0,
        disc_opacity=0.45,
        eclipse_factor=1.0,
    )
    arr = qimage_to_np_rgba(img)

    assert int(arr[20, 80, 3]) == 115
    assert int(arr[20, 80, :3].max()) <= 20


def test_uniform_sky_disc_uses_single_disc_color() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    img = draw_uniform_sky_color_disc(geom, view_center=(0.0, 0.0), edge_fov_deg=90.0, content_fov_deg=90.0)
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
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
        disc_opacity=0.45,
    )
    arr = qimage_to_np_rgba(img)

    assert int(arr[80, 80, 3]) == 115
    assert np.all(np.abs(arr[80, 80, :3].astype(int) - np.array([10, 10, 10])) <= 1)


def test_uniform_sky_disc_content_fov_fills_corner_overscan_area() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)

    default_img = draw_uniform_sky_color_disc(geom, view_center=(90.0, 0.0), edge_fov_deg=90.0, content_fov_deg=90.0)
    overscan_img = draw_uniform_sky_color_disc(geom, view_center=(90.0, 0.0), edge_fov_deg=90.0, content_fov_deg=110.0)

    default_arr = qimage_to_np_rgba(default_img)
    overscan_arr = qimage_to_np_rgba(overscan_img)

    # This sample lies outside the 90-degree inscribed circle but inside a 110-degree square overscan region.
    assert int(default_arr[20, 20, 3]) == 0
    assert int(overscan_arr[20, 20, 3]) == 255


def test_compositor_sky_disc_clip_antialiases_final_resolution_edge() -> None:
    source = QImage(100, 100, QImage.Format.Format_ARGB32_Premultiplied)
    source.fill(0xFFFFFFFF)

    clipped = _clip_sky_image_to_disc(
        source,
        geometry=ScreenGeometry(center=(50, 50), radius=40),
        viewport_offset=(0, 0),
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
    )
    alpha = qimage_to_np_rgba(clipped)[..., 3]

    assert int(alpha[9, 50]) == 0
    assert 0 < int(alpha[10, 50]) < 255
    assert int(alpha[50, 50]) == 255


def test_altitude_rings_dim_sky_disc_before_compositing() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    sky = np.zeros((160, 160, 4), dtype=np.uint8)
    sky[..., :3] = 90
    sky[..., 3] = 255

    compositor = SkyCompositorCache()
    canvas = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    try:
        compositor.draw(
            painter,
            geom,
            np_rgba_to_qimage(sky),
            cloud_alpha=0.0,
            view_center=(90.0, 180.0),
            theme=THEME_STYLES_BY_PRESET["white"],
            edge_fov_deg=90.0,
            content_fov_deg=110.0,
            sky_disc_altaz_rings="dimalt",
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(canvas)
    nx, ny = altaz_to_normalized_xy(
        30.0,
        180.0,
        (90.0, 180.0),
        edge_fov_deg=90.0,
    )
    x = int(round(geom.center[0] + (nx * geom.radius)))
    y = int(round(geom.center[1] + (ny * geom.radius)))
    y0 = max(0, y - 1)
    y1 = min(arr.shape[0], y + 2)
    x0 = max(0, x - 1)
    x1 = min(arr.shape[1], x + 2)
    assert float(arr[y0:y1, x0:x1, :3].mean()) < 90.0


def test_compositor_can_skip_sky_disc_without_black_fallback() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    canvas = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    try:
        SkyCompositorCache().draw(
            painter,
            geom,
            QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied),
            cloud_alpha=1.0,
            ground_reset_rgba=(4, 4, 4, 255),
            draw_sky_disc=False,
            content_fov_deg=90.0,
        )
    finally:
        painter.end()

    assert not np.any(qimage_to_np_rgba(canvas)[..., 3])


def test_ground_reset_keeps_transparent_theme_below_horizon_semi_transparent() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=80)
    base = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    base.fill(0)

    result = _apply_ground_reset(
        base,
        geometry=geom,
        view_center=(0.0, 0.0),
        terrain_profile_altaz=None,
        ground_reset_rgba=(4, 4, 4, 102),
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
    )
    arr = qimage_to_np_rgba(result)

    assert int(arr[80, 80, 3]) == 0
    assert int(arr[140, 80, 3]) == 102
    assert np.array_equal(arr[140, 80, :3], np.array([4, 4, 4]))


def test_radial_background_uses_black_inner_disc_for_all_main_themes() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=60)
    rect = QRectF(0.0, 0.0, 160.0, 160.0)

    for preset in ("white", "day", "night", "black", "transparent"):
        theme = THEME_STYLES_BY_PRESET[preset]
        img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
        img.fill(0)
        painter = QPainter(img)
        draw_radial_background(painter, rect, geom, theme=theme)
        painter.end()

        arr = qimage_to_np_rgba(img)
        center_rgb = arr[80, 80, :3].astype(int)

        if preset == "transparent":
            assert int(center_rgb.max()) <= 5, preset
            assert np.array_equal(arr[80, 80], arr[20, 20]), preset
        else:
            assert np.array_equal(center_rgb, np.array([4, 4, 4])), preset


def test_radial_background_fades_between_content_fov_and_window_edge() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=60)
    rect = QRectF(0.0, 0.0, 160.0, 160.0)

    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_radial_background(
        painter,
        rect,
        geom,
        theme=THEME_STYLES_BY_PRESET["night"],
        content_fov_deg=100.0,
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    # Around the content-FOV boundary there should still be visible alpha.
    assert int(arr[13, 80, 3]) > 0
    # Toward the window corner, the fade should become more transparent.
    assert int(arr[10, 10, 3]) < int(arr[30, 30, 3])


def test_radial_background_alt_rings_dim_background() -> None:
    geom = ScreenGeometry(center=(80, 80), radius=60)
    rect = QRectF(0.0, 0.0, 160.0, 160.0)

    base_img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    base_img.fill(0)
    base_painter = QPainter(base_img)
    draw_radial_background(
        base_painter,
        rect,
        geom,
        theme=THEME_STYLES_BY_PRESET["night"],
        view_center=(90.0, 180.0),
        altaz_rings_mode="off",
    )
    base_painter.end()

    ring_img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    ring_img.fill(0)
    ring_painter = QPainter(ring_img)
    draw_radial_background(
        ring_painter,
        rect,
        geom,
        theme=THEME_STYLES_BY_PRESET["night"],
        view_center=(90.0, 180.0),
        altaz_rings_mode="dimalt",
    )
    ring_painter.end()

    base_arr = qimage_to_np_rgba(base_img)
    ring_arr = qimage_to_np_rgba(ring_img)
    nx, ny = altaz_to_normalized_xy(
        30.0,
        180.0,
        (90.0, 180.0),
        edge_fov_deg=90.0,
    )
    x = int(round(geom.center[0] + (nx * geom.radius)))
    y = int(round(geom.center[1] + (ny * geom.radius)))
    y0 = max(0, y - 1)
    y1 = min(ring_arr.shape[0], y + 2)
    x0 = max(0, x - 1)
    x1 = min(ring_arr.shape[1], x + 2)
    assert float(ring_arr[y0:y1, x0:x1, :3].mean()) > float(base_arr[y0:y1, x0:x1, :3].mean())
    assert int(ring_arr[y0:y1, x0:x1, 3].max()) >= int(base_arr[y0:y1, x0:x1, 3].max())


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
        theme=THEME_STYLES_BY_PRESET["night"],
        content_fov_deg=100.0,
        opaque=True,
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[80, 80, 3]) == 255
    assert int(arr[13, 80, 3]) == 255
    assert int(arr[10, 10, 3]) == 255


def test_window_frame_does_not_draw_outer_border_for_black() -> None:
    theme = THEME_STYLES_BY_PRESET["black"]
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        theme=theme,
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[1, 80, 3]) == 0
    assert int(arr[80, 1, 3]) == 0
    assert int(arr[159, 159, 3]) == 0


def test_window_frame_does_not_draw_outer_border_for_white_and_day() -> None:
    for preset in ("white", "day"):
        theme = THEME_STYLES_BY_PRESET[preset]
        img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
        img.fill(0)
        painter = QPainter(img)
        draw_window_border(
            painter,
            QRectF(0.0, 0.0, 160.0, 160.0),
            theme=theme,
        )
        painter.end()

        arr = qimage_to_np_rgba(img)

        assert int(arr[0, 80, 3]) == 0, preset
        assert int(arr[80, 0, 3]) == 0, preset
        assert int(arr[159, 159, 3]) == 0, preset


def test_window_border_does_not_draw_resize_grip() -> None:
    theme = THEME_STYLES_BY_PRESET["white"]
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        theme=theme,
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[133, 151, 3]) == 0
    assert int(arr[136, 148, 3]) == 0
    assert int(arr[142, 142, 3]) == 0
    assert int(arr[128, 128, 3]) == 0


def test_window_frame_draws_top_right_menu_square_inside_frame() -> None:
    theme = THEME_STYLES_BY_PRESET["white"]
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        theme=theme,
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert 0 < int(arr[0, 145, 3]) < 255
    assert 0 < int(arr[20, 131, 3]) < 255


def test_window_frame_does_not_double_draw_under_menu_panel() -> None:
    theme = THEME_STYLES_BY_PRESET["white"]
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        theme=theme,
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[2, 145, 3]) == int(arr[24, 145, 3])


def test_window_frame_draws_hamburger_icon_lines() -> None:
    theme = THEME_STYLES_BY_PRESET["white"]
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        theme=theme,
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[9, 145, 3]) > 0
    assert int(arr[14, 145, 3]) > 0
    assert int(arr[19, 145, 3]) > 0


def test_window_frame_menu_panel_position_is_consistent_across_presets() -> None:
    night_theme = THEME_STYLES_BY_PRESET["night"]
    black_theme = THEME_STYLES_BY_PRESET["black"]
    night_img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    night_img.fill(0)
    night_painter = QPainter(night_img)
    draw_window_border(
        night_painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        theme=night_theme,
    )
    night_painter.end()

    black_img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    black_img.fill(0)
    black_painter = QPainter(black_img)
    draw_window_border(
        black_painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        theme=black_theme,
    )
    black_painter.end()

    night_arr = qimage_to_np_rgba(night_img)
    black_arr = qimage_to_np_rgba(black_img)

    assert 0 < int(night_arr[0, 145, 3]) < 255
    assert 0 < int(black_arr[0, 145, 3]) < 255
    assert 0 < int(night_arr[20, 131, 3]) < 255
    assert 0 < int(black_arr[20, 131, 3]) < 255


def test_transparent_window_frame_skips_border_and_keeps_menu() -> None:
    theme = THEME_STYLES_BY_PRESET["transparent"]
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        theme=theme,
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[1, 80, 3]) == 0
    assert int(arr[80, 1, 3]) == 0
    assert 0 < int(arr[0, 145, 3]) < 255
    assert 0 < int(arr[20, 131, 3]) < 255


def test_window_border_does_not_draw_resize_grip_for_transparent_theme() -> None:
    theme = THEME_STYLES_BY_PRESET["transparent"]
    img = QImage(160, 160, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_window_border(
        painter,
        QRectF(0.0, 0.0, 160.0, 160.0),
        theme=theme,
    )
    painter.end()

    arr = qimage_to_np_rgba(img)

    assert int(arr[142, 142, 3]) == 0
    assert int(arr[136, 148, 3]) == 0
