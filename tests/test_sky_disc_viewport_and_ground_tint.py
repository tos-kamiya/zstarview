from zstarview.render.draw import get_screen_geometry
from zstarview.render.draw_sky_disc import draw_sky_color_disc
from zstarview.types import ScreenGeometry
from zstarview.utils.qt import qimage_to_np_rgba


def test_screen_geometry_is_fixed_circle_independent_of_view_altitude() -> None:
    width, height = 1000, 600
    g_zenith = get_screen_geometry(width, height, 90.0)
    g_low = get_screen_geometry(width, height, 5.0)

    expected_center = (width // 2, height // 2)
    expected_radius = min((width - 20) // 2, (height - 20) // 2)

    assert g_zenith.center == expected_center
    assert g_low.center == expected_center
    assert g_zenith.radius == expected_radius
    assert g_low.radius == expected_radius


def test_sky_disc_adds_tint_below_horizon() -> None:
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

    # Above horizon remains near black; below horizon gets blue-gray tint.
    assert int(top_rgb.max()) <= 1
    assert int(bottom_rgb.max()) >= 10
