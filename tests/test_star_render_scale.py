from zstarview.ui.window import compute_star_render_surface_size


def test_render_surface_size_is_native_below_or_at_expected_width() -> None:
    assert compute_star_render_surface_size(400, 300, 600) == (400, 300)
    assert compute_star_render_surface_size(600, 450, 600) == (600, 450)


def test_render_surface_size_follows_inverse_sqrt_scale_above_expected_width() -> None:
    assert compute_star_render_surface_size(1200, 600, 600) == (849, 424)
    assert compute_star_render_surface_size(2400, 1200, 600) == (1200, 600)
    assert compute_star_render_surface_size(5400, 2700, 600) == (1800, 900)
