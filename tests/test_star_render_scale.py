from zstarview.gui.window import compute_star_render_surface_size


def test_render_surface_size_is_native_below_or_at_expected_width() -> None:
    assert compute_star_render_surface_size(400, 300, 400, 600) == (400, 300)
    assert compute_star_render_surface_size(600, 450, 600, 600) == (600, 450)


def test_render_surface_size_uses_quantized_disc_width_steps() -> None:
    # Smooth sqrt scaling above threshold.
    assert compute_star_render_surface_size(900, 450, 699, 600) == (834, 417)
    assert compute_star_render_surface_size(900, 450, 700, 600) == (833, 417)
    assert compute_star_render_surface_size(900, 450, 799, 600) == (780, 390)
    # Typical larger discs.
    assert compute_star_render_surface_size(1200, 600, 1200, 600) == (849, 424)
    assert compute_star_render_surface_size(1300, 650, 1299, 600) == (884, 442)
