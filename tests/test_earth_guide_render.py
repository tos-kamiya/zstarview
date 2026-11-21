import math

import numpy as np
import pytest
from PySide6.QtGui import QImage, QPainter

import zstarview.gui.composite as render_composite
from zstarview.paths import EARTH_GUIDE_LINE_COLOR, EARTH_GUIDE_OCEAN_COLOR
from zstarview.render.earth_guide import (
    EARTH_GUIDE_DAY_MIN_LEVEL,
    EARTH_GUIDE_NEAR_CUTOFF_KM,
    EARTH_GUIDE_NEAR_FADE_KM,
    EARTH_GUIDE_NIGHT_OUTLINE_FULL_DEG,
    EARTH_GUIDE_NIGHT_LEVEL,
    EARTH_GUIDE_TWILIGHT_LIMIT_DEG,
    EarthGuideRing,
    _observer_basis,
    altaz_to_earth_fixed,
    _is_seam_segment,
    build_land_mask,
    compute_night_outline_segments,
    draw_earth_guide,
    earth_guide_near_fade,
    earth_guide_night_outline_weight,
    earth_guide_peak_alpha,
    earth_guide_sunlight_level,
    far_surface_points,
    load_coastline_segments,
    load_earth_guide_rings,
    load_land_mask,
    render_earth_guide_rgba,
    sample_land_fraction,
    visible_altitude_limit_deg,
)
from zstarview.render.qt_image import qimage_to_np_rgba
from zstarview.types import ScreenGeometry, ViewerData


def _viewer(
    *,
    view_center: tuple[float, float] = (-60.0, 180.0),
    observer_lat_deg: float = 35.68,
    observer_lon_deg: float = 139.76,
    observer_height_m: float = 10.0,
    edge_fov_deg: float = 95.0,
    content_fov_deg: float = 100.0,
) -> ViewerData:
    return ViewerData(
        location=(observer_lat_deg, observer_lon_deg),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=view_center,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
        observer_height_m=observer_height_m,
    )


def _square_ring(lon0: float, lat0: float, lon1: float, lat1: float) -> EarthGuideRing:
    lonlat = np.asarray(
        [(lon0, lat0), (lon1, lat0), (lon1, lat1), (lon0, lat1)], dtype=np.float64
    )
    return EarthGuideRing(
        source_name="square",
        label_name=None,
        points_lonlat_deg=lonlat,
        points_xyz=np.zeros((4, 3), dtype=np.float64),
    )


def test_load_earth_guide_rings_has_expected_runtime_payload() -> None:
    rings = load_earth_guide_rings()

    assert len(rings) >= 10
    assert sum(len(ring.points_lonlat_deg) for ring in rings) > 200
    assert all(ring.points_xyz.shape[1] == 3 for ring in rings)


def test_load_earth_guide_lods_preserves_increasing_detail() -> None:
    rings_by_lod = [load_earth_guide_rings(lod) for lod in range(4)]
    point_counts = [
        sum(len(ring.points_lonlat_deg) for ring in rings) for rings in rings_by_lod
    ]

    assert point_counts == sorted(point_counts)
    assert point_counts[-1] > point_counts[0]


def test_build_land_mask_fills_polygon_interior() -> None:
    mask = build_land_mask((_square_ring(10.0, 10.0, 30.0, 20.0),), res_deg=1.0)

    assert mask.shape == (180, 360)
    inside = sample_land_fraction(mask, np.array([20.0]), np.array([15.0]))
    outside = sample_land_fraction(mask, np.array([-20.0]), np.array([15.0]))
    assert inside[0] == pytest.approx(1.0)
    assert outside[0] == pytest.approx(0.0)


def test_sample_land_fraction_wraps_longitude() -> None:
    mask = build_land_mask((_square_ring(170.0, -5.0, 180.0, 5.0),), res_deg=1.0)

    east = sample_land_fraction(mask, np.array([179.9]), np.array([0.0]))
    wrapped = sample_land_fraction(mask, np.array([-180.1]), np.array([0.0]))
    assert east[0] == pytest.approx(wrapped[0])
    assert east[0] > 0.5


@pytest.mark.parametrize(
    ("lon", "lat", "expected_land"),
    [
        (20.0, 0.0, True),  # Africa
        (-60.0, -10.0, True),  # South America
        (135.0, -25.0, True),  # Australia
        (0.0, -90.0, True),  # Antarctica
        (-150.0, 0.0, False),  # Pacific
        (-30.0, 30.0, False),  # Atlantic
        (80.0, -20.0, False),  # Indian Ocean
    ],
)
def test_bundled_land_mask_matches_known_places(
    lon: float, lat: float, expected_land: bool
) -> None:
    value = sample_land_fraction(load_land_mask(), np.array([lon]), np.array([lat]))[0]

    assert bool(value > 0.5) is expected_land


def test_sunlight_level_day_twilight_and_night() -> None:
    def level(h_deg: float) -> float:
        return float(earth_guide_sunlight_level(np.array([math.sin(math.radians(h_deg))]))[0])

    assert level(90.0) == pytest.approx(1.0)
    assert level(0.0) == pytest.approx(EARTH_GUIDE_DAY_MIN_LEVEL)
    assert level(EARTH_GUIDE_TWILIGHT_LIMIT_DEG) == pytest.approx(EARTH_GUIDE_NIGHT_LEVEL)
    assert level(-60.0) == pytest.approx(EARTH_GUIDE_NIGHT_LEVEL)
    samples = [level(h) for h in (-20.0, -12.0, -6.0, -1.0, 0.0, 30.0, 60.0)]
    assert samples == sorted(samples)


def test_sun_direction_matches_observer_solar_altitude() -> None:
    _origin, east, north, up = _observer_basis(35.0, 139.0, 0.0)

    for sun_alt in (-30.0, 0.0, 45.0):
        sun_dir = altaz_to_earth_fixed(sun_alt, 123.0, east=east, north=north, up=up)
        assert float(up @ sun_dir) == pytest.approx(math.sin(math.radians(sun_alt)))


def test_nadir_ray_reaches_antipode() -> None:
    origin, east, north, up = _observer_basis(35.0, 139.0, 0.0)
    nadir = altaz_to_earth_fixed(np.array([-90.0]), np.array([0.0]), east=east, north=north, up=up)

    points, valid = far_surface_points(nadir, origin=origin)

    assert bool(valid[0])
    np.testing.assert_allclose(points[0], -up, atol=1.0e-9)


def test_ray_far_point_uses_half_angle_depression() -> None:
    origin, east, north, up = _observer_basis(0.0, 0.0, 0.0)
    dirs = altaz_to_earth_fixed(np.array([-20.0]), np.array([90.0]), east=east, north=north, up=up)

    points, _valid = far_surface_points(dirs, origin=origin)

    central_angle = math.degrees(math.acos(float(points[0] @ up)))
    assert central_angle == pytest.approx(40.0, abs=1.0e-6)


def test_elevated_observer_ray_just_below_horizon_misses_globe() -> None:
    origin, east, north, up = _observer_basis(0.0, 0.0, 10_000.0)
    dirs = altaz_to_earth_fixed(np.array([-0.5]), np.array([0.0]), east=east, north=north, up=up)

    _points, valid = far_surface_points(dirs, origin=origin)

    assert not bool(valid[0])


def test_visible_altitude_limit_uses_horizon_margin_and_terrain() -> None:
    az = np.array([0.0, 90.0], dtype=np.float32)

    flat = visible_altitude_limit_deg(az, 0.0, None)
    np.testing.assert_allclose(flat, [-1.0, -1.0])

    high = visible_altitude_limit_deg(az, 1000.0, None)
    assert float(high[0]) < -1.5

    terrain = visible_altitude_limit_deg(az, 0.0, [(-5.0, 0.0), (-5.0, 180.0)])
    np.testing.assert_allclose(terrain, [-5.0, -5.0], atol=1.0e-5)


def test_peak_alpha_scales_with_opacity_and_boost() -> None:
    assert earth_guide_peak_alpha(0.0) == 0.0
    base = earth_guide_peak_alpha(0.028)
    assert 0.0 < base < 1.0
    assert earth_guide_peak_alpha(0.1) > base
    assert earth_guide_peak_alpha(0.028, visibility_boost=2.0) == pytest.approx(
        min(1.0, base * 2.0)
    )


def _render(sun_altaz: tuple[float, float] | None, **kwargs) -> np.ndarray | None:
    return render_earth_guide_rgba(
        120,
        120,
        geometry=ScreenGeometry(center=(60, 60), radius=55),
        viewer_data=_viewer(**kwargs),
        sun_altaz=sun_altaz,
        earth_guide_opacity=0.028,
    )


def test_render_draws_only_below_horizon() -> None:
    rgba = _render((-40.0, 0.0), view_center=(0.0, 180.0))

    assert rgba is not None
    alpha = rgba[..., 3]
    assert int(np.count_nonzero(alpha[:60, :])) == 0
    assert int(np.count_nonzero(alpha[60:, :])) > 100


def test_render_far_side_is_dark_when_sun_is_overhead() -> None:
    # With the Sun at the observer's zenith, the region around the antipode
    # (seen near the nadir) is in night.
    rgba = _render(
        (90.0, 0.0), view_center=(-90.0, 0.0), edge_fov_deg=20.0, content_fov_deg=20.0
    )

    assert rgba is None or int(rgba[..., 3].max()) == 0


def test_render_far_side_is_lit_at_local_midnight() -> None:
    rgba = _render((-90.0, 0.0), view_center=(-90.0, 0.0))

    assert rgba is not None
    alpha = rgba[..., 3]
    expected_peak = earth_guide_peak_alpha(0.028) * 255.0
    assert int(alpha.max()) == pytest.approx(expected_peak, abs=2.0)


def test_render_uses_land_and_ocean_colors() -> None:
    # Tokyo's antipode lies in the South Atlantic near South America.
    rgba = _render((-90.0, 0.0), view_center=(-90.0, 0.0))

    assert rgba is not None
    visible = rgba[rgba[..., 3] > 0][:, :3].astype(int)
    land = np.asarray(EARTH_GUIDE_LINE_COLOR)
    ocean = np.asarray(EARTH_GUIDE_OCEAN_COLOR)
    near_land = np.all(np.abs(visible - land) <= 2, axis=1)
    near_ocean = np.all(np.abs(visible - ocean) <= 2, axis=1)
    assert int(np.count_nonzero(near_land)) > 20
    assert int(np.count_nonzero(near_ocean)) > 20


def test_render_without_sun_uses_uniform_level() -> None:
    rgba = _render(None, view_center=(-90.0, 0.0))

    assert rgba is not None
    expected = earth_guide_peak_alpha(0.028) * EARTH_GUIDE_DAY_MIN_LEVEL * 255.0
    assert int(rgba[..., 3].max()) == pytest.approx(expected, abs=2.0)


def test_render_skips_zero_opacity() -> None:
    rgba = render_earth_guide_rgba(
        64,
        64,
        geometry=ScreenGeometry(center=(32, 32), radius=30),
        viewer_data=_viewer(),
        sun_altaz=(-90.0, 0.0),
        earth_guide_opacity=0.0,
    )

    assert rgba is None


@pytest.mark.parametrize("fast_mode", [False, True])
def test_draw_earth_guide_paints_below_horizon(fast_mode: bool) -> None:
    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        draw_earth_guide(
            painter,
            geometry=ScreenGeometry(center=(120, 120), radius=100),
            viewer_data=_viewer(view_center=(0.0, 180.0)),
            sun_altaz=(-40.0, 0.0),
            earth_guide_opacity=0.028,
            fast_mode=fast_mode,
        )
    finally:
        painter.end()

    alpha = qimage_to_np_rgba(image)[..., 3]
    assert int(np.count_nonzero(alpha[:110, :])) == 0
    assert int(np.count_nonzero(alpha[130:, :])) > 1000


def test_overlay_earth_guide_forwards_sun_and_fast_mode(monkeypatch) -> None:
    calls: list[dict[str, object]] = []

    monkeypatch.setattr(
        render_composite,
        "draw_earth_guide",
        lambda *_args, **kwargs: calls.append(kwargs),
    )

    image = QImage(16, 16, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    out = render_composite._overlay_earth_guide(
        image,
        geometry=ScreenGeometry(center=(8, 8), radius=6),
        viewer_data=_viewer(view_center=(0.0, 180.0)),
        terrain_profile_altaz=None,
        earth_guide_opacity=0.028,
        fast_mode=True,
        sun_altaz=(-12.0, 270.0),
    )

    assert len(calls) == 1
    assert calls[0]["fast_mode"] is True
    assert calls[0]["sun_altaz"] == (-12.0, 270.0)
    assert calls[0]["target_size"] == (16, 16)
    assert out.size() == image.size()


def test_overlay_earth_guide_skips_zero_opacity(monkeypatch) -> None:
    calls: list[object] = []
    monkeypatch.setattr(
        render_composite,
        "draw_earth_guide",
        lambda *_args, **kwargs: calls.append(kwargs),
    )
    image = QImage(16, 16, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)

    render_composite._overlay_earth_guide(
        image,
        geometry=ScreenGeometry(center=(8, 8), radius=6),
        viewer_data=_viewer(),
        earth_guide_opacity=0.0,
    )

    assert calls == []


def _sin(deg: float) -> np.ndarray:
    return np.array([math.sin(math.radians(deg))])


def test_night_outline_weight_starts_where_twilight_fill_ends() -> None:
    assert float(earth_guide_night_outline_weight(_sin(10.0))[0]) == 0.0
    assert float(
        earth_guide_night_outline_weight(_sin(EARTH_GUIDE_TWILIGHT_LIMIT_DEG))[0]
    ) == pytest.approx(0.0, abs=1.0e-12)
    assert float(earth_guide_sunlight_level(_sin(EARTH_GUIDE_TWILIGHT_LIMIT_DEG))[0]) == pytest.approx(
        EARTH_GUIDE_NIGHT_LEVEL
    )
    assert float(
        earth_guide_night_outline_weight(_sin(EARTH_GUIDE_NIGHT_OUTLINE_FULL_DEG))[0]
    ) == pytest.approx(1.0)
    assert float(earth_guide_night_outline_weight(_sin(-80.0))[0]) == pytest.approx(1.0)


def test_seam_segments_are_not_coastline() -> None:
    assert _is_seam_segment(np.array([180.0, -84.7]), np.array([180.0, -90.0]))
    assert _is_seam_segment(np.array([180.0, -90.0]), np.array([-180.0, -90.0]))
    assert _is_seam_segment(np.array([180.0, 65.0]), np.array([180.0, 70.0]))
    assert not _is_seam_segment(np.array([179.0, 65.0]), np.array([180.0, 66.0]))


def test_coastline_segments_are_short_and_skip_seams() -> None:
    seg_a, seg_b = load_coastline_segments(0)

    assert len(seg_a) == len(seg_b) > 100
    angles = np.degrees(np.arccos(np.clip(np.sum(seg_a * seg_b, axis=1), -1.0, 1.0)))
    # Pieces are split evenly in chord space, so arcs may slightly exceed 0.5 deg.
    assert float(angles.max()) <= 0.55
    # No segment runs along the antimeridian.
    lon_a = np.degrees(np.arctan2(seg_a[:, 1], seg_a[:, 0]))
    lon_b = np.degrees(np.arctan2(seg_b[:, 1], seg_b[:, 0]))
    on_seam = (np.abs(np.abs(lon_a) - 180.0) < 1.0e-6) & (
        np.abs(np.abs(lon_b) - 180.0) < 1.0e-6
    )
    assert not np.any(on_seam)


def test_near_fade_hides_land_close_to_observer() -> None:
    _origin, _east, _north, up = _observer_basis(0.0, 0.0, 0.0)
    km_per_deg = 6371.0088 * math.pi / 180.0

    def fade_at(km: float) -> float:
        angle = math.radians(km / km_per_deg)
        point = np.array([[math.cos(angle), math.sin(angle), 0.0]])
        return float(earth_guide_near_fade(point, up)[0])

    assert fade_at(EARTH_GUIDE_NEAR_CUTOFF_KM - 50.0) == 0.0
    assert 0.0 < fade_at(EARTH_GUIDE_NEAR_CUTOFF_KM + EARTH_GUIDE_NEAR_FADE_KM / 2) < 1.0
    assert fade_at(EARTH_GUIDE_NEAR_CUTOFF_KM + EARTH_GUIDE_NEAR_FADE_KM + 50.0) == pytest.approx(1.0)


def test_night_outlines_drawn_only_on_night_side() -> None:
    geometry = ScreenGeometry(center=(200, 200), radius=190)
    # A narrow view around the nadir shows only the region near the antipode.
    viewer = _viewer(view_center=(-90.0, 0.0), edge_fov_deg=40.0, content_fov_deg=40.0)

    # Sun at the zenith: the far side is in night, so outlines appear.
    _a, _b, night_weight = compute_night_outline_segments(
        geometry=geometry, viewer_data=viewer, sun_altaz=(90.0, 0.0)
    )
    # Sun at the nadir: the far side is in daylight, so there are none.
    _a, _b, day_weight = compute_night_outline_segments(
        geometry=geometry, viewer_data=viewer, sun_altaz=(-90.0, 0.0)
    )
    _a, _b, no_sun_weight = compute_night_outline_segments(
        geometry=geometry, viewer_data=viewer, sun_altaz=None
    )

    assert len(night_weight) > 100
    assert float(night_weight.max()) <= 1.0
    assert len(day_weight) == 0
    assert len(no_sun_weight) == 0


def test_draw_earth_guide_draws_white_night_outlines() -> None:
    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        draw_earth_guide(
            painter,
            geometry=ScreenGeometry(center=(120, 120), radius=110),
            viewer_data=_viewer(
                view_center=(-90.0, 0.0), edge_fov_deg=40.0, content_fov_deg=40.0
            ),
            sun_altaz=(90.0, 0.0),
            earth_guide_opacity=0.028,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    drawn = arr[arr[..., 3] > 0]
    assert len(drawn) > 50
    # Outline pixels are neutral (white), unlike the land/ocean fill colors.
    # Skip faint antialiased pixels whose color is dominated by rounding.
    rgb = drawn[drawn[:, 3] >= 32][:, :3].astype(int)
    assert len(rgb) > 20
    assert np.all(np.abs(rgb[:, 0] - rgb[:, 2]) <= 2)
