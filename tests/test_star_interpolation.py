from __future__ import annotations

from types import SimpleNamespace

import astropy.time
import numpy as np

from zstarview.render.star_interpolation import (
    STAR_INTERPOLATION_COVERAGE,
    STAR_INTERPOLATION_MAX_UPDATE_INTERVAL_SECONDS,
    apply_homography,
    build_star_interpolation_homography,
    should_interpolate_stars,
)
from zstarview.render.zstarview_pipeline import _star_interpolation_matrix


def test_star_interpolation_coverage_leaves_a_small_boundary_step() -> None:
    assert 0.0 < STAR_INTERPOLATION_COVERAGE < 1.0


def test_star_interpolation_is_disabled_above_ninety_second_update_interval() -> None:
    assert STAR_INTERPOLATION_MAX_UPDATE_INTERVAL_SECONDS == 90
    assert should_interpolate_stars(90)
    assert not should_interpolate_stars(90.1)


def test_pipeline_disables_star_interpolation_above_ninety_seconds() -> None:
    snapshot_time = astropy.time.Time("2026-07-19T00:00:00", scale="utc")
    frame = SimpleNamespace(
        sky_update_interval=91,
        time_obj=snapshot_time,
        geometry=SimpleNamespace(center=(800.0, 450.0), radius=450.0),
        viewport_rect=SimpleNamespace(width=lambda: 1600, height=lambda: 900),
        viewer=SimpleNamespace(
            view_center=(45.0, 180.0), location=(35.0, 139.0), edge_fov_deg=90.0
        ),
    )
    scene = SimpleNamespace(
        celestial_data=SimpleNamespace(star_time=snapshot_time, time=snapshot_time)
    )

    assert _star_interpolation_matrix(frame=frame, scene=scene) is None


def test_pipeline_uses_the_configured_interval_as_interpolation_limit() -> None:
    snapshot_time = astropy.time.Time("2026-07-19T00:00:00", scale="utc")
    frame = SimpleNamespace(
        sky_update_interval=90,
        time_obj=snapshot_time + astropy.time.TimeDelta(45.0, format="sec"),
        geometry=SimpleNamespace(center=(800.0, 450.0), radius=450.0),
        viewport_rect=SimpleNamespace(width=lambda: 1600, height=lambda: 900),
        viewer=SimpleNamespace(
            view_center=(45.0, 180.0), location=(35.0, 139.0), edge_fov_deg=90.0
        ),
    )
    scene = SimpleNamespace(
        celestial_data=SimpleNamespace(star_time=snapshot_time, time=snapshot_time)
    )

    matrix = _star_interpolation_matrix(frame=frame, scene=scene)

    assert matrix is not None


def test_star_interpolation_is_identity_at_snapshot_time() -> None:
    matrix = build_star_interpolation_homography(
        width_px=1600,
        height_px=900,
        geometry_center=(800.0, 450.0),
        geometry_radius=450.0,
        view_center_altaz_deg=(45.0, 180.0),
        observer_lat_deg=35.0,
        edge_fov_deg=90.0,
        elapsed_seconds=0.0,
    )

    assert np.allclose(matrix, np.eye(3), atol=1.0e-12)


def test_star_interpolation_moves_a_centered_star_with_sidereal_time() -> None:
    matrix = build_star_interpolation_homography(
        width_px=1600,
        height_px=900,
        geometry_center=(800.0, 450.0),
        geometry_radius=450.0,
        view_center_altaz_deg=(45.0, 180.0),
        observer_lat_deg=35.0,
        edge_fov_deg=90.0,
        elapsed_seconds=30.0,
    )
    center = np.array([[800.0, 450.0]])

    mapped = apply_homography(center, matrix)

    assert np.all(np.isfinite(mapped))
    assert not np.allclose(mapped, center, atol=1.0e-6)


def test_star_interpolation_supports_the_first_half_of_a_centered_interval() -> None:
    matrix = build_star_interpolation_homography(
        width_px=1600,
        height_px=900,
        geometry_center=(800.0, 450.0),
        geometry_radius=450.0,
        view_center_altaz_deg=(45.0, 180.0),
        observer_lat_deg=35.0,
        edge_fov_deg=90.0,
        elapsed_seconds=-30.0,
    )
    center = np.array([[800.0, 450.0]])

    mapped = apply_homography(center, matrix)

    assert np.all(np.isfinite(mapped))
    assert not np.allclose(mapped, center, atol=1.0e-6)


def test_star_interpolation_has_finite_large_view_transform() -> None:
    matrix = build_star_interpolation_homography(
        width_px=2560,
        height_px=1440,
        geometry_center=(1280.0, 720.0),
        geometry_radius=720.0,
        view_center_altaz_deg=(90.0, 180.0),
        observer_lat_deg=35.0,
        edge_fov_deg=90.0,
        elapsed_seconds=60.0,
    )
    points = np.array(
        [
            [100.0, 100.0],
            [1280.0, 720.0],
            [2460.0, 100.0],
            [100.0, 1340.0],
            [2460.0, 1340.0],
        ]
    )

    mapped = apply_homography(points, matrix)

    assert np.all(np.isfinite(matrix))
    assert np.all(np.isfinite(mapped))
