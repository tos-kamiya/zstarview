from __future__ import annotations

from types import SimpleNamespace

import astropy.time
import numpy as np

from zstarview.render.star_interpolation import (
    STAR_INTERPOLATION_COVERAGE,
    STAR_INTERPOLATION_MAX_UPDATE_INTERVAL_SECONDS,
    apply_homography,
    build_star_interpolation_homography,
    build_star_interpolation_mesh,
    should_interpolate_stars,
)
from zstarview.render.star_interpolation import (
    _direction_to_screen,
    _rotate_about_axis,
    _screen_to_direction,
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


def test_star_interpolation_mesh_uses_square_100px_cells() -> None:
    source, target = build_star_interpolation_mesh(
        width_px=1600,
        height_px=900,
        geometry_center=(800.0, 450.0),
        geometry_radius=450.0,
        view_center_altaz_deg=(45.0, 0.0),
        observer_lat_deg=35.0,
        edge_fov_deg=90.0,
        elapsed_seconds=30.0,
    )

    assert source.shape == (17 * 10, 2)
    assert target.shape == source.shape
    assert np.all(np.isfinite(target))
    assert np.allclose(source[0], (0.0, 0.0))
    assert np.allclose(source[-1], (1600.0, 900.0))


def test_star_interpolation_gate_currently_blocks_homography_go(capsys) -> None:
    """Record that the current homography does not pass the 1px gate.

    The reference is the same 3D sidereal rotation and screen projection used
    to construct the homography, evaluated on a denser grid.  The 75% coverage
    factor is intentionally included in the reference elapsed time because it
    is part of the current display policy.
    """
    width, height = 1600.0, 900.0
    center = (800.0, 450.0)
    radius = 450.0
    edge_fov_deg = 90.0
    observer_lat_deg = 35.0
    grid_x, grid_y = np.meshgrid(
        np.linspace(0.04 * width, 0.96 * width, 25),
        np.linspace(0.04 * height, 0.96 * height, 15),
    )
    points = np.column_stack((grid_x.ravel(), grid_y.ravel()))
    normalized_x = (points[:, 0] - center[0]) / radius
    normalized_y = (points[:, 1] - center[1]) / radius
    points = points[np.hypot(normalized_x, normalized_y) <= 1.0]
    cases = (
        (45.0, 0.0, "north"),
        (45.0, 180.0, "south"),
        (80.0, 90.0, "zenith"),
        (10.0, 270.0, "low-altitude"),
    )
    errors: list[float] = []
    report: list[str] = []
    for view_alt, view_az, name in cases:
        for elapsed_seconds in (-30.0, 30.0):
            matrix = build_star_interpolation_homography(
                width_px=int(width),
                height_px=int(height),
                geometry_center=center,
                geometry_radius=radius,
                view_center_altaz_deg=(view_alt, view_az),
                observer_lat_deg=observer_lat_deg,
                edge_fov_deg=edge_fov_deg,
                elapsed_seconds=elapsed_seconds,
            )
            approximate = apply_homography(points, matrix)
            directions = _screen_to_direction(
                points[:, 0],
                points[:, 1],
                width_px=width,
                height_px=height,
                geometry_center=center,
                geometry_radius=radius,
                view_center_altaz_deg=(view_alt, view_az),
                edge_fov_deg=edge_fov_deg,
            )
            latitude = np.radians(observer_lat_deg)
            pole_axis = np.array(
                [np.cos(latitude), 0.0, np.sin(latitude)], dtype=float
            )
            angle = np.radians(
                (360.0 / 86164.0905)
                * elapsed_seconds
                * STAR_INTERPOLATION_COVERAGE
            )
            exact = _direction_to_screen(
                _rotate_about_axis(directions, pole_axis, angle),
                width_px=width,
                height_px=height,
                geometry_center=center,
                geometry_radius=radius,
                view_center_altaz_deg=(view_alt, view_az),
                edge_fov_deg=edge_fov_deg,
            )
            point_errors = np.linalg.norm(approximate - exact, axis=1)
            maximum = float(np.max(point_errors))
            mean = float(np.mean(point_errors))
            errors.append(maximum)
            report.append(f"{name} {elapsed_seconds:+.0f}s mean={mean:.4f}px max={maximum:.4f}px")

    output = "star interpolation gate: " + "; ".join(report)
    print(output)
    # This is deliberately a failing-go gate: mesh implementation may proceed
    # only after the current global homography is shown to exceed this error
    # budget in at least one representative view.
    assert max(errors) >= 1.0, capsys.readouterr().out
