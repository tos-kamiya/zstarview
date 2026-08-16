from __future__ import annotations

from types import SimpleNamespace

import astropy.time
import numpy as np
import pytest

from zstarview.gui.window_render import SkyWindowRenderMixin
from zstarview.render.star_interpolation import (
    STAR_INTERPOLATION_COVERAGE,
    STAR_INTERPOLATION_MAX_UPDATE_INTERVAL_SECONDS,
    StarInterpolationMesh,
    build_star_interpolation_mesh,
    should_interpolate_stars,
)


def test_star_interpolation_coverage_is_configurable() -> None:
    assert 0.0 < STAR_INTERPOLATION_COVERAGE <= 1.0


def test_star_interpolation_is_disabled_above_ninety_second_update_interval() -> None:
    assert STAR_INTERPOLATION_MAX_UPDATE_INTERVAL_SECONDS == 90
    assert should_interpolate_stars(90)
    assert not should_interpolate_stars(90.1)


def test_star_interpolation_mesh_uses_square_100px_cells() -> None:
    mesh = build_star_interpolation_mesh(
        width_px=1600,
        height_px=900,
        geometry_center=(800.0, 450.0),
        geometry_radius=450.0,
        view_center_altaz_deg=(45.0, 0.0),
        observer_lat_deg=35.0,
        edge_fov_deg=90.0,
        elapsed_seconds=30.0,
    )

    assert (mesh.columns, mesh.rows) == (16, 9)
    assert mesh.source_vertices.shape == (17 * 10, 2)
    assert mesh.target_vertices.shape == mesh.source_vertices.shape
    assert np.all(np.isfinite(mesh.target_vertices))
    assert np.allclose(mesh.source_vertices[0], (0.0, 0.0))
    assert np.allclose(mesh.source_vertices[-1], (1600.0, 900.0))


def test_star_mesh_cache_key_tracks_projection_inputs() -> None:
    snapshot_time = astropy.time.Time("2026-08-17T00:00:00", scale="utc")

    def make_frame(
        *,
        center: tuple[float, float] = (800.0, 450.0),
        location: tuple[float, float] = (35.0, 139.0),
    ) -> SimpleNamespace:
        return SimpleNamespace(
            viewport_rect=SimpleNamespace(width=lambda: 1600, height=lambda: 900),
            geometry=SimpleNamespace(center=center, radius=450.0),
            viewer=SimpleNamespace(
                view_center=(45.0, 180.0),
                location=location,
                edge_fov_deg=90.0,
            ),
            sky_update_interval=60,
            time_obj=snapshot_time,
        )

    scene = SimpleNamespace(
        celestial_data=SimpleNamespace(star_time=snapshot_time, time=snapshot_time)
    )
    window = SimpleNamespace()
    baseline = SkyWindowRenderMixin._star_mesh_cache_key(
        window, make_frame(), scene
    )

    assert baseline != SkyWindowRenderMixin._star_mesh_cache_key(
        window, make_frame(center=(801.0, 450.0)), scene
    )
    assert baseline != SkyWindowRenderMixin._star_mesh_cache_key(
        window, make_frame(location=(-35.0, 139.0)), scene
    )


def test_star_interpolation_mesh_maps_points_without_separate_dimensions() -> None:
    source = np.asarray(
        ((0.0, 0.0), (10.0, 0.0), (0.0, 10.0), (10.0, 10.0))
    )
    mesh = StarInterpolationMesh(
        source_vertices=source,
        target_vertices=source + np.asarray((2.0, 3.0)),
        columns=1,
        rows=1,
    )

    assert np.allclose(mesh.map_points(np.asarray(((4.0, 6.0),))), ((6.0, 9.0),))


def test_star_interpolation_mesh_rejects_mismatched_dimensions() -> None:
    with pytest.raises(
        ValueError, match="vertices must match its dimensions"
    ):
        StarInterpolationMesh(
            source_vertices=np.zeros((4, 2)),
            target_vertices=np.zeros((4, 2)),
            columns=2,
            rows=1,
        )
