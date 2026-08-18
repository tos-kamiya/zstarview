from __future__ import annotations

from types import SimpleNamespace

import astropy.time
import numpy as np
import pytest
from PySide6.QtCore import QRect
from PySide6.QtGui import QColor, QImage, QPainter

from zstarview.gui.window_render import SkyWindowRenderMixin
from zstarview.render.pipeline import _draw_mesh_transformed_star_surface
from zstarview.render.qt_image import qimage_to_np_rgba
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


def test_expanded_mesh_maps_viewport_points_through_guard_origin() -> None:
    source = np.asarray(
        ((0.0, 0.0), (30.0, 0.0), (0.0, 30.0), (30.0, 30.0))
    )
    mesh = StarInterpolationMesh(
        source_vertices=source,
        target_vertices=source + np.asarray((2.0, 3.0)),
        columns=1,
        rows=1,
        viewport_origin=(10.0, 10.0),
    )

    assert np.allclose(
        mesh.map_viewport_points(np.asarray(((4.0, 6.0),))),
        ((6.0, 9.0),),
    )


def test_star_interpolation_mesh_scales_to_render_surface_coordinates() -> None:
    mesh = build_star_interpolation_mesh(
        width_px=3840,
        height_px=2160,
        geometry_center=(1920.0, 1080.0),
        geometry_radius=1080.0,
        view_center_altaz_deg=(45.0, 0.0),
        observer_lat_deg=35.0,
        edge_fov_deg=90.0,
        elapsed_seconds=30.0,
    )

    scaled = mesh.scaled(2024.0 / 3840.0, 1138.0 / 2160.0)

    assert (scaled.columns, scaled.rows) == (39, 22)
    assert np.allclose(scaled.source_vertices[-1], (2024.0, 1138.0))
    assert np.allclose(
        scaled.target_vertices,
        mesh.target_vertices * np.asarray((2024.0 / 3840.0, 1138.0 / 2160.0)),
    )


def test_low_resolution_mesh_surface_applies_affine_translation() -> None:
    source_vertices = np.asarray(
        ((0.0, 0.0), (20.0, 0.0), (0.0, 20.0), (20.0, 20.0))
    )
    mesh = StarInterpolationMesh(
        source_vertices=source_vertices,
        target_vertices=source_vertices + np.asarray((2.0, 4.0)),
        columns=1,
        rows=1,
    )
    source = QImage(10, 10, QImage.Format.Format_ARGB32_Premultiplied)
    source.fill(0)
    source.setPixelColor(4, 3, QColor(255, 255, 255, 255))
    destination = QImage(20, 20, QImage.Format.Format_ARGB32_Premultiplied)
    destination.fill(0)
    painter = QPainter(destination)
    try:
        _draw_mesh_transformed_star_surface(
            painter,
            source,
            mesh=mesh,
            viewport_rect=QRect(0, 0, 20, 20),
        )
    finally:
        painter.end()

    alpha = qimage_to_np_rgba(destination)[..., 3]
    assert int(alpha[6, 8]) == 0
    assert int(alpha[10, 10]) > 0


def test_mesh_surface_applies_affine_rotation_about_nonzero_center() -> None:
    source_vertices = np.asarray(
        ((0.0, 0.0), (100.0, 0.0), (0.0, 100.0), (100.0, 100.0))
    )
    target_vertices = np.column_stack(
        (100.0 - source_vertices[:, 1], source_vertices[:, 0])
    )
    mesh = StarInterpolationMesh(
        source_vertices=source_vertices,
        target_vertices=target_vertices,
        columns=1,
        rows=1,
    )
    source = QImage(100, 100, QImage.Format.Format_ARGB32_Premultiplied)
    source.fill(0)
    source.setPixelColor(70, 40, QColor(255, 255, 255, 255))
    destination = QImage(100, 100, QImage.Format.Format_ARGB32_Premultiplied)
    destination.fill(0)
    painter = QPainter(destination)
    try:
        _draw_mesh_transformed_star_surface(
            painter,
            source,
            mesh=mesh,
            viewport_rect=QRect(0, 0, 100, 100),
        )
    finally:
        painter.end()

    alpha = qimage_to_np_rgba(destination)[..., 3]
    assert int(alpha[40, 70]) == 0
    assert np.any(alpha[69:72, 58:62] > 0)


def test_cached_star_surface_keeps_four_k_internal_render_size(monkeypatch) -> None:
    captured_size: list[tuple[int, int]] = []
    sentinel = QImage(1, 1, QImage.Format.Format_ARGB32_Premultiplied)

    def fake_render_cached_image(_self, *, image_size, **_kwargs):
        captured_size.append((int(image_size.width()), int(image_size.height())))
        return sentinel

    monkeypatch.setattr(
        SkyWindowRenderMixin,
        "_render_cached_image",
        fake_render_cached_image,
    )
    frame = SimpleNamespace(
        viewport_rect=QRect(0, 0, 3840, 2160),
        geometry=SimpleNamespace(center=(1920, 1080), radius=1080),
    )
    render_inputs = SimpleNamespace(
        style=SimpleNamespace(
            star_render_expected_width=600,
            vmag_limit=6.0,
        ),
    )

    result = SkyWindowRenderMixin._render_cached_star_surface_image(
        SimpleNamespace(),
        base_frame_key=("frame",),
        frame=frame,
        render_inputs=render_inputs,
        faint_only=True,
    )

    assert result is sentinel
    assert captured_size == [(2058, 1172)]


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
