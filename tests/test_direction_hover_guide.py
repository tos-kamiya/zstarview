from __future__ import annotations

from PySide6.QtCore import QPoint, Qt
from PySide6.QtGui import QImage, QPainter

from zstarview.astro import altaz_to_normalized_xy
from zstarview.render.geometry import normalized_to_screen_xy
from zstarview.render.guides import (
    draw_direction_grid_overlay,
    resolve_direction_marker_hover,
)
from zstarview.types import ScreenGeometry, ViewerData


def test_resolve_direction_marker_hover_targets_marker_not_label() -> None:
    geometry = ScreenGeometry(center=(200, 200), radius=180)
    viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(0.0, 0.0),
        edge_fov_deg=95.0,
        content_fov_deg=100.0,
    )

    marker_nx, marker_ny = altaz_to_normalized_xy(
        0.0,
        0.0,
        viewer_data.view_center,
        edge_fov_deg=95.0,
    )
    marker_x, marker_y = normalized_to_screen_xy(marker_nx, marker_ny, geometry)
    hovered = resolve_direction_marker_hover(
        geometry,
        viewer_data,
        QPoint(int(round(marker_x)), int(round(marker_y))),
    )
    assert hovered is not None
    assert hovered.label == "N"

    hovered_label = resolve_direction_marker_hover(
        geometry,
        viewer_data,
        QPoint(int(round(marker_x + 40.0)), int(round(marker_y))),
    )
    assert hovered_label is None


def test_draw_direction_grid_overlay_draws_major_guides_and_minor_crosses(monkeypatch) -> None:
    geometry = ScreenGeometry(center=(200, 200), radius=180)
    viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(0.0, 0.0),
        edge_fov_deg=95.0,
        content_fov_deg=100.0,
    )

    call_count = 0
    widths: list[float] = []
    cross_count = 0
    cross_widths: list[float] = []
    cross_half_lens: list[float] = []

    def fake_draw_direction_polyline(*_args, **_kwargs) -> None:
        nonlocal call_count
        call_count += 1
        widths.append(float(_kwargs["width"]))

    def fake_draw_direction_cross_marker(*_args, **_kwargs) -> None:
        nonlocal cross_count
        cross_count += 1
        cross_widths.append(float(_kwargs["width"]))
        cross_half_lens.append(float(_kwargs["half_len"]))

    monkeypatch.setattr(
        "zstarview.render.guides._draw_direction_polyline",
        fake_draw_direction_polyline,
    )
    monkeypatch.setattr(
        "zstarview.render.guides._draw_direction_cross_marker",
        fake_draw_direction_cross_marker,
    )

    image = QImage(400, 400, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(Qt.GlobalColor.transparent)
    painter = QPainter(image)
    try:
        draw_direction_grid_overlay(
            painter,
            geometry,
            viewer_data,
            (400, 400),
        )
    finally:
        painter.end()

    assert call_count == 17
    assert len(widths) == 17
    assert all(abs(width - 0.51) < 1.0e-9 for width in widths)
    assert cross_count > 0
    assert len(cross_widths) == cross_count
    assert all(abs(width - 0.51) < 1.0e-9 for width in cross_widths)
    assert all(abs(half_len - 2.4) < 1.0e-9 for half_len in cross_half_lens)
