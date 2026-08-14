
from tests._window_render_support import *


def test_draw_urban_outlines_clips_two_point_outline_out_of_view() -> None:
    class _Painter:
        def __init__(self) -> None:
            self.polylines: list[list[tuple[float, float]]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            if not hasattr(pen, "color") or not hasattr(pen, "widthF"):
                return

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, poly) -> None:
            self.polylines.append(
                [(poly.at(i).x(), poly.at(i).y()) for i in range(poly.count())]
            )

    painter = _Painter()
    render_terrain_module.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[
            UrbanOutlinePolyline(points=[(-10.0, 10.0), (-12.0, 10.3)], height_m=50.0)
        ],
        viewer=_viewer(),
        opacity=0.38,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert painter.polylines == []


def test_draw_urban_outlines_uses_fixed_alpha_and_near_underlay(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )

    class _Painter:
        def __init__(self) -> None:
            self.alpha_values: list[int] = []
            self.width_values: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            if not hasattr(pen, "color") or not hasattr(pen, "widthF"):
                return
            self.alpha_values.append(int(pen.color().alpha()))
            self.width_values.append(float(pen.widthF()))

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, _poly) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[
            UrbanOutlinePolyline(
                points=[(-10.0, 10.0), (-12.0, 12.0)],
                height_m=0.0,
                distance_km=0.01,
            ),
            UrbanOutlinePolyline(
                points=[(-10.0, 20.0), (-12.0, 22.0)],
                height_m=50.0,
                distance_km=0.5,
            ),
        ],
        viewer=_viewer(),
        opacity=0.2,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert painter.alpha_values == [3, 2, 1, 51, 3, 2, 1, 51]
    assert [round(width, 1) for width in painter.width_values] == [
        9.2,
        7.2,
        4.4,
        2.3,
        2.4,
        2.0,
        2.0,
        1.3,
    ]


def test_draw_urban_outlines_can_remove_fill_without_removing_outline(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )
    monkeypatch.setattr(render_terrain_module, "is_in_fov", lambda *_args, **_kwargs: True)

    class _Painter:
        def __init__(self) -> None:
            self.brush_alphas: list[int] = []
            self.pen_alphas: list[int] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            if hasattr(pen, "color") and hasattr(pen, "widthF"):
                self.pen_alphas.append(int(pen.color().alpha()))

        def setBrush(self, brush) -> None:
            if hasattr(brush, "alpha"):
                self.brush_alphas.append(int(brush.alpha()))

        def drawPolyline(self, _poly) -> None:
            pass

        def drawPolygon(self, _polygon) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=100),
        urban_outlines=[
            UrbanOutlinePolyline(
                points=[(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0), (0.0, 0.0)],
                height_m=0.0,
                distance_km=1.0,
            )
        ],
        viewer=_viewer(),
        opacity=0.2,
        fill_opacity_factor=0.0,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert painter.brush_alphas == []
    assert painter.pen_alphas


def test_urban_outline_fill_alpha_applies_brightness_scale() -> None:
    assert render_terrain_module._urban_outline_fill_alpha(0.2) == 23


def test_draw_urban_outlines_allows_sub_unit_width_scale(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )

    class _Painter:
        def __init__(self) -> None:
            self.width_values: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            if not hasattr(pen, "widthF"):
                return
            self.width_values.append(float(pen.widthF()))

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, _poly) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[
            UrbanOutlinePolyline(
                points=[(-10.0, 10.0), (-12.0, 12.0)],
                height_m=50.0,
                distance_km=0.01,
            )
        ],
        viewer=_viewer(),
        opacity=0.38,
        line_width_scale=0.5,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(alt),
            float(az),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert [round(width, 2) for width in painter.width_values[:4]] == [
        4.6,
        3.6,
        2.2,
        1.14,
    ]


def test_draw_urban_outlines_scales_alpha_for_subpixel_widths(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )

    class _Painter:
        def __init__(self) -> None:
            self.alpha_values: list[int] = []
            self.width_values: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            if not hasattr(pen, "color") or not hasattr(pen, "widthF"):
                return
            self.alpha_values.append(int(pen.color().alpha()))
            self.width_values.append(float(pen.widthF()))

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, _poly) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[
            UrbanOutlinePolyline(
                points=[(-10.0, 10.0), (-12.0, 12.0)],
                height_m=0.0,
                distance_km=0.01,
            )
        ],
        viewer=_viewer(),
        opacity=0.2,
        line_width_scale=0.1,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert painter.alpha_values == [3, 1, 0, 12]
    assert [round(width, 3) for width in painter.width_values] == [
        0.92,
        0.72,
        0.44,
        0.228,
    ]


def test_draw_urban_outlines_thickens_tall_buildings(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )

    class _Painter:
        def __init__(self) -> None:
            self.width_values: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            if not hasattr(pen, "widthF"):
                return
            self.width_values.append(float(pen.widthF()))

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, _poly) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[
            UrbanOutlinePolyline(
                points=[(-10.0, 10.0), (-12.0, 12.0)],
                height_m=600.0,
                distance_km=0.01,
            )
        ],
        viewer=_viewer(),
        opacity=0.38,
        line_width_scale=0.5,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(alt),
            float(az),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert [round(width, 2) for width in painter.width_values[:4]] == [
        9.2,
        7.2,
        4.4,
        1.14,
    ]


def test_draw_terrain_horizon_line_scales_line_widths(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    class _Painter:
        def __init__(self) -> None:
            self.pen_widths: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            self.pen_widths.append(float(pen.widthF()))

        def drawPolyline(self, _poly) -> None:
            pass

    painter = _Painter()
    render_terrain_module._draw_terrain_profile_layer(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=_viewer(),
        terrain_profile_altaz=[(0.0, 0.0), (0.1, 0.1)],
        terrain_profile_distances_m=None,
        spec=_terrain_horizon_spec(opacity=0.38, line_width_scale=2.0, fast_mode=True),
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
        split_by_gaps_func=lambda points: [points],
    )

    assert painter.pen_widths == [
        pytest.approx(render_terrain_module.TERRAIN_HORIZON_FAST_WIDTH * 2.0)
    ]


def test_draw_terrain_horizon_line_scales_widths_by_distance(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    class _Painter:
        def __init__(self) -> None:
            self.pen_widths: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            self.pen_widths.append(float(pen.widthF()))

        def drawLine(self, *_args) -> None:
            pass

        def drawPolyline(self, *_args) -> None:
            pass

    painter = _Painter()
    render_terrain_module._draw_terrain_profile_layer(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=_viewer(),
        terrain_profile_altaz=[(0.0, 0.0), (0.0, 0.1), (0.0, 0.2)],
        terrain_profile_distances_m=[1_000.0, 50_000.0, 120_000.0],
        spec=_terrain_horizon_spec(opacity=0.38, line_width_scale=1.0, fast_mode=False),
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
        split_by_gaps_func=lambda points: [points],
    )

    assert len(painter.pen_widths) == 2
    assert painter.pen_widths[0] > painter.pen_widths[1]


def test_draw_terrain_secondary_ridges_use_fixed_widths(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    class _Painter:
        def __init__(self) -> None:
            self.pen_widths: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            self.pen_widths.append(float(pen.widthF()))

        def drawLine(self, *_args) -> None:
            pass

        def drawPolyline(self, *_args) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_terrain_secondary_ridges(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=SimpleNamespace(
            view_center=(45.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        ),
        terrain_secondary_ridges_layers=[
            [(0.0, 0.0), (0.0, 0.1), (0.0, 0.2)],
            [(0.1, 0.0), (0.1, 0.1), (0.1, 0.2)],
            [(0.2, 0.0), (0.2, 0.1), (0.2, 0.2)],
        ],
        terrain_secondary_ridges_distances_m_layers=[
            [1_000.0, 2_000.0, 3_000.0],
            [10_000.0, 12_000.0, 15_000.0],
            [50_000.0, 60_000.0, 70_000.0],
        ],
        opacity=0.38,
        line_width_scale=1.0,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert len(painter.pen_widths) == 6
    assert all(width > 0.0 for width in painter.pen_widths)
    assert painter.pen_widths[0] > painter.pen_widths[1]
    assert painter.pen_widths[2] > painter.pen_widths[3]
    assert painter.pen_widths[4] > painter.pen_widths[5]
    assert painter.pen_widths[0] > painter.pen_widths[2] > painter.pen_widths[4]


def test_draw_terrain_secondary_ridges_keeps_beige_color(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    class _Painter:
        def __init__(self) -> None:
            self.pen_rgbs: list[tuple[int, int, int]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            color = QColor(pen.color())
            self.pen_rgbs.append((color.red(), color.green(), color.blue()))

        def drawLine(self, *_args) -> None:
            pass

        def drawPolyline(self, *_args) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_terrain_secondary_ridges(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=SimpleNamespace(
            view_center=(45.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        ),
        terrain_secondary_ridges_layers=[
            [(20.0, 0.0), (20.0, 0.1)],
            [(10.0, 0.0), (10.0, 0.1)],
        ],
        terrain_secondary_ridges_distances_m_layers=[
            [1_000.0, 2_000.0],
            [10_000.0, 12_000.0],
        ],
        opacity=0.38,
        line_width_scale=1.0,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert len(painter.pen_rgbs) == 4
    assert all(
        rgb == render_terrain_module.TERRAIN_HORIZON_LINE_COLOR
        for rgb in painter.pen_rgbs
    )


def test_secondary_ridge_alpha_base_is_lower() -> None:
    assert render_terrain_module.terrain_secondary_ridge_line_alpha(0.38) < 0.08


def test_secondary_ridge_overlay_alpha_keeps_hidden_runs_weaker(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    class _Painter:
        def __init__(self) -> None:
            self.alphas: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            self.alphas.append(float(pen.color().alphaF()))

        def drawLine(self, *_args) -> None:
            pass

        def drawPolyline(self, *_args) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_terrain_secondary_ridges(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=SimpleNamespace(
            view_center=(45.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        ),
        terrain_secondary_ridges_layers=[
            [(20.0, 0.0), (20.0, 0.1)],
            [(10.0, 0.0), (10.0, 0.1)],
        ],
        terrain_secondary_ridges_distances_m_layers=[
            [1_000.0, 2_000.0],
            [10_000.0, 12_000.0],
        ],
        opacity=0.38,
        line_width_scale=1.0,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert len(painter.alphas) == 4
    assert all(0.0 <= alpha <= 1.0 for alpha in painter.alphas)
    assert painter.alphas[0] > painter.alphas[2]
    assert painter.alphas[1] > painter.alphas[3]
    assert painter.alphas[2] < painter.alphas[0]


def test_draw_terrain_secondary_ridges_uses_polylines_only(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "split_by_gaps",
        lambda points: [points],
    )

    class _Painter:
        def __init__(self) -> None:
            self.lines: list[tuple[tuple[float, float], tuple[float, float]]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def drawLine(self, start, end) -> None:
            self.lines.append(
                ((float(start.x()), float(start.y())), (float(end.x()), float(end.y())))
            )

        def drawPolyline(self, *_args) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_terrain_secondary_ridges(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=SimpleNamespace(
            view_center=(45.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        ),
        terrain_secondary_ridges_layers=[
            [(5.0, 359.0), (5.0, 0.0), (5.0, 1.0)],
        ],
        terrain_secondary_ridges_distances_m_layers=[
            [1_000.0, 2_000.0, 3_000.0],
        ],
        opacity=0.38,
        line_width_scale=1.0,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            math.sin(math.radians(float(az))),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
        split_by_gaps_func=lambda points: [points],
    )

    assert painter.lines == []


def test_draw_terrain_horizon_line_uses_edge_fov_for_projection() -> None:
    class _Painter:
        def __init__(self) -> None:
            self.polylines: list[list[tuple[float, float]]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, poly) -> None:
            self.polylines.append(
                [(poly.at(i).x(), poly.at(i).y()) for i in range(poly.count())]
            )

    def _render(edge_fov_deg: float) -> list[tuple[float, float]]:
        painter = _Painter()
        render_terrain_module._draw_terrain_profile_layer(
            painter,
            geometry=SimpleNamespace(center=(0, 0), radius=1),
            viewer=_viewer(edge_fov_deg=edge_fov_deg, content_fov_deg=180.0),
            terrain_profile_altaz=[(0.0, 180.0), (0.0, 190.0)],
            terrain_profile_distances_m=None,
            spec=_terrain_horizon_spec(opacity=0.38, fast_mode=True),
            is_in_fov_func=lambda *_args, **_kwargs: True,
            altaz_to_normalized_xy_func=render_terrain_module.altaz_to_normalized_xy,
            normalized_to_screen_xy_func=lambda nx, ny, _geometry: (
                float(nx),
                float(ny),
            ),
            split_by_gaps_func=lambda points: [points],
        )
        assert painter.polylines
        return painter.polylines[0]

    points_90 = _render(90.0)
    points_120 = _render(120.0)

    assert points_120[0][1] < points_90[0][1]


def test_draw_terrain_horizon_line_rotates_profile_away_from_north_seam() -> None:
    class _Painter:
        def __init__(self) -> None:
            self.polylines: list[list[tuple[float, float]]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, poly) -> None:
            self.polylines.append(
                [(poly.at(i).x(), poly.at(i).y()) for i in range(poly.count())]
            )

    painter = _Painter()
    render_terrain_module._draw_terrain_profile_layer(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=_viewer(edge_fov_deg=120.0, content_fov_deg=180.0),
        terrain_profile_altaz=[
            (0.0, 0.0),
            (0.0, 10.0),
            (0.0, 180.0),
            (0.0, 190.0),
        ],
        terrain_profile_distances_m=None,
        spec=_terrain_horizon_spec(opacity=0.38, fast_mode=True),
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        split_by_gaps_func=lambda points: [points],
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert painter.polylines
    assert [point[0] for point in painter.polylines[0]] == [0.0, 10.0, 180.0, 190.0]


def test_draw_star_layer_forwards_outline_flag(monkeypatch) -> None:
    captured: dict[str, object] = {}

    def fake_draw_stars(*_args, **kwargs) -> None:
        captured.update(kwargs)

    monkeypatch.setattr(
        pipeline_module.render_stars, "draw_stars_fast", fake_draw_stars
    )

    pipeline_module._draw_star_layer(
        painter=object(),
        geometry=SimpleNamespace(center=(120, 120), radius=80),
        viewport_rect=QRect(0, 0, 240, 240),
        scene=_make_scene(),
        style=_make_style(bright_bodies_mode="outline"),
        star_render_surface_size=(240, 240),
        fast_mode=True,
    )

    assert captured["outline_bright_bodies"] is True
    assert "fast_mode" not in captured


def test_draw_star_layer_disables_twinkle_in_fast_mode(monkeypatch) -> None:
    captured: dict[str, object] = {}

    def fake_draw_stars(*_args, **kwargs) -> None:
        captured.update(kwargs)

    monkeypatch.setattr(
        pipeline_module.render_stars, "draw_stars_fast", fake_draw_stars
    )

    pipeline_module._draw_star_layer(
        painter=object(),
        geometry=SimpleNamespace(center=(120, 120), radius=80),
        viewport_rect=QRect(0, 0, 240, 240),
        scene=_make_scene(),
        style=_make_style(),
        star_render_surface_size=(240, 240),
        fast_mode=True,
        twinkle_targets=((123, 0.5),),
    )

    assert captured["twinkle_targets"] == ()
