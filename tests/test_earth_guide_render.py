import numpy as np
import pytest
from PySide6.QtGui import QImage, QPainter

import zstarview.gui.composite as render_composite
from zstarview.render.earth_guide import (
    EARTH_GUIDE_FOREGROUND_WIDTH,
    EARTH_GUIDE_UNDERLAY_WIDTH,
    EarthGuideRing,
    _build_ring_fill_points,
    _earth_guide_underlay_pass_specs,
    _effective_visible_altitude_limit_deg,
    _observer_dead_zone_km,
    _observer_visible_altitude_limit_deg,
    draw_earth_guide,
    earth_guide_line_alpha,
    earth_guide_underlay_line_alpha,
    load_earth_guide_rings,
)
from zstarview.render.qt_image import qimage_to_np_rgba
from zstarview.types import ScreenGeometry, ViewerData


def _viewer(
    *,
    view_center: tuple[float, float] = (0.0, 180.0),
    observer_lat_deg: float = 35.68,
    observer_lon_deg: float = 139.76,
    observer_height_m: float = 635.0,
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


def test_load_earth_guide_rings_has_expected_runtime_payload() -> None:
    rings = load_earth_guide_rings()

    assert len(rings) >= 10
    assert sum(len(ring.points_lonlat_deg) for ring in rings) > 200
    assert all(ring.points_xyz.shape[1] == 3 for ring in rings)


def test_draw_earth_guide_renders_visible_lines_below_horizon() -> None:
    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        draw_earth_guide(
            painter,
            geometry=ScreenGeometry(center=(120, 120), radius=100),
            viewer_data=_viewer(),
            terrain_profile_altaz=None,
            earth_guide_opacity=0.028,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    alpha = arr[..., 3]

    assert int(np.count_nonzero(alpha)) > 20
    top_half = int(np.count_nonzero(alpha[:120, :]))
    bottom_half = int(np.count_nonzero(alpha[120:, :]))
    assert bottom_half > top_half


def test_earth_guide_underlay_is_thicker_and_softer() -> None:
    opacity = 0.028

    assert EARTH_GUIDE_UNDERLAY_WIDTH > EARTH_GUIDE_FOREGROUND_WIDTH
    assert earth_guide_underlay_line_alpha(opacity) < earth_guide_line_alpha(opacity)


def test_draw_earth_guide_single_line_skips_fill_and_underlay(monkeypatch) -> None:
    monkeypatch.setattr(
        "zstarview.render.earth_guide._draw_fill_segments_for_ring",
        lambda *_args, **_kwargs: pytest.fail("Atlas must not draw Earth guide fill lines"),
    )
    monkeypatch.setattr(
        "zstarview.render.earth_guide._earth_guide_underlay_pass_specs",
        lambda *_args, **_kwargs: pytest.fail("Atlas must not draw Earth guide underlay lines"),
    )
    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        draw_earth_guide(
            painter,
            geometry=ScreenGeometry(center=(120, 120), radius=100),
            viewer_data=_viewer(),
            earth_guide_opacity=0.028,
            single_line=True,
        )
    finally:
        painter.end()


def test_earth_guide_uses_geometric_horizon_by_default() -> None:
    assert _effective_visible_altitude_limit_deg(180.0, None) == 0.0
    assert _effective_visible_altitude_limit_deg(180.0, [(12.0, 180.0)]) == 12.0


def test_earth_guide_dead_zone_and_altitude_clip_scale_with_height() -> None:
    assert _observer_dead_zone_km(0.0) == 20.0
    assert _observer_dead_zone_km(635.0) > _observer_dead_zone_km(0.0)
    assert _observer_visible_altitude_limit_deg(180.0, 0.0, None) == -1.0
    assert _observer_visible_altitude_limit_deg(180.0, 635.0, None) < -1.0
    assert _observer_visible_altitude_limit_deg(180.0, 635.0, [(12.0, 180.0)]) < -1.0
    assert _observer_visible_altitude_limit_deg(180.0, 0.0, [(-3.0, 180.0)]) == -3.0


class _DummyPainter:
    def __init__(self) -> None:
        self.polylines: list[object] = []
        self.lines: list[object] = []
        self.pen_widths: list[float] = []
        self.pen_alphas: list[float] = []

    def save(self) -> None:
        pass

    def restore(self) -> None:
        pass

    def setRenderHint(self, *_args, **_kwargs) -> None:  # noqa: N802 - Qt naming
        pass

    def setPen(self, pen, *_args, **_kwargs) -> None:
        self.pen_widths.append(float(pen.widthF()))
        self.pen_alphas.append(float(pen.color().alphaF()))

    def setBrush(self, *_args, **_kwargs) -> None:
        pass

    def drawPolyline(self, polyline) -> None:
        self.polylines.append(polyline)

    def drawLine(self, *args) -> None:
        self.lines.append(args)


def test_build_ring_fill_points_generates_points_for_simple_landmass() -> None:
    ring = EarthGuideRing(
        source_name="test",
        label_name=None,
        points_lonlat_deg=np.asarray(
            [(-10.0, -10.0), (10.0, -10.0), (10.0, 10.0), (-10.0, 10.0)],
            dtype=np.float64,
        ),
        points_xyz=np.asarray(
            [
                (-10.0, -10.0, 0.0),
                (10.0, -10.0, 0.0),
                (10.0, 10.0, 0.0),
                (-10.0, 10.0, 0.0),
            ],
            dtype=np.float64,
        ),
        approx_area_deg2=400.0,
    )

    lonlat, xyz = _build_ring_fill_points(ring)

    assert len(lonlat) > 0
    assert len(xyz) == len(lonlat)
    assert np.all(np.abs(lonlat[:, 0]) <= 10.0)
    assert np.all(np.abs(lonlat[:, 1]) <= 10.0)


def test_build_ring_fill_points_skips_antarctica_like_rings() -> None:
    ring = EarthGuideRing(
        source_name="Antarctica",
        label_name=None,
        points_lonlat_deg=np.asarray(
            [
                (-58.0, -64.0),
                (-65.0, -68.0),
                (-60.0, -74.0),
                (-77.0, -77.0),
                (-73.0, -79.0),
                (-78.0, -83.0),
                (-28.0, -80.0),
                (-35.0, -78.0),
                (-6.0, -70.0),
                (38.0, -69.0),
                (54.0, -65.0),
                (68.0, -67.0),
            ],
            dtype=np.float64,
        ),
        points_xyz=np.asarray([(0.0, 0.0, 1.0)] * 12, dtype=np.float64),
        approx_area_deg2=6008.0,
    )

    lonlat, xyz = _build_ring_fill_points(ring)

    assert len(lonlat) == 0
    assert len(xyz) == 0


def test_draw_earth_guide_fast_mode_subsamples_rings(monkeypatch) -> None:
    rings = tuple(
        EarthGuideRing(
            source_name=f"ring-{idx}",
            label_name=None,
            points_lonlat_deg=np.asarray(
                [(-10.0 + idx, -5.0), (0.0 + idx, 0.0), (10.0 + idx, 5.0)],
                dtype=np.float64,
            ),
            points_xyz=np.asarray(
                [(1.0, 0.0, 0.0), (0.5, 0.5, 0.7071), (0.0, 1.0, 0.0)],
                dtype=np.float64,
            ),
        )
        for idx in range(4)
    )

    seen: list[str] = []
    polyline_counts: list[int] = []

    monkeypatch.setattr("zstarview.render.earth_guide.load_earth_guide_rings", lambda: rings)
    monkeypatch.setattr(
        "zstarview.render.earth_guide._ring_fragments_altaz",
        lambda ring, **kwargs: seen.append(ring.source_name) or [[(0.0, 0.0), (1.0, 1.0)]],
    )

    painter = _DummyPainter()
    draw_earth_guide(
        painter,
        geometry=ScreenGeometry(center=(120, 120), radius=100),
        viewer_data=_viewer(
            view_center=(0.0, 180.0),
            observer_lat_deg=35.68,
            observer_lon_deg=139.76,
            observer_height_m=635.0,
            edge_fov_deg=90.0,
            content_fov_deg=100.0,
        ),
        terrain_profile_altaz=None,
        earth_guide_opacity=0.028,
    )
    polyline_counts.append(len(painter.polylines))
    assert seen == ["ring-0", "ring-1", "ring-2", "ring-3"]
    assert polyline_counts == [16]

    seen.clear()
    painter = _DummyPainter()
    draw_earth_guide(
        painter,
        geometry=ScreenGeometry(center=(120, 120), radius=100),
        viewer_data=_viewer(),
        terrain_profile_altaz=None,
        earth_guide_opacity=0.028,
        fast_mode=True,
    )
    assert seen == ["ring-0", "ring-2"]
    assert len(painter.polylines) == 2


def test_draw_earth_guide_fast_mode_skips_fill_lines(monkeypatch) -> None:
    ring = EarthGuideRing(
        source_name="ring-fill-fast",
        label_name=None,
        points_lonlat_deg=np.asarray(
            [(-20.0, -10.0), (20.0, -10.0), (20.0, 10.0), (-20.0, 10.0)],
            dtype=np.float64,
        ),
        points_xyz=np.asarray(
            [(-20.0, -10.0, 0.0), (20.0, -10.0, 0.0), (20.0, 10.0, 0.0), (-20.0, 10.0, 0.0)],
            dtype=np.float64,
        ),
        approx_area_deg2=800.0,
        fill_points_lonlat_deg=np.asarray(
            [(0.0, 0.0), (5.0, 0.0), (-5.0, 0.0)],
            dtype=np.float64,
        ),
        fill_points_xyz=np.asarray(
            [(0.5, 0.0, 0.0), (0.49, 0.05, 0.0), (0.49, -0.05, 0.0)],
            dtype=np.float64,
        ),
    )

    fill_calls: list[str] = []

    monkeypatch.setattr("zstarview.render.earth_guide.load_earth_guide_rings", lambda: (ring,))
    monkeypatch.setattr(
        "zstarview.render.earth_guide._ring_fragments_altaz",
        lambda ring, **kwargs: [[(0.0, 0.0), (1.0, 1.0)]],
    )
    monkeypatch.setattr(
        "zstarview.render.earth_guide._draw_fill_segments_for_ring",
        lambda *args, **kwargs: fill_calls.append("called"),
    )

    painter = _DummyPainter()
    draw_earth_guide(
        painter,
        geometry=ScreenGeometry(center=(120, 120), radius=100),
        viewer_data=_viewer(),
        terrain_profile_altaz=None,
        earth_guide_opacity=0.028,
        fast_mode=True,
    )

    assert fill_calls == []
    assert len(painter.polylines) > 0
    assert len(painter.lines) == 0


def test_draw_earth_guide_renders_fill_points_before_outline(monkeypatch) -> None:
    ring = EarthGuideRing(
        source_name="ring-fill",
        label_name=None,
        points_lonlat_deg=np.asarray(
            [(-20.0, -10.0), (20.0, -10.0), (20.0, 10.0), (-20.0, 10.0)],
            dtype=np.float64,
        ),
        points_xyz=np.asarray(
            [(-20.0, -10.0, 0.0), (20.0, -10.0, 0.0), (20.0, 10.0, 0.0), (-20.0, 10.0, 0.0)],
            dtype=np.float64,
        ),
        approx_area_deg2=800.0,
        fill_points_lonlat_deg=np.asarray(
            [(0.0, 0.0), (5.0, 0.0), (-5.0, 0.0)],
            dtype=np.float64,
        ),
        fill_points_xyz=np.asarray(
            [(0.5, 0.0, 0.0), (0.49, 0.05, 0.0), (0.49, -0.05, 0.0)],
            dtype=np.float64,
        ),
    )

    monkeypatch.setattr("zstarview.render.earth_guide.load_earth_guide_rings", lambda: (ring,))
    monkeypatch.setattr(
        "zstarview.render.earth_guide._ring_fragments_altaz",
        lambda ring, **kwargs: [[(0.0, 0.0), (1.0, 1.0)]],
    )

    painter = _DummyPainter()
    draw_earth_guide(
        painter,
        geometry=ScreenGeometry(center=(120, 120), radius=100),
        viewer_data=_viewer(),
        terrain_profile_altaz=None,
        earth_guide_opacity=0.028,
    )

    assert len(painter.lines) > 0
    assert len(painter.polylines) > 0


def test_overlay_earth_guide_forwards_fast_mode(monkeypatch) -> None:
    calls: list[bool] = []

    monkeypatch.setattr(
        render_composite,
        "draw_earth_guide",
        lambda *_args, **_kwargs: calls.append(True),
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
    )

    assert calls == [True]
    assert out.size() == image.size()


def test_draw_earth_guide_scales_fast_mode_lines_with_visibility_boost() -> None:
    painter = _DummyPainter()

    draw_earth_guide(
        painter,
        geometry=ScreenGeometry(center=(120, 120), radius=100),
        viewer_data=_viewer(),
        terrain_profile_altaz=None,
        earth_guide_opacity=0.028,
        visibility_boost=2.0,
        fast_mode=True,
    )

    expected_width = max(0.7, EARTH_GUIDE_FOREGROUND_WIDTH * 0.75)
    expected_alpha = 0.18 + (0.028 * 0.25)
    assert painter.pen_widths
    assert all(width == expected_width for width in painter.pen_widths)
    assert all(abs(alpha - expected_alpha) < 1.0e-3 for alpha in painter.pen_alphas)


def test_draw_earth_guide_boosts_only_thin_underlay_pass(monkeypatch) -> None:
    ring = EarthGuideRing(
        source_name="ring",
        label_name=None,
        points_lonlat_deg=np.asarray([(0.0, 0.0), (1.0, 1.0), (-1.0, 1.0)], dtype=np.float64),
        points_xyz=np.asarray([(0.0, 0.0, 1.0), (0.1, 0.1, 1.0), (-0.1, 0.1, 1.0)], dtype=np.float64),
        approx_area_deg2=1.0,
    )

    monkeypatch.setattr("zstarview.render.earth_guide.load_earth_guide_rings", lambda: (ring,))
    monkeypatch.setattr(
        "zstarview.render.earth_guide._draw_fill_segments_for_ring",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        "zstarview.render.earth_guide._ring_fragments_altaz",
        lambda *args, **kwargs: [[(0.0, 0.0), (1.0, 1.0)]],
    )

    painter = _DummyPainter()
    draw_earth_guide(
        painter,
        geometry=ScreenGeometry(center=(120, 120), radius=100),
        viewer_data=_viewer(),
        terrain_profile_altaz=None,
        earth_guide_opacity=0.028,
        visibility_boost=2.0,
        fast_mode=False,
    )

    base_specs = list(_earth_guide_underlay_pass_specs(0.028))
    assert len(painter.pen_widths) >= 4
    underlay_widths = painter.pen_widths[-4:-1]
    underlay_alphas = painter.pen_alphas[-4:-1]
    assert underlay_widths[0] == pytest.approx(base_specs[0][0])
    assert underlay_widths[1] == pytest.approx(base_specs[1][0])
    assert underlay_widths[2] == pytest.approx(base_specs[2][0] * 2.0)
    assert painter.pen_widths[-1] == pytest.approx(EARTH_GUIDE_FOREGROUND_WIDTH)
    assert abs(underlay_alphas[0] - base_specs[0][1]) < 1.0e-3
    assert abs(underlay_alphas[1] - base_specs[1][1]) < 1.0e-3
    assert abs(underlay_alphas[2] - (base_specs[2][1] * 2.0)) < 1.0e-3
    assert abs(painter.pen_alphas[-1] - earth_guide_line_alpha(0.028)) < 1.0e-3
