import numpy as np
from PySide6.QtGui import QImage, QPainter

import zstarview.gui.composite as render_composite
from zstarview.render.earth_guide import (
    EARTH_GUIDE_FOREGROUND_WIDTH,
    EARTH_GUIDE_UNDERLAY_WIDTH,
    EarthGuideRing,
    _effective_visible_altitude_limit_deg,
    _observer_dead_zone_km,
    _observer_visible_altitude_limit_deg,
    earth_guide_line_alpha,
    earth_guide_underlay_line_alpha,
    draw_earth_guide,
    load_earth_guide_rings,
)
from zstarview.render.qt_image import qimage_to_np_rgba
from zstarview.types import ScreenGeometry


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
            view_center=(0.0, 180.0),
            observer_lat_deg=35.68,
            observer_lon_deg=139.76,
            observer_height_m=635.0,
            terrain_profile_altaz=None,
            earth_guide_opacity=0.028,
            content_fov_deg=100.0,
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

    def save(self) -> None:
        pass

    def restore(self) -> None:
        pass

    def setRenderHint(self, *_args, **_kwargs) -> None:  # noqa: N802 - Qt naming
        pass

    def setPen(self, *_args, **_kwargs) -> None:
        pass

    def drawPolyline(self, polyline) -> None:
        self.polylines.append(polyline)


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

    monkeypatch.setattr("zstarview.render.earth_guide.load_earth_guide_rings", lambda path_str=None: rings)
    monkeypatch.setattr(
        "zstarview.render.earth_guide._ring_fragments_altaz",
        lambda ring, **kwargs: seen.append(ring.source_name) or [[(0.0, 0.0), (1.0, 1.0)]],
    )

    painter = _DummyPainter()
    draw_earth_guide(
        painter,
        geometry=ScreenGeometry(center=(120, 120), radius=100),
        view_center=(0.0, 180.0),
        observer_lat_deg=35.68,
        observer_lon_deg=139.76,
        observer_height_m=635.0,
        terrain_profile_altaz=None,
        earth_guide_opacity=0.028,
        content_fov_deg=100.0,
    )
    polyline_counts.append(len(painter.polylines))
    assert seen == ["ring-0", "ring-1", "ring-2", "ring-3"]
    assert polyline_counts == [16]

    seen.clear()
    painter = _DummyPainter()
    draw_earth_guide(
        painter,
        geometry=ScreenGeometry(center=(120, 120), radius=100),
        view_center=(0.0, 180.0),
        observer_lat_deg=35.68,
        observer_lon_deg=139.76,
        observer_height_m=635.0,
        terrain_profile_altaz=None,
        earth_guide_opacity=0.028,
        content_fov_deg=100.0,
        fast_mode=True,
    )
    assert seen == ["ring-0", "ring-2"]
    assert len(painter.polylines) == 2


def test_overlay_earth_guide_forwards_fast_mode(monkeypatch) -> None:
    calls: list[bool] = []

    monkeypatch.setattr(
        render_composite,
        "draw_earth_guide",
        lambda _painter, **kwargs: calls.append(bool(kwargs["fast_mode"])),
    )

    image = QImage(16, 16, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    out = render_composite._overlay_earth_guide(
        image,
        geometry=ScreenGeometry(center=(8, 8), radius=6),
        view_center=(0.0, 180.0),
        observer_lat_deg=35.68,
        observer_lon_deg=139.76,
        observer_height_m=635.0,
        terrain_profile_altaz=None,
        earth_guide_opacity=0.028,
        content_fov_deg=100.0,
        fast_mode=True,
    )

    assert calls == [True]
    assert out.size() == image.size()
