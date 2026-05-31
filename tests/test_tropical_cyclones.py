from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import pytest
from PySide6.QtCore import QPointF

import zstarview.render.tropical_cyclones as render_tropical_cyclones
from zstarview.tropical_cyclones.cache import (
    TropicalCycloneCacheEntry,
    load_tropical_cyclone_cache,
    save_tropical_cyclone_cache,
)
from zstarview.tropical_cyclones.client import _parse_point
from zstarview.data.import_overture_buildings import iter_download_features
from zstarview.tropical_cyclones.models import (
    TropicalCyclonePoint,
    TropicalCyclonePolygon,
    TropicalCycloneSnapshot,
    project_tropical_cyclone_snapshot,
)
from zstarview.types import ScreenGeometry


def _snapshot() -> TropicalCycloneSnapshot:
    return TropicalCycloneSnapshot(
        storm_name="Jangmi",
        basin="WP",
        advdate_utc=datetime(2026, 5, 30, 2, 10, tzinfo=timezone.utc),
        observed_position=TropicalCyclonePoint(
            lat_deg=12.3,
            lon_deg=145.6,
            valid_time_utc=datetime(2026, 5, 30, 2, 0, tzinfo=timezone.utc),
            label="current",
            maxwind_kt=35.0,
        ),
        forecast_positions=(
            TropicalCyclonePoint(
                lat_deg=13.0,
                lon_deg=146.0,
                valid_time_utc=datetime(2026, 5, 30, 14, 0, tzinfo=timezone.utc),
                label="+12h",
                tau_hr=12,
                maxwind_kt=40.0,
            ),
        ),
        wind_polygons=(
            TropicalCyclonePolygon(
                layer_id=7,
                name="Tropical Storm Force (34kts)",
                rings=(
                    ((12.0, 145.0), (12.0, 146.0), (13.0, 146.0), (13.0, 145.0), (12.0, 145.0)),
                ),
            ),
        ),
        source_url="https://example.invalid/service",
        service_name="Active Hurricanes",
        refreshed_at_utc=datetime(2026, 5, 30, 2, 20, tzinfo=timezone.utc),
        current_storm_id="WP012026",
    )


def test_tropical_cyclone_snapshot_roundtrip() -> None:
    snapshot = _snapshot()
    loaded = TropicalCycloneSnapshot.from_dict(snapshot.to_dict())

    assert loaded == snapshot


def test_tropical_cyclone_cache_roundtrip(tmp_path) -> None:
    snapshot = _snapshot()
    entry = TropicalCycloneCacheEntry(
        snapshot=snapshot,
        cached_at_utc=datetime(2026, 5, 30, 2, 30, tzinfo=timezone.utc),
    )

    save_tropical_cyclone_cache(entry, cache_root=tmp_path)
    loaded = load_tropical_cyclone_cache(cache_root=tmp_path)

    assert loaded == entry
    assert loaded is not None
    assert loaded.cache_version >= 1


def test_tropical_cyclone_cache_entry_loads_legacy_payload() -> None:
    snapshot = _snapshot().to_dict()
    loaded = TropicalCycloneCacheEntry.from_dict(
        {
            "cached_at_utc": "2026-05-30T02:30:00Z",
            "snapshot": snapshot,
        }
    )

    assert loaded is not None
    assert loaded.cache_version == 0
    assert loaded.snapshot == _snapshot()


def test_tropical_cyclone_projection_moves_snapshot_forward() -> None:
    snapshot = _snapshot()
    projected = project_tropical_cyclone_snapshot(
        snapshot,
        datetime(2026, 5, 30, 8, 0, tzinfo=timezone.utc),
    )

    assert projected.observed_position.lat_deg == 12.65
    assert projected.observed_position.lon_deg == 145.8
    assert projected.refreshed_at_utc == datetime(2026, 5, 30, 8, 0, tzinfo=timezone.utc)
    lat_deg, lon_deg = projected.wind_polygons[0].rings[0][0]
    assert lat_deg == pytest.approx(12.35)
    assert lon_deg == pytest.approx(145.2)


def test_tropical_cyclone_projection_uses_advdate_and_tau_when_valid_times_missing() -> None:
    snapshot = TropicalCycloneSnapshot(
        storm_name="Jangmi",
        basin="WP",
        advdate_utc=datetime(2026, 5, 30, 0, 0, tzinfo=timezone.utc),
        observed_position=TropicalCyclonePoint(
            lat_deg=12.3,
            lon_deg=145.6,
        ),
        forecast_positions=(
            TropicalCyclonePoint(
                lat_deg=12.8,
                lon_deg=145.9,
                tau_hr=12,
            ),
        ),
        wind_polygons=(),
    )
    projected = project_tropical_cyclone_snapshot(
        snapshot,
        datetime(2026, 5, 30, 6, 0, tzinfo=timezone.utc),
    )

    assert projected.observed_position.lat_deg == pytest.approx(12.55)
    assert projected.observed_position.lon_deg == pytest.approx(145.75)


def test_missing_overture_download_is_treated_as_empty() -> None:
    assert iter_download_features(
        Path("/tmp/does-not-exist.geojsonseq"),
        fmt="geojsonseq",
    ) == ()


def test_forecast_fldatelbl_is_parsed_as_utc_time() -> None:
    advdate_utc = datetime(2026, 5, 30, 0, 0, tzinfo=timezone.utc)
    point = _parse_point(
        {
            "attributes": {
                "VALIDTIME": None,
                "ADVDATE": int(advdate_utc.timestamp() * 1000.0),
                "FLDATELBL": "2026-05-30 06:00 AM Sat UTC",
                "FCSTPRD": 12,
                "MAXWIND": 35,
            },
            "geometry": {
                "x": 147.7,
                "y": 36.8,
            },
        },
        label_field="FLDATELBL",
    )

    assert point is not None
    assert point.valid_time_utc == datetime(2026, 5, 30, 6, 0, tzinfo=timezone.utc)


def test_forecast_fldatelbl_is_preferred_over_validtime() -> None:
    point = _parse_point(
        {
            "attributes": {
                "VALIDTIME": int(datetime(2026, 5, 30, 7, 40, tzinfo=timezone.utc).timestamp() * 1000.0),
                "FLDATELBL": "2026-05-31 06:00 AM Sun UTC",
                "FCSTPRD": 24,
                "MAXWIND": 35,
            },
            "geometry": {
                "x": 147.7,
                "y": 36.8,
            },
        },
        label_field="FLDATELBL",
    )

    assert point is not None
    assert point.valid_time_utc == datetime(2026, 5, 31, 6, 0, tzinfo=timezone.utc)


def test_tropical_cyclone_projection_is_limited_by_distance_km(monkeypatch) -> None:
    viewer = SimpleNamespace(
        lat_deg=36.75,
        lon_deg=147.65,
        ground_elevation_m=0.0,
        view_center=(45.0, 180.0),
        content_fov_deg=110.0,
        edge_fov_deg=95.0,
    )

    def _fake_project(*_args, **_kwargs):
        return [SimpleNamespace(alt_deg=10.0, az_deg=20.0, distance_km=128.0)]

    monkeypatch.setattr(
        render_tropical_cyclones,
        "project_place_targets_to_altaz",
        _fake_project,
    )
    point = render_tropical_cyclones._project_point(
        36.8,
        147.7,
        viewer=viewer,
        height_m=0.0,
    )

    assert point is not None

    def _too_far(*_args, **_kwargs):
        return [SimpleNamespace(alt_deg=10.0, az_deg=20.0, distance_km=400.0001)]

    monkeypatch.setattr(
        render_tropical_cyclones,
        "project_place_targets_to_altaz",
        _too_far,
    )
    assert (
        render_tropical_cyclones._project_point(
            36.8,
            147.7,
            viewer=viewer,
            height_m=0.0,
        )
        is None
    )


def test_tropical_cyclone_far_marker_projects_without_cutoff(monkeypatch) -> None:
    viewer = SimpleNamespace(
        lat_deg=36.75,
        lon_deg=147.65,
        ground_elevation_m=0.0,
        view_center=(45.0, 180.0),
        content_fov_deg=110.0,
        edge_fov_deg=95.0,
    )

    def _too_far(*_args, **_kwargs):
        return [SimpleNamespace(alt_deg=10.0, az_deg=20.0, distance_km=400.0001)]

    monkeypatch.setattr(
        render_tropical_cyclones,
        "project_place_targets_to_altaz",
        _too_far,
    )
    point = render_tropical_cyclones._project_point_no_cutoff(
        36.8,
        147.7,
        viewer=viewer,
        height_m=0.0,
    )

    assert point is not None
    assert point.distance_km == pytest.approx(400.0001)


def test_tropical_cyclone_far_marker_polygon_is_triangle() -> None:
    point = render_tropical_cyclones._RenderPoint(
        nx=0.25,
        ny=-0.1,
        alt_deg=10.0,
        az_deg=20.0,
        distance_km=450.0,
    )
    geometry = ScreenGeometry(center=(200, 100), radius=80)

    polygon = render_tropical_cyclones._far_marker_polygon(point, geometry=geometry)

    assert polygon.count() == 3
    assert polygon[0].y() == polygon[1].y()
    assert polygon[2].y() > polygon[0].y()
    assert polygon[2].x() == pytest.approx((polygon[0].x() + polygon[1].x()) / 2.0)


def test_tropical_cyclone_filled_marker_tips_at_contact_point() -> None:
    point = render_tropical_cyclones._RenderPoint(
        nx=0.25,
        ny=-0.1,
        alt_deg=10.0,
        az_deg=20.0,
        distance_km=120.0,
    )
    geometry = ScreenGeometry(center=(200, 100), radius=80)

    polygon = render_tropical_cyclones._filled_marker_polygon(point, geometry=geometry)

    assert polygon.count() == 3
    assert polygon[2].x() == pytest.approx((polygon[0].x() + polygon[1].x()) / 2.0)
    assert polygon[2].y() > polygon[0].y()
    assert polygon[2].y() > polygon[1].y()


def test_tropical_cyclone_far_label_uses_lower_alpha() -> None:
    assert render_tropical_cyclones.TROPICAL_CYCLONE_FAR_LABEL_RGBA[3] == 153


def test_tropical_cyclone_falls_back_to_one_pixel_line_when_cylinder_fails(
    monkeypatch,
) -> None:
    viewer = SimpleNamespace(
        lat_deg=36.75,
        lon_deg=147.65,
        ground_elevation_m=0.0,
        view_center=(45.0, 180.0),
        content_fov_deg=110.0,
        edge_fov_deg=95.0,
    )
    snapshot = TropicalCycloneSnapshot(
        storm_name="Jangmi",
        basin="WP",
        advdate_utc=datetime(2026, 5, 30, 2, 10, tzinfo=timezone.utc),
        observed_position=TropicalCyclonePoint(
            lat_deg=12.3,
            lon_deg=145.6,
            valid_time_utc=datetime(2026, 5, 30, 2, 0, tzinfo=timezone.utc),
        ),
        wind_polygons=(),
    )
    calls: list[tuple[str, object]] = []

    def _fake_center(*_args, **_kwargs):
        return render_tropical_cyclones._RenderPoint(
            nx=0.1,
            ny=0.2,
            alt_deg=10.0,
            az_deg=20.0,
            distance_km=399.0,
        )

    def _fake_top(*_args, **_kwargs):
        return render_tropical_cyclones._RenderPoint(
            nx=0.2,
            ny=0.4,
            alt_deg=15.0,
            az_deg=20.0,
            distance_km=399.5,
        )

    def _fake_cylinder(*_args, **_kwargs):
        return None

    def _fake_line(painter, start, end, **kwargs):
        calls.append(("line", (start, end, kwargs)))

    def _fake_far_marker(*_args, **_kwargs):
        calls.append(("far", None))
        return QPointF(12.0, 34.0)

    monkeypatch.setattr(render_tropical_cyclones, "_project_point_no_cutoff", _fake_top)
    monkeypatch.setattr(render_tropical_cyclones, "_project_cyclone_cylinder", _fake_cylinder)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_line", _fake_line)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_far_cyclone_marker", _fake_far_marker)

    class _FakePainter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setRenderHint(self, *_args, **_kwargs) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawText(self, *_args, **_kwargs) -> None:
            pass

    painter = _FakePainter()
    render_tropical_cyclones.draw_tropical_cyclone_overlay(
        painter,
        geometry=ScreenGeometry(center=(100, 100), radius=80),
        viewer=viewer,
        snapshot=snapshot,
        when_utc=datetime(2026, 5, 30, 2, 30, tzinfo=timezone.utc),
        theme=SimpleNamespace(),
        opacity=0.4,
        enabled=True,
    )

    assert calls and calls[0][0] == "line"
    assert not any(kind == "far" for kind, _ in calls)


def test_tropical_cyclone_draws_filled_marker_at_base_when_in_range(monkeypatch) -> None:
    viewer = SimpleNamespace(
        lat_deg=36.75,
        lon_deg=147.65,
        ground_elevation_m=0.0,
        view_center=(45.0, 180.0),
        content_fov_deg=110.0,
        edge_fov_deg=95.0,
    )
    snapshot = TropicalCycloneSnapshot(
        storm_name="Jangmi",
        basin="WP",
        advdate_utc=datetime(2026, 5, 30, 2, 10, tzinfo=timezone.utc),
        observed_position=TropicalCyclonePoint(
            lat_deg=12.3,
            lon_deg=145.6,
            valid_time_utc=datetime(2026, 5, 30, 2, 0, tzinfo=timezone.utc),
        ),
        wind_polygons=(),
    )
    calls: list[tuple[str, object]] = []

    def _fake_center(*_args, **_kwargs):
        return render_tropical_cyclones._RenderPoint(
            nx=0.1,
            ny=0.2,
            alt_deg=10.0,
            az_deg=20.0,
            distance_km=100.0,
        )

    def _fake_cylinder(*_args, **_kwargs):
        return render_tropical_cyclones._CycloneCylinderGeometry(
            top_center_point=render_tropical_cyclones._RenderPoint(
                nx=0.2,
                ny=0.4,
                alt_deg=15.0,
                az_deg=20.0,
                distance_km=100.0,
            ),
            ring_points_by_height=(
                (
                    render_tropical_cyclones._RenderPoint(
                        nx=0.0,
                        ny=0.0,
                        alt_deg=0.0,
                        az_deg=0.0,
                        distance_km=100.0,
                    ),
                    render_tropical_cyclones._RenderPoint(
                        nx=0.1,
                        ny=0.0,
                        alt_deg=0.0,
                        az_deg=0.0,
                        distance_km=100.0,
                    ),
                    render_tropical_cyclones._RenderPoint(
                        nx=0.0,
                        ny=0.1,
                        alt_deg=0.0,
                        az_deg=0.0,
                        distance_km=100.0,
                    ),
                ),
                (
                    render_tropical_cyclones._RenderPoint(
                        nx=0.0,
                        ny=0.0,
                        alt_deg=5.0,
                        az_deg=0.0,
                        distance_km=100.0,
                    ),
                    render_tropical_cyclones._RenderPoint(
                        nx=0.1,
                        ny=0.0,
                        alt_deg=5.0,
                        az_deg=0.0,
                        distance_km=100.0,
                    ),
                    render_tropical_cyclones._RenderPoint(
                        nx=0.0,
                        ny=0.1,
                        alt_deg=5.0,
                        az_deg=0.0,
                        distance_km=100.0,
                    ),
                ),
                (
                    render_tropical_cyclones._RenderPoint(
                        nx=0.0,
                        ny=0.0,
                        alt_deg=10.0,
                        az_deg=0.0,
                        distance_km=100.0,
                    ),
                    render_tropical_cyclones._RenderPoint(
                        nx=0.1,
                        ny=0.0,
                        alt_deg=10.0,
                        az_deg=0.0,
                        distance_km=100.0,
                    ),
                    render_tropical_cyclones._RenderPoint(
                        nx=0.0,
                        ny=0.1,
                        alt_deg=10.0,
                        az_deg=0.0,
                        distance_km=100.0,
                    ),
                ),
                (
                    render_tropical_cyclones._RenderPoint(
                        nx=0.0,
                        ny=0.0,
                        alt_deg=15.0,
                        az_deg=0.0,
                        distance_km=100.0,
                    ),
                    render_tropical_cyclones._RenderPoint(
                        nx=0.1,
                        ny=0.0,
                        alt_deg=15.0,
                        az_deg=0.0,
                        distance_km=100.0,
                    ),
                    render_tropical_cyclones._RenderPoint(
                        nx=0.0,
                        ny=0.1,
                        alt_deg=15.0,
                        az_deg=0.0,
                        distance_km=100.0,
                    ),
                ),
            ),
            radius_km=0.5,
            height_km=15.0,
        )

    def _fake_filled_marker(*_args, **_kwargs):
        calls.append(("filled", None))

    def _fake_far_marker(*_args, **_kwargs):
        calls.append(("far", None))
        return QPointF(12.0, 34.0)

    def _fake_line(*_args, **_kwargs):
        calls.append(("line", None))

    monkeypatch.setattr(render_tropical_cyclones, "_project_point_no_cutoff", _fake_center)
    monkeypatch.setattr(render_tropical_cyclones, "_project_cyclone_cylinder", _fake_cylinder)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_filled_cyclone_marker", _fake_filled_marker)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_far_cyclone_marker", _fake_far_marker)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_line", _fake_line)

    class _FakePainter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setRenderHint(self, *_args, **_kwargs) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawText(self, *_args, **_kwargs) -> None:
            pass

        def drawPolygon(self, *_args, **_kwargs) -> None:
            pass

        def drawLine(self, *_args, **_kwargs) -> None:
            pass

    painter = _FakePainter()
    render_tropical_cyclones.draw_tropical_cyclone_overlay(
        painter,
        geometry=ScreenGeometry(center=(100, 100), radius=80),
        viewer=viewer,
        snapshot=snapshot,
        when_utc=datetime(2026, 5, 30, 2, 30, tzinfo=timezone.utc),
        theme=SimpleNamespace(),
        opacity=0.4,
        enabled=True,
    )

    assert ("filled", None) in calls
    assert not any(kind == "far" for kind, _ in calls)


def test_tropical_cyclone_side_radius_is_fixed() -> None:
    assert render_tropical_cyclones._cyclone_side_radius_km() == pytest.approx(0.5)


def test_tropical_cyclone_cylinder_projects_top_and_base(monkeypatch) -> None:
    viewer = SimpleNamespace(
        lat_deg=36.75,
        lon_deg=147.65,
        ground_elevation_m=0.0,
        view_center=(45.0, 180.0),
        content_fov_deg=110.0,
        edge_fov_deg=95.0,
    )
    calls: list[float] = []

    def _fake_project(*_args, **kwargs):
        height_value = kwargs["target_height_m"]
        if isinstance(height_value, list):
            height_m = float(height_value[0])
        else:
            height_m = float(height_value)
        calls.append(height_m)
        return [
            SimpleNamespace(
                alt_deg=10.0 + (height_m / 10_000.0),
                az_deg=20.0,
                distance_km=120.0 + (height_m / 100_000.0),
            )
        ]

    monkeypatch.setattr(
        render_tropical_cyclones,
        "project_place_targets_to_altaz",
        _fake_project,
    )

    cylinder = render_tropical_cyclones._project_cyclone_cylinder(
        36.8,
        147.7,
        viewer=viewer,
        maxwind_kt=90.0,
    )

    assert cylinder is not None
    assert calls[0] == 0.0
    assert calls[1] == 0.0
    assert calls.count(0.0) == 10
    assert calls.count(5000.0) == 9
    assert calls.count(10000.0) == 9
    assert calls.count(15000.0) == 9
    assert len(calls) == 37
    assert len(cylinder.ring_points_by_height) == 4
    assert all(
        len(ring_points) == render_tropical_cyclones.TROPICAL_CYCLONE_CYLINDER_RING_SAMPLES
        for ring_points in cylinder.ring_points_by_height
    )
    assert cylinder.radius_km == pytest.approx(render_tropical_cyclones._cyclone_side_radius_km())
    assert cylinder.height_km == pytest.approx(render_tropical_cyclones.TROPICAL_CYCLONE_CYLINDER_HEIGHT_KM)


def test_tropical_cyclone_cylinder_side_radius_is_fixed(monkeypatch) -> None:
    viewer = SimpleNamespace(
        lat_deg=36.75,
        lon_deg=147.65,
        ground_elevation_m=0.0,
        view_center=(45.0, 180.0),
        content_fov_deg=110.0,
        edge_fov_deg=95.0,
    )

    def _fake_project(*_args, **kwargs):
        height_value = kwargs["target_height_m"]
        if isinstance(height_value, list):
            height_m = float(height_value[0])
        else:
            height_m = float(height_value)
        return [
            SimpleNamespace(
                alt_deg=10.0 + (height_m / 10_000.0),
                az_deg=20.0,
                distance_km=120.0 + (height_m / 100_000.0),
            )
        ]

    monkeypatch.setattr(
        render_tropical_cyclones,
        "project_place_targets_to_altaz",
        _fake_project,
    )

    weak_cylinder = render_tropical_cyclones._project_cyclone_cylinder(
        36.8,
        147.7,
        viewer=viewer,
        maxwind_kt=35.0,
    )
    strong_cylinder = render_tropical_cyclones._project_cyclone_cylinder(
        36.8,
        147.7,
        viewer=viewer,
        maxwind_kt=120.0,
    )

    assert weak_cylinder is not None
    assert strong_cylinder is not None
    assert weak_cylinder.radius_km == pytest.approx(0.5)
    assert strong_cylinder.radius_km == pytest.approx(0.5)


def test_tropical_cyclone_convex_hull_discards_interior_points() -> None:
    from PySide6.QtCore import QPointF

    hull = render_tropical_cyclones._convex_hull(
        [
            QPointF(0.0, 0.0),
            QPointF(2.0, 0.0),
            QPointF(2.0, 2.0),
            QPointF(0.0, 2.0),
            QPointF(1.0, 1.0),
        ]
    )

    assert len(hull) == 4
    assert {(p.x(), p.y()) for p in hull} == {
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
    }
