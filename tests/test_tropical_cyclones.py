from __future__ import annotations

import logging
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path
from types import SimpleNamespace

import pytest
from PySide6.QtCore import QPointF

import zstarview.gui.tropical_cyclone_controller as tropical_cyclone_controller
import zstarview.render.tropical_cyclones as render_tropical_cyclones
import zstarview.tropical_cyclones.client as cyclone_client
from zstarview.tropical_cyclones.cache import (
    TROPICAL_CYCLONE_CACHE_VERSION,
    TropicalCycloneCacheEntry,
    load_tropical_cyclone_cache,
    save_tropical_cyclone_cache,
)
from zstarview.tropical_cyclones.client import _parse_point
from zstarview.data.import_overture_buildings import iter_download_features
from zstarview.tropical_cyclones.models import (
    TropicalCyclonePoint,
    TropicalCyclonePolygon,
    TropicalCycloneSnapshotCollection,
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


def _collection(*snapshots: TropicalCycloneSnapshot) -> TropicalCycloneSnapshotCollection:
    return TropicalCycloneSnapshotCollection(
        snapshots=tuple(snapshots or (_snapshot(),)),
        source_url="https://example.invalid/service",
        service_name="Active Hurricanes",
        refreshed_at_utc=datetime(2026, 5, 30, 2, 20, tzinfo=timezone.utc),
    )


def test_tropical_cyclone_snapshot_roundtrip() -> None:
    snapshot = _snapshot()
    loaded = TropicalCycloneSnapshot.from_dict(snapshot.to_dict())

    assert loaded == snapshot


def test_tropical_cyclone_summary_text_omits_question_mark_when_advdate_is_unneeded() -> None:
    snapshot = TropicalCycloneSnapshot(
        storm_name="Amanda",
        basin="ep",
        advdate_utc=None,
        observed_position=TropicalCyclonePoint(
            lat_deg=13.0,
            lon_deg=-134.2,
            valid_time_utc=datetime(2026, 6, 6, 0, 0, tzinfo=timezone.utc),
        ),
        forecast_positions=(),
        wind_polygons=(),
        source_url="https://example.invalid/service",
        service_name="FeatureServer",
        refreshed_at_utc=datetime(2026, 6, 6, 4, 50, tzinfo=timezone.utc),
        current_storm_id="1",
    )
    collection = TropicalCycloneSnapshotCollection(
        snapshots=(snapshot,),
        source_url="https://example.invalid/service",
        service_name="FeatureServer",
        refreshed_at_utc=datetime(2026, 6, 6, 4, 50, tzinfo=timezone.utc),
    )

    assert snapshot.advdate_is_required_for_projection() is False
    assert collection.summary_text() == "Amanda"


def test_tropical_cyclone_summary_text_shows_count_only_for_multiple_storms() -> None:
    collection = TropicalCycloneSnapshotCollection(
        snapshots=(
            _snapshot(),
            TropicalCycloneSnapshot(
                storm_name="One-e",
                basin="wp",
                advdate_utc=datetime(2026, 5, 30, 8, 0, tzinfo=timezone.utc),
                observed_position=TropicalCyclonePoint(
                    lat_deg=18.4,
                    lon_deg=132.1,
                    valid_time_utc=datetime(2026, 5, 30, 8, 0, tzinfo=timezone.utc),
                ),
                forecast_positions=(),
                wind_polygons=(),
                source_url="https://example.invalid/service",
                service_name="FeatureServer",
                refreshed_at_utc=datetime(2026, 5, 30, 2, 20, tzinfo=timezone.utc),
                current_storm_id="WP022026",
            ),
        ),
        source_url="https://example.invalid/service",
        service_name="FeatureServer",
        refreshed_at_utc=datetime(2026, 5, 30, 2, 20, tzinfo=timezone.utc),
    )

    assert collection.summary_text() == "2 storms"


def test_tropical_cyclone_cache_roundtrip(tmp_path) -> None:
    collection = _collection(_snapshot())
    entry = TropicalCycloneCacheEntry(
        snapshot_collection=collection,
        cached_at_utc=datetime(2026, 5, 30, 2, 30, tzinfo=timezone.utc),
    )

    save_tropical_cyclone_cache(entry, cache_root=tmp_path)
    loaded = load_tropical_cyclone_cache(cache_root=tmp_path)

    assert loaded == entry
    assert loaded is not None
    assert loaded.cache_version == TROPICAL_CYCLONE_CACHE_VERSION
    assert loaded.snapshot_collection == collection


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
    assert loaded.snapshot_collection.snapshots == (_snapshot(),)


def test_tropical_cyclone_cache_loader_discards_old_version(tmp_path) -> None:
    cache_path = tmp_path / "active_hurricanes.json"
    cache_path.write_text(
        json.dumps(
            {
                "cache_version": 3,
                "cached_at_utc": "2026-05-30T02:30:00Z",
                "snapshot": _snapshot().to_dict(),
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    loaded = load_tropical_cyclone_cache(cache_root=tmp_path)

    assert loaded is None
    assert not cache_path.exists()


def test_tropical_cyclone_snapshot_collection_roundtrip() -> None:
    collection = _collection(
        _snapshot(),
        TropicalCycloneSnapshot(
            storm_name="One-e",
            basin="WP",
            advdate_utc=datetime(2026, 5, 30, 8, 0, tzinfo=timezone.utc),
            observed_position=TropicalCyclonePoint(
                lat_deg=18.4,
                lon_deg=132.1,
                valid_time_utc=datetime(2026, 5, 30, 8, 0, tzinfo=timezone.utc),
            ),
            forecast_positions=(),
            wind_polygons=(),
            source_url="https://example.invalid/service",
            service_name="Active Hurricanes",
            refreshed_at_utc=datetime(2026, 5, 30, 2, 20, tzinfo=timezone.utc),
            current_storm_id="WP022026",
        ),
    )

    loaded = TropicalCycloneSnapshotCollection.from_dict(collection.to_dict())

    assert loaded == collection


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


def test_fetch_active_hurricanes_snapshot_returns_all_active_storms(monkeypatch) -> None:
    obs_time_a = datetime(2026, 5, 30, 2, 0, tzinfo=timezone.utc)
    obs_time_b = datetime(2026, 5, 30, 8, 0, tzinfo=timezone.utc)
    advdate_a = datetime(2026, 5, 30, 0, 0, tzinfo=timezone.utc)
    advdate_b = datetime(2026, 5, 30, 6, 0, tzinfo=timezone.utc)

    def _feature(storm_name: str, basin: str, storm_id: str, dtg: datetime, advdate: datetime, x: float, y: float) -> dict:
        return {
            "attributes": {
                "STORMNAME": storm_name,
                "BASIN": basin,
                "ATCFID": storm_id,
                "DTG": int(dtg.timestamp() * 1000.0),
                "ADVDATE": int(advdate.timestamp() * 1000.0),
                "DATELBL": dtg.strftime("%Y-%m-%d %H:%M UTC"),
                "FCSTPRD": 0,
            },
            "geometry": {
                "x": x,
                "y": y,
            },
        }

    observed_features = [
        _feature("One-e", "WP", "WP022026", obs_time_b, advdate_b, 132.1, 18.4),
        _feature("Jangmi", "WP", "WP012026", obs_time_a, advdate_a, 145.6, 12.3),
    ]
    forecast_features = [
        {
            "attributes": {
                "STORMNAME": "One-e",
                "BASIN": "WP",
                "ATCFID": "WP022026",
                "ADVDATE": int(advdate_b.timestamp() * 1000.0),
                "VALIDTIME": int((advdate_b + timedelta(hours=12)).timestamp() * 1000.0),
                "FLDATELBL": "2026-05-30 06:00 PM Sat UTC",
                "TAU": 12,
                "FCSTPRD": 12,
                "MAXWIND": 40,
            },
            "geometry": {
                "x": 133.0,
                "y": 18.9,
            },
        },
        {
            "attributes": {
                "STORMNAME": "Jangmi",
                "BASIN": "WP",
                "ATCFID": "WP012026",
                "ADVDATE": int(advdate_a.timestamp() * 1000.0),
                "VALIDTIME": int((advdate_a + timedelta(hours=12)).timestamp() * 1000.0),
                "FLDATELBL": "2026-05-30 12:00 PM Sat UTC",
                "TAU": 12,
                "FCSTPRD": 12,
                "MAXWIND": 35,
            },
            "geometry": {
                "x": 146.0,
                "y": 13.0,
            },
        },
    ]
    wind_feature = {
        "attributes": {
            "STORMNAME": "Jangmi",
            "BASIN": "WP",
            "ATCFID": "WP012026",
            "ADVDATE": int(advdate_a.timestamp() * 1000.0),
        },
        "geometry": {
            "rings": [
                [[145.0, 12.0], [146.0, 12.0], [146.0, 13.0], [145.0, 13.0], [145.0, 12.0]]
            ],
        },
    }
    wind_feature_b = {
        "attributes": {
            "STORMNAME": "One-e",
            "BASIN": "WP",
            "ATCFID": "WP022026",
            "ADVDATE": int(advdate_b.timestamp() * 1000.0),
        },
        "geometry": {
            "rings": [
                [[132.0, 18.0], [133.0, 18.0], [133.0, 19.0], [132.0, 19.0], [132.0, 18.0]]
            ],
        },
    }

    def _fake_fetch_json(url: str, *, timeout_s: float, user_agent: str) -> dict[str, object]:
        if url.endswith("?f=json"):
            return {
                "name": "Active Hurricanes",
                "layers": [
                    {"id": 0, "name": "Forecast Position"},
                    {"id": 1, "name": "Observed Position"},
                    {"id": 7, "name": "Tropical Storm Force (34kts)"},
                    {"id": 8, "name": "Strong Tropical Storm (50kts)"},
                    {"id": 9, "name": "Hurricane Force (64kts+)"},
                    {"id": 11, "name": "Observed Wind Swath"},
                ],
            }
        raise AssertionError(f"unexpected fetch_json url: {url}")

    def _fake_query_json(base_url: str, path: str, *, params, timeout_s: float, user_agent: str) -> dict[str, object]:
        if path == "1/query":
            return {"features": observed_features}
        if path == "0/query":
            return {"features": forecast_features}
        if path in {"7/query", "8/query", "9/query", "11/query"}:
            layer_id = int(path.split("/", 1)[0])
            feature = wind_feature if layer_id != 7 else wind_feature_b
            return {"features": [feature]}
        raise AssertionError(f"unexpected query path: {path}")

    monkeypatch.setattr(cyclone_client, "fetch_json", _fake_fetch_json)
    monkeypatch.setattr(cyclone_client, "query_json", _fake_query_json)

    collection = cyclone_client.fetch_active_hurricanes_snapshot()

    assert collection.service_name == "Active Hurricanes"
    assert len(collection.snapshots) == 2
    assert collection.snapshots[0].storm_name == "One-e"
    assert collection.snapshots[1].storm_name == "Jangmi"
    assert collection.summary_text().startswith("2 storms")


def test_fetch_latest_observed_feature_returns_none_when_query_is_empty(monkeypatch) -> None:
    def _fake_query_layer(*_args, **_kwargs) -> dict[str, object]:
        return {"features": []}

    monkeypatch.setattr(cyclone_client, "_query_layer", _fake_query_layer)

    assert cyclone_client.fetch_latest_observed_feature() is None


def test_tropical_cyclone_controller_treats_empty_observed_query_as_empty_overlay(
    monkeypatch,
    tmp_path: Path,
) -> None:
    controller = tropical_cyclone_controller.TropicalCycloneController(cache_root=tmp_path)
    ready_payloads: list[dict[str, object]] = []
    failed_payloads: list[str] = []
    saved_entries: list[TropicalCycloneCacheEntry] = []

    monkeypatch.setattr(
        tropical_cyclone_controller,
        "load_tropical_cyclone_cache",
        lambda _cache_root: None,
    )
    monkeypatch.setattr(
        tropical_cyclone_controller,
        "fetch_latest_observed_feature",
        lambda **_kwargs: None,
    )
    monkeypatch.setattr(
        tropical_cyclone_controller,
        "fetch_active_hurricanes_snapshot",
        lambda **_kwargs: pytest.fail("full snapshot fetch should be skipped when observed query is empty"),
    )
    monkeypatch.setattr(
        tropical_cyclone_controller,
        "save_tropical_cyclone_cache",
        lambda entry, *, cache_root: saved_entries.append(entry),
    )
    monkeypatch.setattr(
        controller,
        "_emit_ready",
        lambda payload, *, request_id: ready_payloads.append(payload),
    )
    monkeypatch.setattr(
        controller,
        "_emit_failed",
        lambda banner, *, request_id: failed_payloads.append(banner),
    )

    controller._run_update(reason="manual", request_id=1)  # noqa: SLF001

    assert failed_payloads == []
    assert len(ready_payloads) == 1
    assert ready_payloads[0]["snapshot_collection"]["snapshots"] == []
    assert ready_payloads[0]["banner"] == "Typhoon: none"
    assert len(saved_entries) == 1
    assert saved_entries[0].snapshot_collection.snapshots == ()


def test_tropical_cyclone_controller_logs_empty_observed_error_without_traceback(
    monkeypatch,
    tmp_path: Path,
    caplog,
) -> None:
    controller = tropical_cyclone_controller.TropicalCycloneController(cache_root=tmp_path)
    ready_payloads: list[dict[str, object]] = []
    saved_entries: list[TropicalCycloneCacheEntry] = []

    monkeypatch.setattr(
        tropical_cyclone_controller,
        "load_tropical_cyclone_cache",
        lambda _cache_root: None,
    )
    monkeypatch.setattr(
        tropical_cyclone_controller,
        "fetch_latest_observed_feature",
        lambda **_kwargs: {"attributes": {"STORMNAME": "Foo", "BASIN": "WP", "ADVDATE": 1}},
    )

    def _raise_empty_observed_error(**_kwargs):
        raise tropical_cyclone_controller.TropicalCycloneFetchError(
            "No observed position features returned"
        )

    monkeypatch.setattr(
        tropical_cyclone_controller,
        "fetch_active_hurricanes_snapshot",
        _raise_empty_observed_error,
    )
    monkeypatch.setattr(
        tropical_cyclone_controller,
        "save_tropical_cyclone_cache",
        lambda entry, *, cache_root: saved_entries.append(entry),
    )
    monkeypatch.setattr(
        controller,
        "_emit_ready",
        lambda payload, *, request_id: ready_payloads.append(payload),
    )

    with caplog.at_level(logging.WARNING, logger="zstarview.gui.tropical_cyclone_controller"):
        controller._run_update(reason="manual", request_id=1)  # noqa: SLF001

    assert ready_payloads
    assert ready_payloads[0]["banner"] == "Typhoon: none"
    assert len(saved_entries) == 1
    assert "No observed tropical cyclone positions returned; treating overlay as empty." in caplog.text
    assert "Traceback" not in caplog.text


def test_tropical_cyclone_controller_logs_generic_failure_without_traceback(
    monkeypatch,
    tmp_path: Path,
    caplog,
) -> None:
    controller = tropical_cyclone_controller.TropicalCycloneController(cache_root=tmp_path)
    failed_payloads: list[str] = []

    monkeypatch.setattr(
        tropical_cyclone_controller,
        "load_tropical_cyclone_cache",
        lambda _cache_root: None,
    )
    monkeypatch.setattr(
        tropical_cyclone_controller,
        "fetch_latest_observed_feature",
        lambda **_kwargs: {"attributes": {"STORMNAME": "Foo", "BASIN": "WP", "ADVDATE": 1}},
    )
    monkeypatch.setattr(
        tropical_cyclone_controller,
        "fetch_active_hurricanes_snapshot",
        lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("network unreachable")),
    )
    monkeypatch.setattr(
        controller,
        "_emit_failed",
        lambda banner, *, request_id: failed_payloads.append(banner),
    )

    with caplog.at_level(logging.WARNING, logger="zstarview.gui.tropical_cyclone_controller"):
        controller._run_update(reason="manual", request_id=1)  # noqa: SLF001

    assert failed_payloads == ["Typhoon: unavailable"]
    assert "Tropical cyclone update failed: network unreachable" in caplog.text
    assert "Traceback" not in caplog.text


def test_tropical_cyclone_draw_uses_far_marker_beyond_distance_limit(monkeypatch) -> None:
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
            distance_km=400.0001,
        )

    def _fake_far_marker(*_args, **_kwargs):
        calls.append(("far", None))
        return QPointF(12.0, 34.0)

    monkeypatch.setattr(render_tropical_cyclones, "_project_point_no_cutoff", _fake_center)
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

        def drawText(self, *_args, **_kwargs) -> None:
            pass

    render_tropical_cyclones.draw_tropical_cyclone_overlay(
        _FakePainter(),
        geometry=ScreenGeometry(center=(100, 100), radius=80),
        viewer=viewer,
        snapshot=snapshot,
        when_utc=datetime(2026, 5, 30, 2, 30, tzinfo=timezone.utc),
        theme=SimpleNamespace(),
        opacity=0.4,
        enabled=True,
    )

    assert calls == [("far", None)]


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


def test_tropical_cyclone_far_label_is_hidden_without_hover(monkeypatch) -> None:
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
    calls: list[str] = []

    def _fake_center(*_args, **_kwargs):
        return render_tropical_cyclones._RenderPoint(
            nx=0.1,
            ny=0.2,
            alt_deg=10.0,
            az_deg=20.0,
            distance_km=400.0001,
        )

    def _fake_far_marker(*_args, **_kwargs):
        calls.append("far")
        return QPointF(12.0, 34.0)

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

        def drawPolygon(self, *_args, **_kwargs) -> None:
            pass

        def drawText(self, *_args, **_kwargs) -> None:
            raise AssertionError("far label should be hidden when not hovered")

    monkeypatch.setattr(render_tropical_cyclones, "_project_point_no_cutoff", _fake_center)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_far_cyclone_marker", _fake_far_marker)

    render_tropical_cyclones.draw_tropical_cyclone_overlay(
        _FakePainter(),
        geometry=ScreenGeometry(center=(100, 100), radius=80),
        viewer=viewer,
        snapshot=snapshot,
        when_utc=datetime(2026, 5, 30, 2, 30, tzinfo=timezone.utc),
        theme=SimpleNamespace(),
        opacity=0.4,
        highlighted=False,
        enabled=True,
    )

    assert calls == ["far"]


def test_tropical_cyclone_far_label_is_drawn_when_hovered(monkeypatch) -> None:
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
    calls: list[str] = []

    def _fake_center(*_args, **_kwargs):
        return render_tropical_cyclones._RenderPoint(
            nx=0.1,
            ny=0.2,
            alt_deg=10.0,
            az_deg=20.0,
            distance_km=400.0001,
        )

    def _fake_far_marker(*_args, **_kwargs):
        calls.append("far")
        return QPointF(12.0, 34.0)

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

        def drawPolygon(self, *_args, **_kwargs) -> None:
            pass

        def drawText(self, pos, text) -> None:
            calls.append(f"text:{text}")
            assert isinstance(pos, QPointF)

    monkeypatch.setattr(render_tropical_cyclones, "_project_point_no_cutoff", _fake_center)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_far_cyclone_marker", _fake_far_marker)

    render_tropical_cyclones.draw_tropical_cyclone_overlay(
        _FakePainter(),
        geometry=ScreenGeometry(center=(100, 100), radius=80),
        viewer=viewer,
        snapshot=snapshot,
        when_utc=datetime(2026, 5, 30, 2, 30, tzinfo=timezone.utc),
        theme=SimpleNamespace(),
        opacity=0.4,
        highlighted=True,
        enabled=True,
    )

    assert calls == ["far", "text:Jangmi"]


def test_tropical_cyclone_far_marker_uses_thinner_pen_when_not_highlighted() -> None:
    point = render_tropical_cyclones._RenderPoint(
        nx=0.1,
        ny=0.2,
        alt_deg=10.0,
        az_deg=20.0,
        distance_km=400.0001,
    )
    pen_widths: list[float] = []

    class _FakePainter:
        def setPen(self, pen, *_args, **_kwargs) -> None:
            pen_widths.append(float(pen.widthF()))

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolygon(self, *_args, **_kwargs) -> None:
            pass

    render_tropical_cyclones._draw_far_cyclone_marker(
        _FakePainter(),
        point,
        geometry=ScreenGeometry(center=(100, 100), radius=80),
        color_rgba=(240, 122, 122, 200),
        highlighted=False,
    )

    assert pen_widths == [1.0]


def test_tropical_cyclone_far_marker_uses_wider_pen_when_highlighted() -> None:
    point = render_tropical_cyclones._RenderPoint(
        nx=0.1,
        ny=0.2,
        alt_deg=10.0,
        az_deg=20.0,
        distance_km=400.0001,
    )
    pen_widths: list[float] = []

    class _FakePainter:
        def setPen(self, pen, *_args, **_kwargs) -> None:
            pen_widths.append(float(pen.widthF()))

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolygon(self, *_args, **_kwargs) -> None:
            pass

    render_tropical_cyclones._draw_far_cyclone_marker(
        _FakePainter(),
        point,
        geometry=ScreenGeometry(center=(100, 100), radius=80),
        color_rgba=(240, 122, 122, 200),
        highlighted=True,
    )

    assert pen_widths == [2.0]


def test_find_highlighted_tropical_cyclone_prefers_nearest_far_marker(monkeypatch) -> None:
    snapshots = (
        TropicalCycloneSnapshot(
            storm_name="Near",
            basin="WP",
            advdate_utc=datetime(2026, 5, 30, 2, 10, tzinfo=timezone.utc),
            observed_position=TropicalCyclonePoint(lat_deg=10.0, lon_deg=20.0),
            wind_polygons=(),
        ),
        TropicalCycloneSnapshot(
            storm_name="Far",
            basin="WP",
            advdate_utc=datetime(2026, 5, 30, 2, 10, tzinfo=timezone.utc),
            observed_position=TropicalCyclonePoint(lat_deg=30.0, lon_deg=40.0),
            wind_polygons=(),
        ),
    )

    def _identity_project(snapshot, _when_utc):
        return snapshot

    def _fake_center(lat_deg, lon_deg, *, viewer, height_m):
        del viewer, height_m
        if lat_deg == 10.0 and lon_deg == 20.0:
            return render_tropical_cyclones._RenderPoint(
                nx=10.0,
                ny=10.0,
                alt_deg=10.0,
                az_deg=20.0,
                distance_km=450.0,
            )
        return render_tropical_cyclones._RenderPoint(
            nx=40.0,
            ny=40.0,
            alt_deg=10.0,
            az_deg=20.0,
            distance_km=450.0,
        )

    def _identity_screen(nx, ny, geometry):
        del geometry
        return float(nx), float(ny)

    monkeypatch.setattr(render_tropical_cyclones, "project_tropical_cyclone_snapshot", _identity_project)
    monkeypatch.setattr(render_tropical_cyclones, "_project_point_no_cutoff", _fake_center)
    monkeypatch.setattr(render_tropical_cyclones, "normalized_to_screen_xy", _identity_screen)

    highlighted = render_tropical_cyclones.find_highlighted_tropical_cyclone(
        snapshots,
        QPointF(12.0, 12.0),
        ScreenGeometry(center=(100, 100), radius=80),
        viewer_data=SimpleNamespace(),
        time_obj=datetime(2026, 5, 30, 2, 30, tzinfo=timezone.utc),
    )

    assert highlighted is not None
    assert highlighted[0].storm_name == "Near"


def test_tropical_cyclone_tether_matches_asterism_hover_passes() -> None:
    tether_points = (
        render_tropical_cyclones._RenderPoint(
            nx=0.1,
            ny=0.2,
            alt_deg=10.0,
            az_deg=20.0,
            distance_km=120.0,
        ),
        render_tropical_cyclones._RenderPoint(
            nx=0.12,
            ny=0.23,
            alt_deg=11.0,
            az_deg=20.0,
            distance_km=120.0,
        ),
        render_tropical_cyclones._RenderPoint(
            nx=0.14,
            ny=0.26,
            alt_deg=12.0,
            az_deg=20.0,
            distance_km=120.0,
        ),
        render_tropical_cyclones._RenderPoint(
            nx=0.16,
            ny=0.29,
            alt_deg=13.0,
            az_deg=20.0,
            distance_km=120.0,
        ),
        render_tropical_cyclones._RenderPoint(
            nx=0.18,
            ny=0.32,
            alt_deg=14.0,
            az_deg=20.0,
            distance_km=120.0,
        ),
        render_tropical_cyclones._RenderPoint(
            nx=0.2,
            ny=0.35,
            alt_deg=15.0,
            az_deg=20.0,
            distance_km=120.0,
        ),
    )
    pen_widths: list[float] = []
    pen_alphas: list[int] = []
    line_calls: list[tuple[QPointF, QPointF]] = []

    class _FakePainter:
        def setPen(self, pen, *_args, **_kwargs) -> None:
            pen_widths.append(float(pen.widthF()))
            pen_alphas.append(int(pen.color().alpha()))

        def drawLine(self, start, end) -> None:
            line_calls.append((start, end))

    render_tropical_cyclones._draw_marker_tether(
        _FakePainter(),
        tether_points,
        geometry=ScreenGeometry(center=(100, 100), radius=80),
        color_rgba=(240, 122, 122, 200),
    )

    assert pen_widths == [1.0]
    assert pen_alphas == [204]
    assert len(line_calls) == 5


def test_tropical_cyclone_draws_filled_marker_at_5km_when_in_range(monkeypatch) -> None:
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
    filled_colors: list[tuple[int, int, int, int]] = []

    def _fake_project(*_args, **kwargs):
        height_m = float(kwargs["height_m"])
        if height_m == 0.0:
            return render_tropical_cyclones._RenderPoint(
                nx=0.1,
                ny=0.2,
                alt_deg=10.0,
                az_deg=20.0,
                distance_km=399.0,
            )
        return render_tropical_cyclones._RenderPoint(
            nx=0.2,
            ny=0.35,
            alt_deg=15.0,
            az_deg=20.0,
            distance_km=399.0,
        )

    def _fake_tether(_painter, tether_points, **_kwargs):
        calls.append(("tether", tuple(tether_points)))

    def _fake_far_marker(*_args, **_kwargs):
        calls.append(("far", None))
        return QPointF(12.0, 34.0)

    def _fake_filled_marker(_painter, point, **_kwargs):
        calls.append(("filled", point))
        filled_colors.append(tuple(_kwargs["color_rgba"]))

    monkeypatch.setattr(render_tropical_cyclones, "_project_point_no_cutoff", _fake_project)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_marker_tether", _fake_tether)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_far_cyclone_marker", _fake_far_marker)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_filled_cyclone_marker", _fake_filled_marker)

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

    assert [kind for kind, _ in calls] == ["tether", "filled"]
    tether_points = calls[0][1]
    ground_point = tether_points[0]
    marker_point = tether_points[-1]
    assert ground_point.alt_deg == pytest.approx(10.0)
    assert marker_point.alt_deg == pytest.approx(15.0)
    assert calls[1][1] == marker_point
    assert not any(kind == "far" for kind, _ in calls)
    assert filled_colors == [
        (
            render_tropical_cyclones.TROPICAL_CYCLONE_COLOR_RGB[0],
            render_tropical_cyclones.TROPICAL_CYCLONE_COLOR_RGB[1],
            render_tropical_cyclones.TROPICAL_CYCLONE_COLOR_RGB[2],
            round(255 * render_tropical_cyclones.TROPICAL_CYCLONE_MARKER_ALPHA_SCALE),
        )
    ]


def test_tropical_cyclone_in_range_marker_projects_5km_height(monkeypatch) -> None:
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
    requested_heights: list[float] = []

    def _fake_project(*_args, **kwargs):
        height_m = float(kwargs["height_m"])
        requested_heights.append(height_m)
        return render_tropical_cyclones._RenderPoint(
            nx=0.1,
            ny=0.2,
            alt_deg=10.0,
            az_deg=20.0,
            distance_km=120.0,
        )

    monkeypatch.setattr(render_tropical_cyclones, "_project_point_no_cutoff", _fake_project)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_marker_tether", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        render_tropical_cyclones,
        "_draw_filled_cyclone_marker",
        lambda *_args, **_kwargs: None,
    )

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

    render_tropical_cyclones.draw_tropical_cyclone_overlay(
        _FakePainter(),
        geometry=ScreenGeometry(center=(100, 100), radius=80),
        viewer=viewer,
        snapshot=snapshot,
        when_utc=datetime(2026, 5, 30, 2, 30, tzinfo=timezone.utc),
        theme=SimpleNamespace(),
        opacity=0.4,
        enabled=True,
    )

    assert requested_heights == [
        0.0,
        render_tropical_cyclones.TROPICAL_CYCLONE_MARKER_HEIGHT_M,
        1000.0,
        2000.0,
        3000.0,
        4000.0,
    ]


def test_tropical_cyclone_draws_filled_marker_and_tether_when_in_range(monkeypatch) -> None:
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

    def _fake_project(*_args, **kwargs):
        height_m = float(kwargs["height_m"])
        if height_m == 0.0:
            return render_tropical_cyclones._RenderPoint(
                nx=0.1,
                ny=0.2,
                alt_deg=10.0,
                az_deg=20.0,
                distance_km=100.0,
            )
        return render_tropical_cyclones._RenderPoint(
            nx=0.2,
            ny=0.35,
            alt_deg=15.0,
            az_deg=20.0,
            distance_km=100.0,
        )

    def _fake_filled_marker(*_args, **_kwargs):
        calls.append(("filled", None))

    def _fake_far_marker(*_args, **_kwargs):
        calls.append(("far", None))
        return QPointF(12.0, 34.0)

    def _fake_tether(*_args, **_kwargs):
        calls.append(("tether", None))

    monkeypatch.setattr(render_tropical_cyclones, "_project_point_no_cutoff", _fake_project)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_filled_cyclone_marker", _fake_filled_marker)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_far_cyclone_marker", _fake_far_marker)
    monkeypatch.setattr(render_tropical_cyclones, "_draw_marker_tether", _fake_tether)

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

    assert [kind for kind, _ in calls] == ["tether", "filled"]
    assert not any(kind == "far" for kind, _ in calls)
