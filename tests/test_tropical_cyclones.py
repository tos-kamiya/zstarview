from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import pytest

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
        return [SimpleNamespace(alt_deg=10.0, az_deg=20.0, distance_km=128.0001)]

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
