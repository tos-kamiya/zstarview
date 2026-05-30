from __future__ import annotations

from datetime import datetime, timezone

from zstarview.tropical_cyclones.cache import (
    TropicalCycloneCacheEntry,
    load_tropical_cyclone_cache,
    save_tropical_cyclone_cache,
)
from zstarview.tropical_cyclones.models import (
    TropicalCyclonePoint,
    TropicalCyclonePolygon,
    TropicalCycloneSnapshot,
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
