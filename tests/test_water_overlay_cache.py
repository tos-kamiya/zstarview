from __future__ import annotations

import json

from zstarview.gui.water_overlay_cache import (
    WATER_OVERLAY_CACHE_FORMAT_VERSION,
    WaterOverlayCacheSnapshot,
    load_water_overlay_cache,
    save_water_overlay_cache,
)
from zstarview.water_overlay import WaterPolygonFootprint


def _footprint() -> WaterPolygonFootprint:
    return WaterPolygonFootprint(
        water_id="way/1",
        kind="natural_water",
        outer_rings_lonlat=(((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)),),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water"},
    )


def test_load_water_overlay_cache_rejects_old_format(tmp_path) -> None:
    path = tmp_path / "earth_+00.0000_+000.0000_r2.00.json"
    path.write_text(
        json.dumps(
            {
                "cache_format_version": WATER_OVERLAY_CACHE_FORMAT_VERSION - 1,
                "scope_key": "earth_+00.0000_+000.0000_r2.00",
                "fetched_at_utc": "2026-05-13T00:00:00+00:00",
                "water_polygon_count": 1,
                "footprints": [],
            }
        ),
        encoding="utf-8",
    )

    assert load_water_overlay_cache("earth_+00.0000_+000.0000_r2.00", cache_root=tmp_path) is None


def test_save_and_load_water_overlay_cache_roundtrip(tmp_path) -> None:
    snapshot = WaterOverlayCacheSnapshot(
        footprints=(_footprint(),),
        water_polygon_count=1,
        fetched_at_utc=None,
    )
    saved_path = save_water_overlay_cache("earth_+00.0000_+000.0000_r2.00", snapshot, cache_root=tmp_path)
    loaded = load_water_overlay_cache("earth_+00.0000_+000.0000_r2.00", cache_root=tmp_path)

    assert saved_path.exists()
    assert loaded is not None
    assert loaded.water_polygon_count == 1
    assert len(loaded.footprints) == 1
