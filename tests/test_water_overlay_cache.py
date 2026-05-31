from __future__ import annotations

import json
import os

from zstarview.gui.water_overlay_cache import (
    WATER_OVERLAY_CACHE_FORMAT_VERSION,
    WaterOverlayCacheSnapshot,
    load_water_overlay_cache,
    save_water_overlay_cache,
    water_overlay_cache_legacy_path,
    water_overlay_cache_scope_key,
    water_overlay_cache_path,
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
    path = water_overlay_cache_path("earth_+00.0000_+000.0000_r2.00", cache_root=tmp_path)
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
    assert saved_path.name.endswith("_simplified.json")
    assert loaded is not None
    assert loaded.water_polygon_count == 1
    assert len(loaded.footprints) == 1


def test_water_overlay_cache_scope_key_adds_density_suffix_for_non_default_steps() -> None:
    base = water_overlay_cache_scope_key(
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        radius_km=2.0,
    )
    dense = water_overlay_cache_scope_key(
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        radius_km=2.0,
        azimuth_step_deg=1.0,
    )

    assert base == "earth_+00.0000_+000.0000_r2.00"
    assert dense == "earth_+00.0000_+000.0000_r2.00_a1.00"


def test_load_water_overlay_cache_migrates_newer_legacy_cache(tmp_path, caplog) -> None:
    scope_key = "earth_+00.0000_+000.0000_r2.00"
    simplified_path = water_overlay_cache_path(scope_key, cache_root=tmp_path)
    legacy_path = water_overlay_cache_legacy_path(scope_key, cache_root=tmp_path)

    save_water_overlay_cache(scope_key, WaterOverlayCacheSnapshot(footprints=(), water_polygon_count=0), cache_root=tmp_path)
    legacy_path.write_text(
        json.dumps(
            {
                "cache_format_version": WATER_OVERLAY_CACHE_FORMAT_VERSION,
                "scope_key": scope_key,
                "fetched_at_utc": None,
                "water_polygon_count": 1,
                "footprints": [
                    {
                        "water_id": "river",
                        "kind": "natural_water",
                        "outer_rings_lonlat": [
                            [
                                [0.0200, 0.0000],
                                [0.0201, 0.0000],
                                [0.0202, 0.0000],
                                [0.0203, 0.0000],
                                [0.0210, 0.0000],
                                [0.0210, 0.0010],
                                [0.0200, 0.0010],
                                [0.0200, 0.0000],
                            ]
                        ],
                        "inner_rings_lonlat": [],
                        "source": "way",
                        "tags": {"natural": "water", "water": "river"},
                    }
                ],
            },
            separators=(",", ":"),
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    now = 1_000_000_000
    os.utime(simplified_path, (now, now))
    os.utime(legacy_path, (now + 10, now + 10))

    with caplog.at_level("INFO", logger="zstarview.gui.water_overlay_cache"):
        loaded = load_water_overlay_cache(
            scope_key,
            cache_root=tmp_path,
            observer_lat_deg=0.0,
            observer_lon_deg=0.0,
        )

    assert loaded is not None
    assert len(loaded.footprints[0].outer_rings_lonlat[0]) < 8
    assert simplified_path.exists()
    assert simplified_path.stat().st_mtime_ns >= legacy_path.stat().st_mtime_ns
    assert "Water overlay cache simplification started" in caplog.text
    assert "Water overlay cache simplification finished" in caplog.text
