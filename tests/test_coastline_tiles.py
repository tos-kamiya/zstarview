from __future__ import annotations

import gzip
import json
from pathlib import Path

from zstarview.coastline_tiles import (
    PREVIEW_ROOT_ENV,
    _tile_roots,
    _clip_line_to_radius,
    _line_coordinates,
    load_coastline_overlay_polylines,
)
from zstarview.water_overlay import (
    DEFAULT_WATER_BOUNDARY_RADIUS_KM,
    DEFAULT_WATER_RADIUS_KM,
    resolve_water_scan_radius_km,
)


def _write_tile(root: Path, name: str, geometry: dict[str, object], *, gzip_output: bool = False) -> None:
    payload = {
        "type": "FeatureCollection",
        "features": [{"type": "Feature", "properties": {}, "geometry": geometry}],
    }
    path = root / f"{name}.geojson"
    if gzip_output:
        with gzip.open(f"{path}.gz", "wt", encoding="utf-8") as handle:
            json.dump(payload, handle)
    else:
        path.write_text(json.dumps(payload), encoding="utf-8")


def test_line_coordinates_supports_lines_and_geometry_collections() -> None:
    geometry = {
        "type": "GeometryCollection",
        "geometries": [
            {
                "type": "LineString",
                "coordinates": [[1, 2], [3, 4]],
            },
            {
                "type": "MultiLineString",
                "coordinates": [[[-1, -2], [-3, -4]]],
            },
        ],
    }

    assert tuple(_line_coordinates(geometry)) == (
        ((1.0, 2.0), (3.0, 4.0)),
        ((-1.0, -2.0), (-3.0, -4.0)),
    )


def test_clip_line_to_radius_interpolates_boundary() -> None:
    fragments = _clip_line_to_radius(((0.0, 0.0), (20_000.0, 0.0)), 10_000.0)

    assert len(fragments) == 1
    assert fragments[0][0] == (0.0, 0.0)
    assert fragments[0][-1] == (10_000.0, 0.0)


def test_loader_reads_normal_tile_and_clips_to_ten_km(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setenv(PREVIEW_ROOT_ENV, str(tmp_path))
    _write_tile(
        tmp_path,
        "tile_y7_x16",
        {
            "type": "LineString",
            "coordinates": [[0.01, 0.0], [0.2, 0.0]],
        },
        gzip_output=True,
    )

    polylines = load_coastline_overlay_polylines(
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=10.0,
    )

    assert polylines
    assert all(point.distance_km <= 10.001 for line in polylines for point in line.points)
    assert all(line.water_category == "coastline" for line in polylines)


def test_loader_reads_downloaded_nested_tile(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv(PREVIEW_ROOT_ENV, str(tmp_path))
    tile_dir = tmp_path / "y07" / "x16"
    tile_dir.mkdir(parents=True)
    payload = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {},
                "geometry": {"type": "LineString", "coordinates": [[0.01, 0.0], [0.02, 0.0]]},
            }
        ],
    }
    with gzip.open(tile_dir / "tile.geojson.gz", "wt", encoding="utf-8") as handle:
        json.dump(payload, handle)

    polylines = load_coastline_overlay_polylines(
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=10.0,
    )

    assert polylines


def test_cache_root_requires_ready_marker(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.delenv(PREVIEW_ROOT_ENV, raising=False)
    monkeypatch.setattr("zstarview.coastline_tiles.CACHE_PATH", str(tmp_path))
    cache_root = tmp_path / "coastline" / "osm-water-polygons" / "20260725" / "schema-1"
    grid_root = cache_root / "grid-32x16"
    grid_root.mkdir(parents=True)
    assert grid_root not in _tile_roots()
    (cache_root / "READY").write_text("ready\n", encoding="ascii")
    assert grid_root in _tile_roots()


def test_loader_prefers_split_children_over_parent(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv(PREVIEW_ROOT_ENV, str(tmp_path))
    (tmp_path / "tile_y7_x16.geojson").write_text("not valid json", encoding="utf-8")
    _write_tile(
        tmp_path,
        "tile_y7_x16_q00",
        {
            "type": "LineString",
            "coordinates": [[0.01, 0.0], [0.02, 0.0]],
        },
    )

    polylines = load_coastline_overlay_polylines(
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=10.0,
    )

    assert polylines
    assert all("q00" in line.water_id for line in polylines)


def test_loader_treats_zero_marker_as_empty(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv(PREVIEW_ROOT_ENV, str(tmp_path))
    (tmp_path / "tile_y7_x16.0").write_bytes(b"")

    assert load_coastline_overlay_polylines(
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=10.0,
    ) == ()


def test_loader_ignores_invalid_geojson(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv(PREVIEW_ROOT_ENV, str(tmp_path))
    (tmp_path / "tile_y7_x16.geojson").write_text("not valid json", encoding="utf-8")

    assert load_coastline_overlay_polylines(
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=10.0,
    ) == ()


def test_dot_radius_remains_separate_from_boundary_radius() -> None:
    assert DEFAULT_WATER_RADIUS_KM == 2.0
    assert DEFAULT_WATER_BOUNDARY_RADIUS_KM == 10.0
    assert resolve_water_scan_radius_km(1_000.0) > DEFAULT_WATER_BOUNDARY_RADIUS_KM
