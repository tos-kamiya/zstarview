from __future__ import annotations

import math
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


def _load_probe_module():
    module_path = Path(__file__).resolve().parents[1] / "dev-samples" / "overpass_water_overlay_probe.py"
    spec = spec_from_file_location("overpass_water_overlay_probe", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("failed to load overpass_water_overlay_probe.py")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_expanded_query_bbox_from_point_scales_radius_by_20_percent() -> None:
    probe = _load_probe_module()
    view_bbox = probe.bbox_from_point(34.6825, 135.1867, 5.0)
    query_bbox = probe.expanded_query_bbox_from_point(34.6825, 135.1867, 5.0)

    view_width = view_bbox[2] - view_bbox[0]
    view_height = view_bbox[3] - view_bbox[1]
    query_width = query_bbox[2] - query_bbox[0]
    query_height = query_bbox[3] - query_bbox[1]

    assert math.isclose(query_width, view_width * 1.2, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(query_height, view_height * 1.2, rel_tol=0.0, abs_tol=1e-12)


def test_expanded_query_bbox_from_point_keeps_center() -> None:
    probe = _load_probe_module()
    view_bbox = probe.bbox_from_point(34.6825, 135.1867, 5.0)
    query_bbox = probe.expanded_query_bbox_from_point(34.6825, 135.1867, 5.0)

    view_center_lon = (view_bbox[0] + view_bbox[2]) * 0.5
    view_center_lat = (view_bbox[1] + view_bbox[3]) * 0.5
    query_center_lon = (query_bbox[0] + query_bbox[2]) * 0.5
    query_center_lat = (query_bbox[1] + query_bbox[3]) * 0.5

    assert math.isclose(query_center_lon, view_center_lon, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(query_center_lat, view_center_lat, rel_tol=0.0, abs_tol=1e-12)


def test_polygon_overlaps_bbox_rejects_disjoint_polygons() -> None:
    probe = _load_probe_module()
    polygon = probe.WaterPolygon(
        osm_id="coastline/0",
        kind="coastline",
        outer_rings=(((20.0, 20.0), (21.0, 20.0), (21.0, 21.0), (20.0, 21.0), (20.0, 20.0)),),
        inner_rings=(),
        source="coastline",
        tags={"natural": "coastline"},
    )

    assert not probe._polygon_overlaps_bbox(polygon, (0.0, 0.0, 10.0, 10.0))


def test_filter_polygon_by_size_drops_small_outer_ring_for_other() -> None:
    probe = _load_probe_module()
    polygon = probe.WaterPolygon(
        osm_id="relation/1",
        kind="natural_water",
        outer_rings=(
            ((0.0500, 0.0500), (0.1500, 0.0500), (0.1500, 0.1500), (0.0500, 0.1500), (0.0500, 0.0500)),
            ((0.0500, 0.0500), (0.05005, 0.0500), (0.05005, 0.05005), (0.0500, 0.05005), (0.0500, 0.0500)),
        ),
        inner_rings=(),
        source="relation",
        tags={"natural": "water"},
    )

    filtered, removed_rings = probe._filter_polygon_by_size(
        polygon,
        center_lon_deg=0.0,
        center_lat_deg=0.0,
    )

    assert filtered is not None
    assert len(filtered.outer_rings) == 1
    assert removed_rings == 1


def test_load_request_bbox_from_cache_reads_query_bbox(tmp_path) -> None:
    probe = _load_probe_module()
    cache_path = tmp_path / "cache.json"
    cache_path.write_text(
        """{
  "cache_format_version": 2,
  "query": {
    "bbox": {
      "west": 1.0,
      "south": 2.0,
      "east": 3.0,
      "north": 4.0
    }
  },
  "footprints": []
}
""",
        encoding="utf-8",
    )

    assert probe.load_request_bbox_from_cache(cache_path) == (1.0, 2.0, 3.0, 4.0)
