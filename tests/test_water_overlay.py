from __future__ import annotations

from zstarview.water_overlay import (
    WaterSurfacePatch,
    assemble_rings_from_segments,
    classify_water_surface_mode,
    extract_water_polygons,
)


def _node(node_id: int, lon: float, lat: float) -> dict[str, object]:
    return {
        "type": "node",
        "id": node_id,
        "lon": lon,
        "lat": lat,
    }


def _way(way_id: int, node_ids: list[int], tags: dict[str, str]) -> dict[str, object]:
    return {
        "type": "way",
        "id": way_id,
        "nodes": node_ids,
        "tags": tags,
    }


def test_assemble_rings_from_segments_reconstructs_closed_ring() -> None:
    ring = assemble_rings_from_segments(
        [
            ((0.0, 0.0), (1.0, 0.0)),
            ((1.0, 0.0), (1.0, 1.0)),
            ((1.0, 1.0), (0.0, 1.0), (0.0, 0.0)),
        ]
    )
    assert ring == [((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0))]


def test_extract_water_polygons_keeps_outer_and_inner_rings() -> None:
    elements: list[dict[str, object]] = [
        _node(1, 0.0, 0.0),
        _node(2, 4.0, 0.0),
        _node(3, 4.0, 4.0),
        _node(4, 0.0, 4.0),
        _node(5, 1.0, 1.0),
        _node(6, 2.0, 1.0),
        _node(7, 2.0, 2.0),
        _node(8, 1.0, 2.0),
        _way(10, [1, 2, 3, 4, 1], {"natural": "water"}),
        _way(20, [5, 6, 7, 8, 5], {}),
        {
            "type": "relation",
            "id": 100,
            "tags": {"natural": "water", "type": "multipolygon"},
            "members": [
                {"type": "way", "ref": 10, "role": "outer"},
                {"type": "way", "ref": 20, "role": "inner"},
            ],
        },
    ]

    polygons = extract_water_polygons(elements)
    assert len(polygons) == 2
    assert polygons[0].water_id == "way/10"
    assert polygons[0].outer_rings_lonlat == (((0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)),)
    assert polygons[0].inner_rings_lonlat == ()
    assert polygons[1].water_id == "relation/100"
    assert len(polygons[1].outer_rings_lonlat) == 1
    assert len(polygons[1].inner_rings_lonlat) == 1


def test_extract_water_polygons_builds_coastline_water_with_island_hole() -> None:
    elements: list[dict[str, object]] = [
        _node(1, 3.0, 3.0),
        _node(2, 7.0, 3.0),
        _node(3, 7.0, 7.0),
        _node(4, 3.0, 7.0),
        _node(5, 1.0, 1.0),
        _node(6, 2.0, 1.0),
        _node(7, 2.0, 2.0),
        _node(8, 1.0, 2.0),
        _way(10, [1, 2, 3, 4, 1], {"natural": "coastline"}),
        _way(20, [5, 6, 7, 8, 5], {"natural": "coastline"}),
    ]

    polygons = extract_water_polygons(
        elements,
        bbox=(0.0, 0.0, 10.0, 10.0),
    )
    coastline_polygons = [polygon for polygon in polygons if polygon.source == "coastline"]

    assert len(coastline_polygons) == 1
    assert coastline_polygons[0].kind == "coastline"
    assert len(coastline_polygons[0].inner_rings_lonlat) == 2


def test_water_surface_patch_classifies_flat_and_sloped() -> None:
    flat_patch = WaterSurfacePatch(
        patch_id="flat",
        polygon_id="relation/1",
        anchor_points_lonlat=((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
        anchor_elevations_m=(100.0, 100.4, 100.9),
        flat_threshold_m=1.0,
    )
    sloped_patch = WaterSurfacePatch(
        patch_id="sloped",
        polygon_id="relation/2",
        anchor_points_lonlat=((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
        anchor_elevations_m=(100.0, 100.5, 102.2),
        flat_threshold_m=1.0,
    )

    assert flat_patch.surface_mode == "flat"
    assert sloped_patch.surface_mode == "sloped"
    assert classify_water_surface_mode((100.0, 101.0, 102.0), flat_threshold_m=3.0) == "flat"
