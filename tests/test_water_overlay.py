from __future__ import annotations

from zstarview.water_overlay import (
    WaterPolygonFootprint,
    WaterSurfacePatch,
    assemble_rings_from_segments,
    WaterOverlayPoint,
    classify_water_surface_mode,
    build_geometric_distance_samples,
    extract_water_polygons,
    sample_water_overlay_points,
)
from zstarview.render.terrain import _thin_water_overlay_points_pairwise


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


def test_sample_water_overlay_points_uses_fallback_surface_height() -> None:
    footprint = WaterPolygonFootprint(
        water_id="lake",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (-0.01, -0.01),
                (0.01, -0.01),
                (0.01, 0.01),
                (-0.01, 0.01),
                (-0.01, -0.01),
            ),
        ),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water"},
    )

    sea_level_points = sample_water_overlay_points(
        (footprint,),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=100.0,
        fallback_surface_height_m=0.0,
        max_distance_km=0.2,
        sample_step_m=100.0,
        azimuth_step_deg=90.0,
    )
    local_surface_points = sample_water_overlay_points(
        (footprint,),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=100.0,
        fallback_surface_height_m=100.0,
        max_distance_km=0.2,
        sample_step_m=100.0,
        azimuth_step_deg=90.0,
    )

    assert sea_level_points
    assert local_surface_points
    assert max(point.alt_deg for point in local_surface_points) > -0.01
    assert max(point.alt_deg for point in sea_level_points) < -20.0


def test_sample_water_overlay_points_keeps_coastline_at_sea_level() -> None:
    footprint = WaterPolygonFootprint(
        water_id="coast",
        kind="coastline",
        outer_rings_lonlat=(
            (
                (-0.01, -0.01),
                (0.01, -0.01),
                (0.01, 0.01),
                (-0.01, 0.01),
                (-0.01, -0.01),
            ),
        ),
        inner_rings_lonlat=(),
        source="coastline",
        tags={"natural": "coastline"},
    )

    points = sample_water_overlay_points(
        (footprint,),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=100.0,
        fallback_surface_height_m=100.0,
        max_distance_km=0.2,
        sample_step_m=100.0,
        azimuth_step_deg=90.0,
    )

    assert points
    assert max(point.alt_deg for point in points) < -20.0


def test_build_geometric_distance_samples_stays_dense_farther_out() -> None:
    samples = build_geometric_distance_samples(2.0, 1.25**5)

    assert len(samples) > 40
    assert samples[-1] <= 2000.0
    assert (samples[-1] - samples[-2]) < 300.0


def test_thin_water_overlay_points_pairwise_keeps_one_point_per_pair() -> None:
    points = [
        WaterOverlayPoint("water", 1.0, 10.0, 0.1, scan_azimuth_index=0, scan_distance_index=0),
        WaterOverlayPoint("water", 2.0, 20.0, 0.2, scan_azimuth_index=0, scan_distance_index=1),
        WaterOverlayPoint("water", 3.0, 30.0, 0.3, scan_azimuth_index=0, scan_distance_index=2),
        WaterOverlayPoint("water", 4.0, 40.0, 0.4, scan_azimuth_index=0, scan_distance_index=3),
        WaterOverlayPoint("water", 6.0, 60.0, 0.6, scan_azimuth_index=1, scan_distance_index=0),
        WaterOverlayPoint("water", 7.0, 70.0, 0.7, scan_azimuth_index=1, scan_distance_index=1),
        WaterOverlayPoint("water", 8.0, 80.0, 0.8, scan_azimuth_index=1, scan_distance_index=2),
        WaterOverlayPoint("water", 9.0, 90.0, 0.9, scan_azimuth_index=1, scan_distance_index=3),
        WaterOverlayPoint("water", 5.0, 50.0, 0.5),
    ]

    got = _thin_water_overlay_points_pairwise(points)

    assert [point.scan_distance_index for point in got if point.scan_azimuth_index == 0] == [0, 2]
    assert [point.scan_distance_index for point in got if point.scan_azimuth_index == 1] == [1, 3]
    assert any(point.scan_distance_index is None for point in got)
