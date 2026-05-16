from __future__ import annotations

import math
import urllib.error
import threading
from unittest.mock import Mock

import pytest

from zstarview.render.terrain import _thin_water_overlay_points_pairwise
from zstarview.render.terrain import _water_overlay_distance_alpha_scale
from zstarview.render.terrain import _water_overlay_point_color_rgb
from zstarview.water_overlay import (
    WaterOverlayPoint,
    WaterPolygonFootprint,
    WaterSurfacePatch,
    assemble_rings_from_segments,
    build_geometric_distance_samples,
    build_overpass_query,
    classify_water_surface_category,
    classify_water_surface_mode,
    extract_water_polygons,
    horizon_distance_km_from_height,
    expanded_query_bbox_from_point,
    resolve_water_scan_radius_km,
    sample_water_overlay_points,
    simplify_water_footprints_for_observer,
)
from zstarview.clouddisc.types import DownloadCancelledError


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


def test_classify_water_surface_category_uses_tags_and_kind() -> None:
    assert classify_water_surface_category({"natural": "coastline"}) == "sea"
    assert classify_water_surface_category({"water": "lake"}) == "lake"
    assert classify_water_surface_category({"water": "river"}) == "river"
    assert classify_water_surface_category({"waterway": "riverbank"}) == "river"


def test_water_overlay_point_color_rgb_distinguishes_sea_and_inland_water() -> None:
    sea_color = _water_overlay_point_color_rgb(
        WaterOverlayPoint("sea", 0.0, 0.0, 0.0, water_category="sea")
    )
    river_color = _water_overlay_point_color_rgb(
        WaterOverlayPoint("river", 0.0, 0.0, 0.0, water_category="river")
    )
    lake_color = _water_overlay_point_color_rgb(
        WaterOverlayPoint("lake", 0.0, 0.0, 0.0, water_category="lake")
    )

    assert sea_color != river_color
    assert sea_color != lake_color
    assert river_color != lake_color


def test_water_overlay_distance_alpha_scale_decays_with_distance() -> None:
    assert _water_overlay_distance_alpha_scale(0.0) == 1.0
    assert _water_overlay_distance_alpha_scale(128.0) == pytest.approx(0.0625, rel=1e-6)
    assert _water_overlay_distance_alpha_scale(256.0) == pytest.approx(0.00390625, rel=1e-6)


def test_water_surface_height_selection_prefers_explicit_level() -> None:
    from zstarview import water_overlay

    sea = WaterPolygonFootprint(
        water_id="sea",
        kind="coastline",
        outer_rings_lonlat=(((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)),),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "coastline", "water_level": "12.5"},
    )
    river = WaterPolygonFootprint(
        water_id="river",
        kind="natural_water",
        outer_rings_lonlat=(((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)),),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water", "water": "river", "ele": "7.25"},
    )

    assert (
        water_overlay._water_surface_height_m(  # noqa: SLF001
            sea,
            fallback_surface_height_m=3.0,
            latitude_deg=0.0,
            longitude_deg=0.0,
        )
        == 12.5
    )
    assert (
        water_overlay._water_surface_height_m(  # noqa: SLF001
            river,
            fallback_surface_height_m=3.0,
            target_ground_elevation_m_sampler=lambda *_args: 99.0,
            latitude_deg=0.0,
            longitude_deg=0.0,
        )
        == 7.25
    )


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


def test_simplify_water_footprints_for_observer_thins_dense_far_ring() -> None:
    footprint = WaterPolygonFootprint(
        water_id="river",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (0.0200, 0.0000),
                (0.0201, 0.0000),
                (0.0202, 0.0000),
                (0.0203, 0.0000),
                (0.0210, 0.0000),
                (0.0210, 0.0010),
                (0.0200, 0.0010),
                (0.0200, 0.0000),
            ),
        ),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water", "water": "river"},
    )
    reversed_footprint = WaterPolygonFootprint(
        water_id="river-reversed",
        kind="natural_water",
        outer_rings_lonlat=(tuple(reversed(footprint.outer_rings_lonlat[0])),),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water", "water": "river"},
    )

    simplified = simplify_water_footprints_for_observer(
        (footprint,),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
    )
    simplified_reversed = simplify_water_footprints_for_observer(
        (reversed_footprint,),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
    )

    assert len(simplified) == 1
    assert len(simplified[0].outer_rings_lonlat[0]) < len(footprint.outer_rings_lonlat[0])
    assert simplified[0].outer_rings_lonlat[0][0] == simplified[0].outer_rings_lonlat[0][-1]
    assert len(simplified_reversed[0].outer_rings_lonlat[0]) == len(simplified[0].outer_rings_lonlat[0])


def test_simplify_water_footprints_for_observer_drops_small_far_water() -> None:
    far_small_lake = WaterPolygonFootprint(
        water_id="far-lake",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (0.2000, 0.2000),
                (0.20005, 0.2000),
                (0.20005, 0.20005),
                (0.2000, 0.20005),
                (0.2000, 0.2000),
            ),
        ),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water", "water": "lake"},
    )
    far_other = WaterPolygonFootprint(
        water_id="far-other",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (0.2000, 0.2000),
                (0.20005, 0.2000),
                (0.20005, 0.20005),
                (0.2000, 0.20005),
                (0.2000, 0.2000),
            ),
        ),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water"},
    )
    far_river = WaterPolygonFootprint(
        water_id="far-river",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (0.2000, 0.2000),
                (0.2015, 0.2000),
                (0.2015, 0.2015),
                (0.2000, 0.2015),
                (0.2000, 0.2000),
            ),
        ),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water", "water": "river"},
    )
    near_kept_lake = WaterPolygonFootprint(
        water_id="near-lake",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (0.0100, 0.0100),
                (0.0500, 0.0100),
                (0.0500, 0.0500),
                (0.0100, 0.0500),
                (0.0100, 0.0100),
            ),
        ),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water", "water": "lake"},
    )

    simplified = simplify_water_footprints_for_observer(
        (near_kept_lake, far_small_lake, far_other, far_river),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
    )

    assert [footprint.water_id for footprint in simplified] == ["near-lake", "far-river"]


def test_water_vertex_spacing_threshold_grows_with_distance() -> None:
    from zstarview import water_overlay

    near = water_overlay._water_vertex_spacing_threshold_m(  # noqa: SLF001
        1.0,
    )
    far = water_overlay._water_vertex_spacing_threshold_m(  # noqa: SLF001
        4.0,
    )

    assert near == 50.0
    assert far == 200.0
    assert far > near


def test_water_vertex_spacing_threshold_uses_nearer_pair_distance() -> None:
    from zstarview import water_overlay

    assert (
        water_overlay._water_vertex_spacing_threshold_for_pair_m(  # noqa: SLF001
            0.5,
            4.0,
        )
        == 25.0
    )


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


def test_sample_water_overlay_points_uses_ground_sampler_for_river_like_only() -> None:
    lake = WaterPolygonFootprint(
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
        tags={"natural": "water", "water": "lake"},
    )
    river = WaterPolygonFootprint(
        water_id="river",
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
        tags={"natural": "water", "water": "river"},
    )

    lake_sampler = Mock(return_value=123.0)
    river_sampler = Mock(return_value=123.0)

    sample_water_overlay_points(
        (lake,),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=100.0,
        fallback_surface_height_m=50.0,
        target_ground_elevation_m_sampler=lake_sampler,
        max_distance_km=0.2,
        sample_step_m=100.0,
        azimuth_step_deg=90.0,
    )
    sample_water_overlay_points(
        (river,),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=100.0,
        fallback_surface_height_m=50.0,
        target_ground_elevation_m_sampler=river_sampler,
        max_distance_km=0.2,
        sample_step_m=100.0,
        azimuth_step_deg=90.0,
    )

    assert lake_sampler.call_count == 0
    assert river_sampler.call_count > 0


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


def test_resolve_water_scan_radius_scales_with_height() -> None:
    low = resolve_water_scan_radius_km(0.0)
    high = resolve_water_scan_radius_km(500.0)

    assert low == 2.0
    assert high > low
    assert high == horizon_distance_km_from_height(500.0) + 1.0


def test_expanded_query_bbox_from_point_scales_by_20_percent() -> None:
    view_bbox = expanded_query_bbox_from_point(35.0, 139.0, 5.0, scale=1.0)
    query_bbox = expanded_query_bbox_from_point(35.0, 139.0, 5.0)

    view_width = view_bbox[2] - view_bbox[0]
    view_height = view_bbox[3] - view_bbox[1]
    query_width = query_bbox[2] - query_bbox[0]
    query_height = query_bbox[3] - query_bbox[1]

    assert math.isclose(query_width, view_width * 1.2, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(query_height, view_height * 1.2, rel_tol=0.0, abs_tol=1e-12)


def test_fetch_overpass_json_reports_compact_http_error(monkeypatch) -> None:
    from zstarview import water_overlay

    def fake_urlopen(*_args, **_kwargs):
        raise urllib.error.HTTPError(
            url="https://overpass-api.de/api/interpreter",
            code=504,
            msg="Gateway Timeout",
            hdrs=None,
            fp=Mock(),
        )

    monkeypatch.setattr(water_overlay.urllib.request, "urlopen", fake_urlopen)

    with pytest.raises(RuntimeError, match=r"^HTTP 504$"):
        water_overlay.fetch_overpass_json(bbox=(0.0, 0.0, 1.0, 1.0))


def test_fetch_overpass_json_reports_timeout_as_compact_error(monkeypatch) -> None:
    from zstarview import water_overlay

    def fake_urlopen(*_args, **_kwargs):
        raise urllib.error.URLError(TimeoutError("timed out"))

    monkeypatch.setattr(water_overlay.urllib.request, "urlopen", fake_urlopen)

    with pytest.raises(RuntimeError, match=r"^timeout$"):
        water_overlay.fetch_overpass_json(bbox=(0.0, 0.0, 1.0, 1.0))


def test_fetch_overpass_json_can_be_cancelled(monkeypatch) -> None:
    from zstarview import water_overlay

    abort_event = threading.Event()
    abort_event.set()

    with pytest.raises(DownloadCancelledError):
        water_overlay.fetch_overpass_json(
            bbox=(0.0, 0.0, 1.0, 1.0),
            abort_event=abort_event,
        )


def test_sample_water_overlay_points_can_be_cancelled() -> None:
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
    abort_event = threading.Event()
    abort_event.set()

    with pytest.raises(DownloadCancelledError):
        sample_water_overlay_points(
            (footprint,),
            observer_lat_deg=0.0,
            observer_lon_deg=0.0,
            observer_height_m=100.0,
            fallback_surface_height_m=0.0,
            max_distance_km=0.2,
            sample_step_m=100.0,
            azimuth_step_deg=90.0,
            abort_event=abort_event,
        )


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


def test_water_overlay_point_color_rgb_is_unified_for_sea_bands() -> None:
    assert _water_overlay_point_color_rgb(  # noqa: SLF001
        WaterOverlayPoint("a", 0.0, 0.0, 0.0, water_category="sea-125")
    ) == _water_overlay_point_color_rgb(  # noqa: SLF001
        WaterOverlayPoint("b", 0.0, 0.0, 0.0, water_category="sea-250")
    )
    assert _water_overlay_point_color_rgb(  # noqa: SLF001
        WaterOverlayPoint("c", 0.0, 0.0, 0.0, water_category="sea-500")
    ) == _water_overlay_point_color_rgb(  # noqa: SLF001
        WaterOverlayPoint("d", 0.0, 0.0, 0.0, water_category="sea")
    )


def test_water_overlay_point_color_rgb_distinguishes_sea_and_inland_water() -> None:
    sea_color = _water_overlay_point_color_rgb(  # noqa: SLF001
        WaterOverlayPoint("sea", 0.0, 0.0, 0.0, water_category="sea-125")
    )
    river_color = _water_overlay_point_color_rgb(  # noqa: SLF001
        WaterOverlayPoint("river", 0.0, 0.0, 0.0, water_category="river")
    )
    lake_color = _water_overlay_point_color_rgb(  # noqa: SLF001
        WaterOverlayPoint("lake", 0.0, 0.0, 0.0, water_category="lake")
    )

    assert river_color != sea_color
    assert lake_color != sea_color
    assert river_color != lake_color


def test_build_overpass_query_excludes_coastline() -> None:
    query = build_overpass_query((0.0, 1.0, 2.0, 3.0))

    assert 'natural"="coastline' not in query
    assert 'natural"="water' in query
    assert 'waterway"="riverbank' in query
