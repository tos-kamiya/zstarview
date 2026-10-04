from __future__ import annotations

import math
import threading
import urllib.error
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from zstarview.clouddisc.types import DownloadCancelledError
from zstarview.render.geometry import ScreenGeometry
from zstarview.render.terrain import (
    WATER_SURFACE_WAVE_FAR_DECAY_DISTANCE_KM,
    WATER_SURFACE_WAVE_MAX_AMPLITUDE_PX,
    WATER_SURFACE_WAVE_NEAREST_ENVELOPE_DISTANCE_KM,
    _terrain_occlusion_alpha_scale,
    _water_surface_wave_offset_px,
    _water_surface_wave_phases,
    apply_terrain_occlusion_to_water_points,
    draw_water_overlay_polylines,
)
from zstarview.types import ViewerData
from zstarview.water_overlay import (
    WaterOverlayPoint,
    WaterOverlayPolyline,
    WaterPolygonFootprint,
    assemble_rings_from_segments,
    build_geometric_distance_samples,
    build_overpass_query,
    build_water_overlay_polylines,
    classify_water_surface_category,
    expanded_query_bbox_from_point,
    extract_water_polygons,
    resolve_water_scan_radius_km,
    resolve_water_surface_azimuth_step_deg,
    sample_water_overlay_points,
    simplify_water_footprints_for_observer,
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


def test_classify_water_surface_category_uses_tags_and_kind() -> None:
    assert classify_water_surface_category({"natural": "coastline"}) == "sea"
    assert classify_water_surface_category({"water": "lake"}) == "lake"
    assert classify_water_surface_category({"water": "river"}) == "river"
    assert classify_water_surface_category({"waterway": "riverbank"}) == "river"


def test_terrain_occlusion_fades_water_behind_nearer_higher_terrain() -> None:
    assert _terrain_occlusion_alpha_scale(
        1.0,
        90.0,
        20_000.0,
        [(2.0, 90.0)],
        [10_000.0],
    ) == pytest.approx(0.48)


def test_terrain_occlusion_keeps_water_in_front_or_above_terrain_visible() -> None:
    profile = [(2.0, 90.0)]
    distances = [10_000.0]
    assert _terrain_occlusion_alpha_scale(1.0, 90.0, 5_000.0, profile, distances) == 1.0
    assert _terrain_occlusion_alpha_scale(3.0, 90.0, 20_000.0, profile, distances) == 1.0
    assert _terrain_occlusion_alpha_scale(1.0, 92.5, 20_000.0, profile, distances) == 1.0


def test_runtime_water_points_store_terrain_occlusion_alpha() -> None:
    points = apply_terrain_occlusion_to_water_points(
        (WaterOverlayPoint("water", 1.0, 90.0, 20.0, scan_distance_m=20_000.0),),
        [(2.0, 90.0)],
        [10_000.0],
    )

    assert points[0].terrain_occlusion_alpha_scale == pytest.approx(0.48)


def test_water_polyline_keeps_segments_connected_across_terrain_alpha_changes() -> None:
    class PainterStub:
        def __init__(self) -> None:
            self.polylines: list[list[object]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, _pen) -> None:
            pass

        def setBrush(self, _brush) -> None:
            pass

        def drawPolyline(self, polyline) -> None:
            self.polylines.append(list(polyline))

    painter = PainterStub()
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Test",
        view_center=(0.0, 0.0),
        edge_fov_deg=180.0,
        content_fov_deg=180.0,
    )
    points = (
        WaterOverlayPoint("water", 0.0, 0.0, 1.0, terrain_occlusion_alpha_scale=1.0),
        WaterOverlayPoint("water", 0.0, 1.0, 1.0, terrain_occlusion_alpha_scale=0.48),
        WaterOverlayPoint("water", 0.0, 2.0, 1.0, terrain_occlusion_alpha_scale=1.0),
    )

    draw_water_overlay_polylines(
        painter,
        ScreenGeometry(center=(100, 100), radius=100),
        viewer,
        [WaterOverlayPolyline("water", "lake", points)],
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda _alt, az, *_args, **_kwargs: (az, 0.0),
        normalized_to_screen_xy_func=lambda x, y, _geometry: (x, y),
    )

    assert [len(polyline) for polyline in painter.polylines] == [2, 2]


def test_water_line_width_and_alpha_follow_distance_curve() -> None:
    class PainterStub:
        def __init__(self) -> None:
            self.pen = None
            self.strokes: list[tuple[float, int]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            self.pen = pen

        def setBrush(self, _brush) -> None:
            pass

        def drawPolyline(self, _polyline) -> None:
            self.strokes.append((self.pen.widthF(), self.pen.color().alpha()))

    painter = PainterStub()
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Test",
        view_center=(0.0, 0.0),
        edge_fov_deg=180.0,
        content_fov_deg=180.0,
    )
    distances_m = (500.0, 2_500.0, 4_500.0, 40_500.0)
    rings = [
        WaterOverlayPolyline(
            f"ring-{index}",
            "lake",
            (
                WaterOverlayPoint(
                    f"ring-{index}",
                    0.0,
                    0.0,
                    distance_m / 1000.0,
                    scan_distance_m=distance_m,
                ),
                WaterOverlayPoint(
                    f"ring-{index}",
                    0.0,
                    1.0,
                    distance_m / 1000.0,
                    scan_distance_m=distance_m,
                ),
            ),
        )
        for index, distance_m in enumerate(distances_m)
    ]

    draw_water_overlay_polylines(
        painter,
        ScreenGeometry(center=(100, 100), radius=100),
        viewer,
        rings,
        opacity=0.8,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda _alt, az, *_args, **_kwargs: (az, 0.0),
        normalized_to_screen_xy_func=lambda x, y, _geometry: (x, y),
    )

    assert [alpha for _width, alpha in painter.strokes] == [41, 122, 163, 204]
    expected_widths = [
        1.35 * 0.8 * (500.0 / distance_m) ** 0.25
        for distance_m in distances_m
    ]
    assert [width for width, _alpha in painter.strokes] == pytest.approx(
        expected_widths
    )


def test_water_surface_wave_is_bounded_periodic_and_fades_with_distance() -> None:
    phases = _water_surface_wave_phases("sea/1", "sea", 500.0)
    near_offset = _water_surface_wave_offset_px(
        6.0, 0.5, phases, nearest_distance_km=0.5
    )

    assert abs(near_offset) <= WATER_SURFACE_WAVE_MAX_AMPLITUDE_PX * math.exp(
        -0.5 / WATER_SURFACE_WAVE_NEAREST_ENVELOPE_DISTANCE_KM
    )
    assert _water_surface_wave_offset_px(
        366.0, 0.5, phases, nearest_distance_km=0.5
    ) == pytest.approx(near_offset)
    assert _water_surface_wave_offset_px(
        6.0, 1.5, phases, nearest_distance_km=0.5
    ) == pytest.approx(
        near_offset * math.exp(-1.0 / WATER_SURFACE_WAVE_FAR_DECAY_DISTANCE_KM)
    )
    one_km_max_offset = _water_surface_wave_offset_px(
        0.0,
        1.0,
        (math.pi / 2.0, math.pi / 2.0),
        nearest_distance_km=1.0,
    )
    assert one_km_max_offset == pytest.approx(
        WATER_SURFACE_WAVE_MAX_AMPLITUDE_PX
        * math.exp(-1.0 / WATER_SURFACE_WAVE_NEAREST_ENVELOPE_DISTANCE_KM)
    )
    assert _water_surface_wave_phases("sea/1", "sea", 500.0) == phases
    assert _water_surface_wave_phases("sea/1", "sea", 2_500.0) != phases


def test_water_surface_wave_offsets_rings_only_outside_fast_mode() -> None:
    class PainterStub:
        def __init__(self) -> None:
            self.lines: list[list[object]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, _pen) -> None:
            pass

        def setBrush(self, _brush) -> None:
            pass

        def drawPolyline(self, line) -> None:
            self.lines.append(list(line))

    painter = PainterStub()
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Test",
        view_center=(0.0, 0.0),
        edge_fov_deg=180.0,
        content_fov_deg=180.0,
    )
    azimuths = (0.0, 2.0, 4.0, 6.0, 8.0)
    ring_points = tuple(
        WaterOverlayPoint(
            "sea/1",
            0.0,
            azimuth,
            0.5,
            scan_distance_m=500.0,
        )
        for azimuth in azimuths
    )
    outline_points = tuple(
        WaterOverlayPoint("outline/1", 0.0, azimuth, 0.5)
        for azimuth in azimuths
    )

    draw_water_overlay_polylines(
        painter,
        ScreenGeometry(center=(100, 100), radius=100),
        viewer,
        [
            WaterOverlayPolyline("sea/1", "sea", ring_points),
            WaterOverlayPolyline("outline/1", "lake", outline_points),
        ],
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda _alt, az, *_args, **_kwargs: (az, 0.0),
        normalized_to_screen_xy_func=lambda x, y, _geometry: (x, 100.0 + y),
    )

    ring_y = [point.y() for point in painter.lines[0]]
    outline_y = [point.y() for point in painter.lines[1]]
    phases = _water_surface_wave_phases("sea/1", "sea", 500.0)
    expected_ring_y = [
        100.0
        + _water_surface_wave_offset_px(
            azimuth,
            0.5,
            phases,
            nearest_distance_km=0.5,
        )
        for azimuth in azimuths
    ]
    assert ring_y == pytest.approx(expected_ring_y)
    assert outline_y == [100.0] * len(azimuths)

    fast_painter = PainterStub()
    draw_water_overlay_polylines(
        fast_painter,
        ScreenGeometry(center=(100, 100), radius=100),
        viewer,
        [WaterOverlayPolyline("sea/1", "sea", ring_points)],
        fast_mode=True,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda _alt, az, *_args, **_kwargs: (az, 0.0),
        normalized_to_screen_xy_func=lambda x, y, _geometry: (x, 100.0 + y),
    )
    assert [point.y() for point in fast_painter.lines[0]] == [100.0] * len(
        azimuths
    )


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
        water_overlay._water_surface_height_m(
            sea,
            fallback_surface_height_m=3.0,
            latitude_deg=0.0,
            longitude_deg=0.0,
        )
        == 12.5
    )
    assert (
        water_overlay._water_surface_height_m(
            river,
            fallback_surface_height_m=3.0,
            target_ground_elevation_m_sampler=lambda *_args: 99.0,
            latitude_deg=0.0,
            longitude_deg=0.0,
        )
        == 7.25
    )


def test_build_water_overlay_polylines_projects_simplified_ring() -> None:
    footprint = WaterPolygonFootprint(
        water_id="river/1",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (139.0000, 35.0000),
                (139.0040, 35.0000),
                (139.0040, 35.0040),
                (139.0000, 35.0040),
                (139.0000, 35.0000),
            ),
        ),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water", "water": "river"},
    )

    polylines = build_water_overlay_polylines(
        (footprint,),
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        observer_height_m=100.0,
        max_distance_km=2.0,
    )

    assert len(polylines) == 1
    assert polylines[0].water_id == "river/1/ring-0-0"
    assert polylines[0].water_category == "river"
    assert len(polylines[0].points) == 5
    assert all(point.distance_km <= 2.0 for point in polylines[0].points)


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


def test_simplify_water_footprints_for_observer_is_direction_stable() -> None:
    footprint = WaterPolygonFootprint(
        water_id="river",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (0.0000, 0.0000),
                (0.0004, 0.0000),
                (0.0008, 0.0000),
                (0.0012, 0.0000),
                (0.0012, 0.0004),
                (0.0008, 0.0004),
                (0.0004, 0.0004),
                (0.0000, 0.0004),
                (0.0000, 0.0000),
            ),
        ),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water", "water": "river"},
    )
    simplified_forward = simplify_water_footprints_for_observer(
        (footprint,),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
    )
    reversed_footprint = WaterPolygonFootprint(
        water_id="river-reversed",
        kind="natural_water",
        outer_rings_lonlat=(tuple(reversed(footprint.outer_rings_lonlat[0])),),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water", "water": "river"},
    )
    simplified_reversed = simplify_water_footprints_for_observer(
        (reversed_footprint,),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
    )

    assert simplified_forward
    assert simplified_reversed
    assert simplified_forward[0].outer_rings_lonlat[0] == tuple(
        reversed(simplified_reversed[0].outer_rings_lonlat[0])
    )


def test_simplify_water_footprints_for_observer_drops_zero_area_ring() -> None:
    footprint = WaterPolygonFootprint(
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

    simplified = simplify_water_footprints_for_observer(
        (footprint,),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
    )

    assert simplified == ()


def test_water_simplification_grid_size_grows_in_powers_of_two() -> None:
    from zstarview import water_overlay

    assert water_overlay._grid_size_for_distance_m(0.0) == 1
    assert water_overlay._grid_size_for_distance_m(150.0) == 1
    assert water_overlay._grid_size_for_distance_m(250.0) == 2
    assert water_overlay._grid_size_for_distance_m(500.0) == 4


def test_resolve_water_surface_azimuth_step_deg_scales_with_surface_size() -> None:
    assert resolve_water_surface_azimuth_step_deg() == 2.0


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


def test_sample_water_overlay_points_uses_ground_sampler_for_inland_water() -> None:
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

    assert lake_sampler.call_count > 0
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


def test_sample_water_overlay_points_uses_fixed_inland_minimum_distance(monkeypatch) -> None:
    from zstarview import water_overlay

    captured: dict[str, np.ndarray] = {}
    ray_scan = SimpleNamespace(
        azimuths_deg=np.empty(0, dtype=np.float64),
        distance_grid_m=np.empty((0, 0), dtype=np.float64),
        ray_lon_deg=np.empty((0, 0), dtype=np.float64),
        ray_lat_deg=np.empty((0, 0), dtype=np.float64),
    )

    monkeypatch.setattr(
        water_overlay,
        "build_geometric_distance_samples",
        lambda *_args, **_kwargs: np.asarray((1.0, 4.999, 5.0, 10.0)),
    )

    def _build_ray_scan_grid(**kwargs):
        captured["distance_samples_m"] = kwargs["distance_samples_m"]
        return ray_scan

    monkeypatch.setattr(water_overlay, "build_ray_scan_grid", _build_ray_scan_grid)

    points = sample_water_overlay_points(
        (),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=0.0,
        max_distance_km=0.1,
        sample_step_m=1.0,
        azimuth_step_deg=90.0,
    )

    assert points == ()
    assert water_overlay.DEFAULT_WATER_INLAND_SAMPLE_MIN_DISTANCE_M == 5.0
    np.testing.assert_array_equal(captured["distance_samples_m"], (5.0, 10.0))


def test_resolve_water_scan_radius_scales_with_height() -> None:
    low = resolve_water_scan_radius_km(0.0)
    high = resolve_water_scan_radius_km(500.0)
    capped = resolve_water_scan_radius_km(5000.0)

    assert low == 2.0
    assert high > low
    assert high == 128.0
    assert capped == 128.0


def test_resolve_water_scan_radius_uses_power_of_two_tiers() -> None:
    assert resolve_water_scan_radius_km(0.0) == 2.0
    assert resolve_water_scan_radius_km(50.0) == 32.0


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


def test_build_overpass_query_excludes_coastline() -> None:
    query = build_overpass_query((0.0, 1.0, 2.0, 3.0))

    assert 'natural"="coastline' not in query
    assert 'natural"="water' in query
    assert 'waterway"="riverbank' in query


def test_sample_water_overlay_points_can_cull_back_half_rows(monkeypatch) -> None:
    from zstarview import water_overlay

    footprint = WaterPolygonFootprint(
        water_id="water",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (-0.0010, -0.0010),
                (0.0010, -0.0010),
                (0.0010, 0.0010),
                (-0.0010, 0.0010),
                (-0.0010, -0.0010),
            ),
        ),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water"},
    )

    ray_scan = SimpleNamespace(
        azimuths_deg=np.array([0.0, 100.0, 180.0], dtype=np.float64),
        distance_grid_m=np.array([[10.0], [10.0], [10.0]], dtype=np.float64),
        ray_lon_deg=np.array([[0.0], [0.0], [0.0]], dtype=np.float64),
        ray_lat_deg=np.array([[0.0], [0.0], [0.0]], dtype=np.float64),
    )
    project_calls: list[tuple[tuple[float, ...], tuple[float, ...]]] = []

    monkeypatch.setattr(water_overlay, "build_ray_scan_grid", lambda **_kwargs: ray_scan)
    monkeypatch.setattr(water_overlay, "_point_in_footprint", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(
        water_overlay,
        "project_place_targets_to_altaz",
        lambda **kwargs: project_calls.append(
            (
                tuple(kwargs["target_latitude_deg"]),
                tuple(kwargs["target_longitude_deg"]),
            )
        )
        or [
            SimpleNamespace(alt_deg=0.0, az_deg=0.0, distance_km=1.0)
            for _ in kwargs["target_latitude_deg"]
        ],
    )

    points = water_overlay.sample_water_overlay_points(
        (footprint,),
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=0.0,
        max_distance_km=1.0,
        sample_step_m=1.0,
        azimuth_step_deg=1.0,
        front_hemisphere_view_center=(0.0, 0.0),
        front_hemisphere_fov_deg=110.0,
    )

    assert len(points) == 2
    assert len(project_calls) == 2
