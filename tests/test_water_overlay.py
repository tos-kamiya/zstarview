from __future__ import annotations

import math
import threading
import urllib.error
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from PySide6.QtCore import Qt

from zstarview.clouddisc.types import DownloadCancelledError
from zstarview.render.geometry import ScreenGeometry
from zstarview.render.terrain import (
    _terrain_occlusion_alpha_scale,
    _thin_water_overlay_dots_pairwise,
    _water_overlay_distance_alpha_scale,
    _water_overlay_marker_geometry,
    _water_overlay_point_color_rgb,
    apply_terrain_occlusion_to_water_points,
    draw_water_overlay_dots,
    draw_water_overlay_polylines,
)
from zstarview.types import ViewerData
from zstarview.water_overlay import (
    WaterOverlayPoint,
    WaterOverlayPolyline,
    WaterPolygonFootprint,
    WaterSurfacePatch,
    assemble_rings_from_segments,
    build_geometric_distance_samples,
    build_overpass_query,
    build_water_overlay_polylines,
    classify_water_surface_category,
    classify_water_surface_mode,
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


def test_water_overlay_point_color_rgb_distinguishes_sea_125_and_inland_water() -> None:
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
    assert river_color == lake_color


def test_water_overlay_marker_geometry_flattens_the_marker() -> None:
    major_radius, minor_radius, pen_width = _water_overlay_marker_geometry(1.0)

    assert major_radius > minor_radius
    assert pen_width >= 1.0


def test_water_overlay_marker_geometry_shrinks_with_distance() -> None:
    near_major, near_minor, _ = _water_overlay_marker_geometry(1.0, distance_m=0.0)
    far_major, far_minor, _ = _water_overlay_marker_geometry(1.0, distance_m=256_000.0)

    assert near_major > far_major
    assert near_minor > far_minor
    assert (near_minor / near_major) < 0.5
    assert far_minor == pytest.approx(0.6)


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


def test_draw_water_overlay_dots_uses_unfilled_ellipse_marker() -> None:
    class PainterStub:
        def __init__(self) -> None:
            self.calls: list[tuple[object, ...]] = []

        def save(self) -> None:
            self.calls.append(("save",))

        def restore(self) -> None:
            self.calls.append(("restore",))

        def setPen(self, pen) -> None:
            self.calls.append(("setPen", pen))

        def setBrush(self, brush) -> None:
            self.calls.append(("setBrush", brush))

        def translate(self, x, y) -> None:
            self.calls.append(("translate", x, y))

        def drawEllipse(self, center, rx, ry) -> None:
            self.calls.append(("drawEllipse", center, rx, ry))

        def drawLine(self, start, end) -> None:
            self.calls.append(("drawLine", start, end))

    painter = PainterStub()
    geometry = ScreenGeometry(center=(100, 80), radius=60)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Test",
        view_center=(45.0, 90.0),
        edge_fov_deg=95.0,
        content_fov_deg=110.0,
    )
    water_dots = [
        WaterOverlayPoint("near", 10.0, 20.0, 0.5, scan_distance_m=500.0, water_category="lake"),
        WaterOverlayPoint("far", 12.0, 40.0, 256.0, scan_distance_m=256_000.0, water_category="lake"),
    ]

    draw_water_overlay_dots(
        painter,
        geometry,
        viewer,
        water_dots,
        opacity=0.5,
        line_width_scale=1.0,
        pairwise_thinning=False,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, *_args, **_kwargs: (
            (0.3, -0.1) if float(az) < 30.0 else (-0.4, 0.2)
        ),
        normalized_to_screen_xy_func=lambda nx, ny, geometry: (
            geometry.center[0] + nx * geometry.radius,
            geometry.center[1] + ny * geometry.radius,
        ),
    )

    draw_ellipse_calls = [call for call in painter.calls if call[0] == "drawEllipse"]
    draw_line_calls = [call for call in painter.calls if call[0] == "drawLine"]
    translate_calls = [call for call in painter.calls if call[0] == "translate"]
    set_brush_calls = [call for call in painter.calls if call[0] == "setBrush"]
    set_pen_calls = [call for call in painter.calls if call[0] == "setPen"]

    assert draw_ellipse_calls
    assert not draw_line_calls
    assert len(draw_ellipse_calls) == 2
    assert len(translate_calls) == 2
    assert translate_calls[0][1:] == (118.0, 74.0)
    assert translate_calls[1][1:] == (76.0, 92.0)
    for _, center, rx, ry in draw_ellipse_calls:
        assert center.x() == pytest.approx(0.0)
        assert center.y() == pytest.approx(0.0)
        assert rx > ry
    assert draw_ellipse_calls[0][2] > draw_ellipse_calls[1][2]
    assert draw_ellipse_calls[0][3] > draw_ellipse_calls[1][3]
    assert set_brush_calls[0][1] == Qt.BrushStyle.NoBrush
    assert set_pen_calls[0][1].color().alpha() > 0


def test_draw_water_overlay_dots_uses_fast_mode_unfilled_ellipse() -> None:
    class PainterStub:
        def __init__(self) -> None:
            self.calls: list[tuple[object, ...]] = []

        def save(self) -> None:
            self.calls.append(("save",))

        def restore(self) -> None:
            self.calls.append(("restore",))

        def setPen(self, pen) -> None:
            self.calls.append(("setPen", pen))

        def setBrush(self, brush) -> None:
            self.calls.append(("setBrush", brush))

        def translate(self, x, y) -> None:
            self.calls.append(("translate", x, y))

        def drawEllipse(self, center, rx, ry) -> None:
            self.calls.append(("drawEllipse", center, rx, ry))

        def drawLine(self, start, end) -> None:
            self.calls.append(("drawLine", start, end))

    painter = PainterStub()
    geometry = ScreenGeometry(center=(100, 80), radius=60)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Test",
        view_center=(45.0, 90.0),
        edge_fov_deg=95.0,
        content_fov_deg=110.0,
    )
    water_dots = [
        WaterOverlayPoint("near", 10.0, 20.0, 0.5, scan_distance_m=500.0, water_category="lake"),
    ]

    draw_water_overlay_dots(
        painter,
        geometry,
        viewer,
        water_dots,
        opacity=0.5,
        line_width_scale=1.0,
        pairwise_thinning=False,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, *_args, **_kwargs: (0.3, -0.1),
        normalized_to_screen_xy_func=lambda nx, ny, geometry: (
            geometry.center[0] + nx * geometry.radius,
            geometry.center[1] + ny * geometry.radius,
        ),
    )

    draw_ellipse_calls = [call for call in painter.calls if call[0] == "drawEllipse"]
    draw_line_calls = [call for call in painter.calls if call[0] == "drawLine"]
    set_brush_calls = [call for call in painter.calls if call[0] == "setBrush"]
    set_pen_calls = [call for call in painter.calls if call[0] == "setPen"]

    assert draw_ellipse_calls
    assert not draw_line_calls
    assert len(draw_ellipse_calls) == 1
    _, center, rx, ry = draw_ellipse_calls[0]
    assert center.x() == pytest.approx(0.0)
    assert center.y() == pytest.approx(0.0)
    assert rx > ry
    assert set_brush_calls[-1][1] == Qt.BrushStyle.NoBrush
    assert set_pen_calls[-1][1].color().alpha() > 0


def test_water_overlay_distance_alpha_scale_decays_with_distance() -> None:
    assert _water_overlay_distance_alpha_scale(0.0) == 1.0
    assert _water_overlay_distance_alpha_scale(128_000.0) == pytest.approx(0.0625, rel=1e-6)
    assert _water_overlay_distance_alpha_scale(256_000.0) == pytest.approx(0.00390625, rel=1e-6)


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

    got = _thin_water_overlay_dots_pairwise(points)

    assert [point.scan_distance_index for point in got if point.scan_azimuth_index == 0] == [0, 2]
    assert [point.scan_distance_index for point in got if point.scan_azimuth_index == 1] == [1, 3]
    assert any(point.scan_distance_index is None for point in got)


def test_water_overlay_point_color_rgb_is_unified_for_sea_bands() -> None:
    assert _water_overlay_point_color_rgb(
        WaterOverlayPoint("a", 0.0, 0.0, 0.0, water_category="sea-125")
    ) == _water_overlay_point_color_rgb(
        WaterOverlayPoint("b", 0.0, 0.0, 0.0, water_category="sea-250")
    )
    assert _water_overlay_point_color_rgb(
        WaterOverlayPoint("c", 0.0, 0.0, 0.0, water_category="sea-500")
    ) == _water_overlay_point_color_rgb(
        WaterOverlayPoint("d", 0.0, 0.0, 0.0, water_category="sea")
    )


def test_water_overlay_point_color_rgb_distinguishes_sea_and_inland_water() -> None:
    sea_color = _water_overlay_point_color_rgb(
        WaterOverlayPoint("sea", 0.0, 0.0, 0.0, water_category="sea-125")
    )
    river_color = _water_overlay_point_color_rgb(
        WaterOverlayPoint("river", 0.0, 0.0, 0.0, water_category="river")
    )
    lake_color = _water_overlay_point_color_rgb(
        WaterOverlayPoint("lake", 0.0, 0.0, 0.0, water_category="lake")
    )

    assert river_color != sea_color
    assert lake_color != sea_color
    assert river_color == lake_color


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
