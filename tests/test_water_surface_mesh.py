from __future__ import annotations

from zstarview.water_overlay import WaterPolygonFootprint, WaterSurfacePatch
from zstarview.water_surface_mesh import build_water_surface_mesh, make_local_transformer, project_ring_xy


def _ring_area(points: list[tuple[float, float]]) -> float:
    if not points:
        return 0.0
    ring = list(points)
    if ring[0] != ring[-1]:
        ring.append(ring[0])
    area = 0.0
    for index in range(len(ring) - 1):
        x0, y0 = ring[index]
        x1, y1 = ring[index + 1]
        area += x0 * y1 - x1 * y0
    return abs(area) * 0.5


def _triangle_area(triangle: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]) -> float:
    (x0, y0), (x1, y1), (x2, y2) = triangle
    return abs((x0 * (y1 - y2)) + (x1 * (y2 - y0)) + (x2 * (y0 - y1))) * 0.5


def test_build_water_surface_mesh_creates_triangles_for_simple_polygon() -> None:
    footprint = WaterPolygonFootprint(
        water_id="way/1",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (133.0, 35.0),
                (133.01, 35.0),
                (133.01, 35.01),
                (133.0, 35.01),
                (133.0, 35.0),
            ),
        ),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water"},
    )
    patch = WaterSurfacePatch(
        patch_id="patch/1",
        polygon_id="way/1",
        anchor_points_lonlat=((133.0, 35.0), (133.01, 35.0), (133.0, 35.01)),
        anchor_elevations_m=(10.0, 10.2, 10.1),
    )

    mesh = build_water_surface_mesh(
        footprint,
        center_lat_deg=35.005,
        center_lon_deg=133.005,
        patch=patch,
        grid_m=1.0,
        simplify_tolerance_m=1.0,
    )

    assert mesh is not None
    assert mesh.water_id == "way/1"
    assert mesh.surface_mode == "flat"
    assert mesh.split_cell_m == 300.0
    assert len(mesh.triangles_xy_m) > 2

    transformer = make_local_transformer(35.005, 133.005)
    outer_xy = project_ring_xy(footprint.outer_rings_lonlat[0], transformer)
    expected_area = _ring_area(outer_xy)
    got_area = sum(_triangle_area(triangle) for triangle in mesh.triangles_xy_m)

    assert abs(got_area - expected_area) / expected_area < 0.001


def test_build_water_surface_mesh_preserves_hole_area() -> None:
    footprint = WaterPolygonFootprint(
        water_id="relation/2",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (133.0, 35.0),
                (133.02, 35.0),
                (133.02, 35.02),
                (133.0, 35.02),
                (133.0, 35.0),
            ),
        ),
        inner_rings_lonlat=(
            (
                (133.005, 35.005),
                (133.01, 35.005),
                (133.01, 35.01),
                (133.005, 35.01),
                (133.005, 35.005),
            ),
        ),
        source="relation",
        tags={"natural": "water"},
    )

    mesh = build_water_surface_mesh(
        footprint,
        center_lat_deg=35.01,
        center_lon_deg=133.01,
        grid_m=1.0,
        simplify_tolerance_m=1.0,
    )

    assert mesh is not None
    assert len(mesh.triangles_xy_m) >= 2

    transformer = make_local_transformer(35.01, 133.01)
    outer_xy = project_ring_xy(footprint.outer_rings_lonlat[0], transformer)
    hole_xy = project_ring_xy(footprint.inner_rings_lonlat[0], transformer)
    expected_area = _ring_area(outer_xy) - _ring_area(hole_xy)
    got_area = sum(_triangle_area(triangle) for triangle in mesh.triangles_xy_m)

    assert abs(got_area - expected_area) / expected_area < 0.01


def test_build_water_surface_mesh_can_split_triangles_by_grid() -> None:
    footprint = WaterPolygonFootprint(
        water_id="way/3",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (133.0, 35.0),
                (133.02, 35.0),
                (133.02, 35.02),
                (133.0, 35.02),
                (133.0, 35.0),
            ),
        ),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water"},
    )

    mesh = build_water_surface_mesh(
        footprint,
        center_lat_deg=35.01,
        center_lon_deg=133.01,
        grid_m=1.0,
        simplify_tolerance_m=1.0,
    )

    assert mesh is not None
    assert mesh.split_cell_m == 300.0
    assert len(mesh.triangles_xy_m) > 2

    transformer = make_local_transformer(35.01, 133.01)
    outer_xy = project_ring_xy(footprint.outer_rings_lonlat[0], transformer)
    expected_area = _ring_area(outer_xy)
    got_area = sum(_triangle_area(triangle) for triangle in mesh.triangles_xy_m)

    assert abs(got_area - expected_area) / expected_area < 0.01
