from __future__ import annotations

from zstarview.water_overlay import WaterPolygonFootprint, WaterSurfacePatch
from zstarview.water_surface_mesh import build_water_surface_mesh


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
    assert len(mesh.triangles_xy_m) >= 1


def test_build_water_surface_mesh_preserves_hole_triangulation() -> None:
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
    assert len(mesh.triangles_xy_m) >= 1
