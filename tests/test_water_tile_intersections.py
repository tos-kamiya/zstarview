from __future__ import annotations

from zstarview.water_tile_intersections import (
    DEFAULT_WATER_TILES_ROOT_125M,
    DEFAULT_WATER_TILES_ROOT_250M,
    DEFAULT_WATER_TILES_ROOT_500M,
    find_common_boundary_intersections,
    load_tile_grid_spec,
)


def test_find_common_boundary_intersections() -> None:
    specs = [
        load_tile_grid_spec("125m", DEFAULT_WATER_TILES_ROOT_125M),
        load_tile_grid_spec("250m", DEFAULT_WATER_TILES_ROOT_250M),
        load_tile_grid_spec("500m", DEFAULT_WATER_TILES_ROOT_500M),
    ]
    matches = find_common_boundary_intersections(specs)

    assert [(match.latitude_deg, match.longitude_deg) for match in matches] == [
        (45.0, 0.0),
        (45.0, 45.0),
        (0.0, -90.0),
        (0.0, 45.0),
        (0.0, 135.0),
    ]
    for match in matches:
        assert set(match.tile_paths_by_root) == {"125m", "250m", "500m"}
        for quad_paths in match.tile_paths_by_root.values():
            assert len(quad_paths) == 4
            assert all(path.suffix == ".tif" for path in quad_paths)
