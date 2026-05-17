from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


COARSE_GRID_ROWS = 4
COARSE_GRID_COLS = 8

DEFAULT_WATER_TILES_ROOT_125M = Path(__file__).resolve().parent / "data" / "water_tiles_125m"
DEFAULT_WATER_TILES_ROOT_250M = Path(__file__).resolve().parent / "data" / "water_tiles_250m"
DEFAULT_WATER_TILES_ROOT_500M = Path(__file__).resolve().parent / "data" / "water_tiles_500m"


@dataclass(frozen=True, slots=True)
class TileGridSpec:
    label: str
    root: Path
    factor: int
    rows: int
    cols: int
    tiles: dict[tuple[int, int], Path]


@dataclass(frozen=True, slots=True)
class BoundaryIntersection:
    coarse_row_index: int
    coarse_col_index: int
    latitude_deg: float
    longitude_deg: float
    tile_paths_by_root: dict[str, tuple[Path, Path, Path, Path]]


def _parse_tile_key(tile_path: Path) -> tuple[int, int] | None:
    parts = tile_path.stem.split("_")
    if len(parts) != 3:
        return None
    row_part = parts[1]
    col_part = parts[2]
    if not row_part.startswith("y") or not col_part.startswith("x"):
        return None
    try:
        row = int(row_part[1:])
        col = int(col_part[1:])
    except ValueError:
        return None
    return row, col


def _load_tile_map(tile_root: Path) -> tuple[dict[tuple[int, int], Path], int, int]:
    if not tile_root.exists():
        return {}, 0, 0

    tiles: dict[tuple[int, int], Path] = {}
    for tile_path in sorted(tile_root.glob("tile_y*_x*")):
        if tile_path.suffix not in {".tif", ".0", ".1"}:
            continue
        tile_key = _parse_tile_key(tile_path)
        if tile_key is None:
            continue
        tiles[tile_key] = tile_path

    if not tiles:
        return {}, 0, 0

    max_row = max(row for row, _col in tiles)
    max_col = max(col for _row, col in tiles)
    return tiles, max_row + 1, max_col + 1


def load_tile_grid_spec(label: str, tile_root: Path) -> TileGridSpec:
    tiles, rows, cols = _load_tile_map(tile_root)
    if rows == 0 or cols == 0:
        raise ValueError(f"No water tiles found in {tile_root}")
    if rows % COARSE_GRID_ROWS != 0 or cols % COARSE_GRID_COLS != 0:
        raise ValueError(
            f"Tile grid in {tile_root} does not align with the coarse 4x8 boundary grid"
        )
    row_factor = rows // COARSE_GRID_ROWS
    col_factor = cols // COARSE_GRID_COLS
    if row_factor != col_factor:
        raise ValueError(
            f"Tile grid in {tile_root} is not square enough for the shared boundary check"
        )
    return TileGridSpec(
        label=label,
        root=tile_root,
        factor=row_factor,
        rows=rows,
        cols=cols,
        tiles=tiles,
    )


def _coarse_boundary_coordinate(coarse_row_index: int, coarse_col_index: int) -> tuple[float, float]:
    latitude_deg = 90.0 - (float(coarse_row_index) * 45.0)
    longitude_deg = -180.0 + (float(coarse_col_index) * 45.0)
    return latitude_deg, longitude_deg


def _boundary_quad_paths(
    spec: TileGridSpec,
    *,
    coarse_row_index: int,
    coarse_col_index: int,
) -> tuple[Path, Path, Path, Path] | None:
    row = coarse_row_index * spec.factor
    col = coarse_col_index * spec.factor
    if row <= 0 or row >= spec.rows:
        return None
    if col <= 0 or col >= spec.cols:
        return None

    tiles = spec.tiles
    quad_keys = ((row - 1, col - 1), (row - 1, col), (row, col - 1), (row, col))
    quad_paths: list[Path] = []
    for key in quad_keys:
        tile_path = tiles.get(key)
        if tile_path is None or tile_path.suffix != ".tif":
            return None
        quad_paths.append(tile_path)
    return tuple(quad_paths)


def find_common_boundary_intersections(
    specs: Sequence[TileGridSpec] | None = None,
) -> list[BoundaryIntersection]:
    grid_specs = list(specs) if specs is not None else [
        load_tile_grid_spec("125m", DEFAULT_WATER_TILES_ROOT_125M),
        load_tile_grid_spec("250m", DEFAULT_WATER_TILES_ROOT_250M),
        load_tile_grid_spec("500m", DEFAULT_WATER_TILES_ROOT_500M),
    ]

    matches: list[BoundaryIntersection] = []
    for coarse_row_index in range(1, COARSE_GRID_ROWS):
        for coarse_col_index in range(1, COARSE_GRID_COLS):
            tile_paths_by_root: dict[str, tuple[Path, Path, Path, Path]] = {}
            for spec in grid_specs:
                quad_paths = _boundary_quad_paths(
                    spec,
                    coarse_row_index=coarse_row_index,
                    coarse_col_index=coarse_col_index,
                )
                if quad_paths is None:
                    break
                tile_paths_by_root[spec.label] = quad_paths
            else:
                latitude_deg, longitude_deg = _coarse_boundary_coordinate(
                    coarse_row_index,
                    coarse_col_index,
                )
                matches.append(
                    BoundaryIntersection(
                        coarse_row_index=coarse_row_index,
                        coarse_col_index=coarse_col_index,
                        latitude_deg=latitude_deg,
                        longitude_deg=longitude_deg,
                        tile_paths_by_root=tile_paths_by_root,
                    )
                )
    return matches

