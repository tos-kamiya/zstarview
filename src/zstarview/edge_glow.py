from __future__ import annotations

from typing import Sequence

import numpy as np

def _terrain_sample_edge_strength_rows(
    *,
    terrain_sample_distances_m: Sequence[float] | Sequence[Sequence[float]] | np.ndarray | None,
    terrain_sample_terrain_elevation_m: Sequence[Sequence[float]] | np.ndarray | None,
    source_distances_m: np.ndarray,
) -> np.ndarray | None:
    if terrain_sample_distances_m is None or terrain_sample_terrain_elevation_m is None:
        return None
    terrain_distances_arr = np.asarray(terrain_sample_distances_m, dtype=np.float64)
    terrain_elevation = np.asarray(terrain_sample_terrain_elevation_m, dtype=np.float64)
    if terrain_distances_arr.size == 0 or terrain_elevation.size == 0:
        return None
    if terrain_elevation.ndim != 2:
        raise ValueError("terrain_sample_terrain_elevation_m must be a 2D array")
    if terrain_distances_arr.ndim == 1:
        terrain_distances = terrain_distances_arr.reshape(-1)
        if terrain_elevation.shape[1] != terrain_distances.size:
            raise ValueError(
                "terrain_sample_terrain_elevation_m must match terrain_sample_distances_m"
            )
    elif terrain_distances_arr.ndim == 2:
        if terrain_distances_arr.shape != terrain_elevation.shape:
            raise ValueError(
                "terrain_sample_distances_m must match terrain_sample_terrain_elevation_m"
            )
        terrain_distances = np.asarray(terrain_distances_arr[0], dtype=np.float64).reshape(-1)
        if not np.allclose(terrain_distances_arr, terrain_distances[np.newaxis, :], equal_nan=True):
            terrain_distances = None
    else:
        raise ValueError("terrain_sample_distances_m must be 1D or 2D")
    source_distances = np.asarray(source_distances_m, dtype=np.float64).reshape(-1)
    if source_distances.size == 0:
        return np.zeros((terrain_elevation.shape[0], 0), dtype=np.float64)
    terrain_rel_elevation = np.ones_like(np.asarray(terrain_elevation, dtype=np.float64), dtype=np.float64)

    if terrain_distances is not None and np.array_equal(terrain_distances, source_distances):
        return terrain_rel_elevation
    rows = [
        np.interp(
            source_distances,
            terrain_distances if terrain_distances is not None else np.asarray(terrain_distances_arr[row_index], dtype=np.float64).reshape(-1),
            np.asarray(row, dtype=np.float64),
            left=float(row[0]),
            right=float(row[-1]),
        )
        for row_index, row in enumerate(np.asarray(terrain_rel_elevation, dtype=np.float64))
    ]
    return np.asarray(rows, dtype=np.float64)



