#!/usr/bin/env python3
"""Measure piecewise-linear star-surface mesh interpolation.

The reference transform rotates each screen sample as a 3D sky direction and
projects it back to the screen.  The candidate transform first transforms the
mesh vertices, then linearly interpolates the containing triangle.  This is a
Gate 2 experiment for the proposed mesh warp; it does not alter the renderer.

Example:
  uv run -p .venv/bin/python dev-samples/measure_star_mesh_gate.py
  uv run -p .venv/bin/python dev-samples/measure_star_mesh_gate.py \
      --columns 12 --rows 8 --fail-on-not-go
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np

from zstarview.render.star_interpolation import (
    STAR_INTERPOLATION_COVERAGE,
    _direction_to_screen,
    _rotate_about_axis,
    _screen_to_direction,
)


WIDTH_PX = 1600.0
HEIGHT_PX = 900.0
GEOMETRY_CENTER = (800.0, 450.0)
GEOMETRY_RADIUS = 450.0
EDGE_FOV_DEG = 90.0
OBSERVER_LAT_DEG = 35.0
DEFAULT_ELAPSED_SECONDS = (-30.0, 30.0)
SIDEREAL_DEG_PER_SECOND = 360.0 / 86164.0905


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--columns", type=int, default=12, help="mesh cell count across")
    parser.add_argument("--rows", type=int, default=8, help="mesh cell count down")
    parser.add_argument(
        "--max-error-px",
        type=float,
        default=1.0,
        help="maximum error allowed for GO (default: 1.0)",
    )
    parser.add_argument(
        "--fail-on-not-go",
        action="store_true",
        help="return status 1 when the measured result is NOT GO",
    )
    return parser.parse_args(argv)


def _screen_grid(columns: int, rows: int) -> np.ndarray:
    x, y = np.meshgrid(
        np.linspace(0.04 * WIDTH_PX, 0.96 * WIDTH_PX, 25),
        np.linspace(0.04 * HEIGHT_PX, 0.96 * HEIGHT_PX, 15),
    )
    points = np.column_stack((x.ravel(), y.ravel()))
    cell_width = WIDTH_PX / columns
    cell_height = HEIGHT_PX / rows
    cell_x = np.minimum(columns - 1, (points[:, 0] / cell_width).astype(int))
    cell_y = np.minimum(rows - 1, (points[:, 1] / cell_height).astype(int))
    corner_x = np.column_stack(
        (
            cell_x * cell_width,
            (cell_x + 1) * cell_width,
            cell_x * cell_width,
            (cell_x + 1) * cell_width,
        )
    )
    corner_y = np.column_stack(
        (
            cell_y * cell_height,
            cell_y * cell_height,
            (cell_y + 1) * cell_height,
            (cell_y + 1) * cell_height,
        )
    )
    normalized_x = (corner_x - GEOMETRY_CENTER[0]) / GEOMETRY_RADIUS
    normalized_y = (corner_y - GEOMETRY_CENTER[1]) / GEOMETRY_RADIUS
    inside_cell = np.all(np.hypot(normalized_x, normalized_y) <= 1.0, axis=1)
    # The viewport is rectangular, while the sky-disc projection is only
    # meaningful inside the disc. Boundary clipping is measured separately.
    return points[inside_cell]


def _reference_positions(
    points: np.ndarray,
    *,
    view_center_alt_deg: float,
    view_center_az_deg: float,
    elapsed_seconds: float,
) -> np.ndarray:
    directions = _screen_to_direction(
        points[:, 0],
        points[:, 1],
        width_px=WIDTH_PX,
        height_px=HEIGHT_PX,
        geometry_center=GEOMETRY_CENTER,
        geometry_radius=GEOMETRY_RADIUS,
        view_center_altaz_deg=(view_center_alt_deg, view_center_az_deg),
        edge_fov_deg=EDGE_FOV_DEG,
    )
    latitude = math.radians(OBSERVER_LAT_DEG)
    pole_axis = np.array([math.cos(latitude), 0.0, math.sin(latitude)])
    angle = math.radians(
        SIDEREAL_DEG_PER_SECOND * elapsed_seconds * STAR_INTERPOLATION_COVERAGE
    )
    return _direction_to_screen(
        _rotate_about_axis(directions, pole_axis, angle),
        width_px=WIDTH_PX,
        height_px=HEIGHT_PX,
        geometry_center=GEOMETRY_CENTER,
        geometry_radius=GEOMETRY_RADIUS,
        view_center_altaz_deg=(view_center_alt_deg, view_center_az_deg),
        edge_fov_deg=EDGE_FOV_DEG,
    )


def _mesh_vertices(columns: int, rows: int) -> np.ndarray:
    x = np.linspace(0.0, WIDTH_PX, columns + 1)
    y = np.linspace(0.0, HEIGHT_PX, rows + 1)
    grid_x, grid_y = np.meshgrid(x, y)
    return np.column_stack((grid_x.ravel(), grid_y.ravel()))


def _mesh_interpolate(
    points: np.ndarray,
    transformed_vertices: np.ndarray,
    *,
    columns: int,
    rows: int,
) -> np.ndarray:
    cell_width = WIDTH_PX / columns
    cell_height = HEIGHT_PX / rows
    output = np.empty_like(points)
    for index, (x, y) in enumerate(points):
        cell_x = min(columns - 1, max(0, int(x / cell_width)))
        cell_y = min(rows - 1, max(0, int(y / cell_height)))
        u = (x - cell_x * cell_width) / cell_width
        v = (y - cell_y * cell_height) / cell_height
        top_left = cell_y * (columns + 1) + cell_x
        top_right = top_left + 1
        bottom_left = top_left + columns + 1
        bottom_right = bottom_left + 1
        if u + v <= 1.0:
            weights = (1.0 - u - v, u, v)
            output[index] = (
                weights[0] * transformed_vertices[top_left]
                + weights[1] * transformed_vertices[top_right]
                + weights[2] * transformed_vertices[bottom_left]
            )
        else:
            weights = (u + v - 1.0, 1.0 - u, 1.0 - v)
            output[index] = (
                weights[0] * transformed_vertices[bottom_right]
                + weights[1] * transformed_vertices[bottom_left]
                + weights[2] * transformed_vertices[top_right]
            )
    return output


def _measure_case(
    points: np.ndarray,
    *,
    view_center_alt_deg: float,
    view_center_az_deg: float,
    elapsed_seconds: float,
    columns: int,
    rows: int,
) -> tuple[float, float]:
    vertices = _mesh_vertices(columns, rows)
    transformed_vertices = _reference_positions(
        vertices,
        view_center_alt_deg=view_center_alt_deg,
        view_center_az_deg=view_center_az_deg,
        elapsed_seconds=elapsed_seconds,
    )
    mesh_positions = _mesh_interpolate(
        points, transformed_vertices, columns=columns, rows=rows
    )
    reference = _reference_positions(
        points,
        view_center_alt_deg=view_center_alt_deg,
        view_center_az_deg=view_center_az_deg,
        elapsed_seconds=elapsed_seconds,
    )
    errors = np.linalg.norm(mesh_positions - reference, axis=1)
    return float(np.mean(errors)), float(np.max(errors))


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    if args.columns < 1 or args.rows < 1:
        print("columns and rows must be positive", file=sys.stderr)
        return 2
    points = _screen_grid(args.columns, args.rows)
    cases = (
        (45.0, 0.0, "north"),
        (45.0, 180.0, "south"),
        (80.0, 90.0, "zenith"),
        (10.0, 270.0, "low-altitude"),
    )
    maximum_error = 0.0
    print(f"mesh={args.columns}x{args.rows} samples={len(points)}")
    for alt_deg, az_deg, name in cases:
        for elapsed_seconds in DEFAULT_ELAPSED_SECONDS:
            mean, maximum = _measure_case(
                points,
                view_center_alt_deg=alt_deg,
                view_center_az_deg=az_deg,
                elapsed_seconds=elapsed_seconds,
                columns=args.columns,
                rows=args.rows,
            )
            maximum_error = max(maximum_error, maximum)
            print(
                f"{name:12s} elapsed={elapsed_seconds:+.0f}s "
                f"mean={mean:.4f}px max={maximum:.4f}px"
            )
    result = "GO" if maximum_error < args.max_error_px else "NOT GO"
    print(f"gate_max_error={maximum_error:.4f}px threshold={args.max_error_px:.4f}px result={result}")
    if args.fail_on_not_go and result == "NOT GO":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
