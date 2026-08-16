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
    parser.add_argument("--columns", type=int, default=16, help="mesh cell count across")
    parser.add_argument("--rows", type=int, default=9, help="mesh cell count down")
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
    parser.add_argument(
        "--coverage",
        type=float,
        default=1.0,
        help="fraction of elapsed sidereal motion to apply (default: 1.0)",
    )
    return parser.parse_args(argv)


def _mesh_edges(cell_size: float, phase: float, extent: float) -> np.ndarray:
    edges = np.arange(-phase * cell_size, extent + cell_size, cell_size)
    edges = np.unique(np.clip(np.append(edges, (0.0, extent)), 0.0, extent))
    return edges


def _screen_grid(columns: int, rows: int, phase_x: float, phase_y: float) -> np.ndarray:
    x, y = np.meshgrid(
        np.linspace(0.04 * WIDTH_PX, 0.96 * WIDTH_PX, 25),
        np.linspace(0.04 * HEIGHT_PX, 0.96 * HEIGHT_PX, 15),
    )
    points = np.column_stack((x.ravel(), y.ravel()))
    cell_width = WIDTH_PX / columns
    cell_height = HEIGHT_PX / rows
    x_edges = _mesh_edges(cell_width, phase_x, WIDTH_PX)
    y_edges = _mesh_edges(cell_height, phase_y, HEIGHT_PX)
    cell_x = np.clip(np.searchsorted(x_edges, points[:, 0], side="right") - 1, 0, len(x_edges) - 2)
    cell_y = np.clip(np.searchsorted(y_edges, points[:, 1], side="right") - 1, 0, len(y_edges) - 2)
    corner_x = np.column_stack(
        (
            x_edges[cell_x], x_edges[cell_x + 1],
            x_edges[cell_x], x_edges[cell_x + 1],
        )
    )
    corner_y = np.column_stack(
        (
            y_edges[cell_y], y_edges[cell_y],
            y_edges[cell_y + 1], y_edges[cell_y + 1],
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
    coverage: float,
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
        SIDEREAL_DEG_PER_SECOND * elapsed_seconds * coverage
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


def _mesh_vertices(columns: int, rows: int, phase_x: float, phase_y: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = _mesh_edges(WIDTH_PX / columns, phase_x, WIDTH_PX)
    y = _mesh_edges(HEIGHT_PX / rows, phase_y, HEIGHT_PX)
    grid_x, grid_y = np.meshgrid(x, y)
    return np.column_stack((grid_x.ravel(), grid_y.ravel())), x, y


def _mesh_interpolate(
    points: np.ndarray,
    transformed_vertices: np.ndarray,
    *,
    columns: int,
    rows: int,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
) -> np.ndarray:
    output = np.empty_like(points)
    for index, (x, y) in enumerate(points):
        cell_x = min(len(x_edges) - 2, max(0, int(np.searchsorted(x_edges, x, side="right") - 1)))
        cell_y = min(len(y_edges) - 2, max(0, int(np.searchsorted(y_edges, y, side="right") - 1)))
        u = (x - x_edges[cell_x]) / (x_edges[cell_x + 1] - x_edges[cell_x])
        v = (y - y_edges[cell_y]) / (y_edges[cell_y + 1] - y_edges[cell_y])
        stride = len(x_edges)
        top_left = cell_y * stride + cell_x
        top_right = top_left + 1
        bottom_left = top_left + stride
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
    phase_x: float,
    phase_y: float,
    coverage: float,
) -> tuple[float, float]:
    vertices, x_edges, y_edges = _mesh_vertices(columns, rows, phase_x, phase_y)
    transformed_vertices = _reference_positions(
        vertices,
        view_center_alt_deg=view_center_alt_deg,
        view_center_az_deg=view_center_az_deg,
        elapsed_seconds=elapsed_seconds,
        coverage=coverage,
    )
    mesh_positions = _mesh_interpolate(
        points, transformed_vertices, columns=columns, rows=rows,
        x_edges=x_edges, y_edges=y_edges
    )
    reference = _reference_positions(
        points,
        view_center_alt_deg=view_center_alt_deg,
        view_center_az_deg=view_center_az_deg,
        elapsed_seconds=elapsed_seconds,
        coverage=coverage,
    )
    errors = np.linalg.norm(mesh_positions - reference, axis=1)
    return float(np.mean(errors)), float(np.max(errors))


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    if args.columns < 1 or args.rows < 1 or not 0.0 <= args.coverage <= 1.0:
        print("columns and rows must be positive", file=sys.stderr)
        return 2
    cases = (
        (45.0, 0.0, "north"),
        (45.0, 180.0, "south"),
        (80.0, 90.0, "zenith"),
        (10.0, 270.0, "low-altitude"),
    )
    maximum_error = 0.0
    for phase_x, phase_y in ((0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5)):
        points = _screen_grid(args.columns, args.rows, phase_x, phase_y)
        print(f"mesh={args.columns}x{args.rows} phase={phase_x:.1f},{phase_y:.1f} samples={len(points)}")
        for alt_deg, az_deg, name in cases:
            for elapsed_seconds in DEFAULT_ELAPSED_SECONDS:
                mean, maximum = _measure_case(
                    points,
                    view_center_alt_deg=alt_deg,
                    view_center_az_deg=az_deg,
                    elapsed_seconds=elapsed_seconds,
                    columns=args.columns,
                    rows=args.rows,
                    phase_x=phase_x,
                    phase_y=phase_y,
                    coverage=args.coverage,
                )
                maximum_error = max(maximum_error, maximum)
                print(
                    f"{name:12s} elapsed={elapsed_seconds:+.0f}s "
                    f"mean={mean:.4f}px max={maximum:.4f}px"
                )
    result = "GO" if maximum_error < args.max_error_px else "NOT GO"
    print(f"coverage={args.coverage:.2f}")
    print(f"gate_max_error={maximum_error:.4f}px threshold={args.max_error_px:.4f}px result={result}")
    if args.fail_on_not_go and result == "NOT GO":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
