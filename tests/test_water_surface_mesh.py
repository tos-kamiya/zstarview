from __future__ import annotations

from zstarview.water_surface_mesh import split_local_polygon_by_grid


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


def test_split_local_polygon_by_grid_splits_shells() -> None:
    shell = [
        (0.0, 0.0),
        (900.0, 0.0),
        (900.0, 900.0),
        (0.0, 900.0),
        (0.0, 0.0),
    ]

    pieces = split_local_polygon_by_grid(shell, [], cell_m=300.0)

    assert len(pieces) > 1
    expected_area = 900.0 * 900.0
    got_area = sum(_ring_area(piece_shell) for piece_shell, _ in pieces)
    assert abs(got_area - expected_area) / expected_area < 0.01


def test_split_local_polygon_by_grid_keeps_hole_assignment() -> None:
    shell = [
        (0.0, 0.0),
        (900.0, 0.0),
        (900.0, 900.0),
        (0.0, 900.0),
        (0.0, 0.0),
    ]
    hole = [
        (300.0, 300.0),
        (600.0, 300.0),
        (600.0, 600.0),
        (300.0, 600.0),
        (300.0, 300.0),
    ]

    pieces = split_local_polygon_by_grid(shell, [hole], cell_m=300.0)

    assert pieces
    assert any(piece_holes for _, piece_holes in pieces)
