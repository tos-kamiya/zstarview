"""Compare the reference and Numba alt/az cloud-cell renderers.

Run with:
    uv run -p .venv/bin/python scripts/benchmark_altaz_cloud_renderer.py
"""

from __future__ import annotations

import argparse
import datetime as dt
import time

import numpy as np

from zstarview.clouddisc.altaz_grid import CloudAltAzGrid
from zstarview.clouddisc.altaz_projection import altaz_to_screen_coords
from zstarview.clouddisc.altaz_render import (
    render_altaz_grid_circles,
    render_altaz_grid_circles_reference,
)
from zstarview.clouddisc.types import SourceKey


def _make_grid(active_fraction: float, seed: int) -> CloudAltAzGrid:
    rng = np.random.default_rng(seed)
    amount = np.zeros((90, 720), dtype=np.float32)
    active_count = int(amount.size * active_fraction)
    selected = rng.choice(amount.size, size=active_count, replace=False)
    amount.flat[selected] = rng.uniform(0.04, 1.0, size=active_count).astype(np.float32)
    missing = np.zeros(amount.shape, dtype=np.uint8)
    timestamp = dt.datetime(2026, 10, 2, 0, tzinfo=dt.timezone.utc)
    source_key = SourceKey(
        satellite="G19",
        provider="GOES",
        timeslot_utc=timestamp,
        sat_priority=("AUTO",),
    )
    return CloudAltAzGrid(
        amount=amount,
        missing_mask=missing,
        alt_min_deg=0.0,
        alt_max_deg=90.0,
        az_min_deg=0.0,
        az_max_deg=360.0,
        observer_lat=35.0,
        observer_lon=135.0,
        satellite="G19",
        product="CMIPF-C13",
        time_utc=timestamp,
        shells_km=(6374.0, 6376.0, 6378.0),
        source_key=source_key,
        coverage_ratio=1.0,
    )


def _render(renderer, grid: CloudAltAzGrid, size: int) -> np.ndarray:
    return renderer(
        grid,
        size,
        size,
        center_alt_deg=45.0,
        center_az_deg=180.0,
        edge_fov_deg=117.0,
        mask_fov_deg=117.0,
    )


def _visible_active_count(grid: CloudAltAzGrid, size: int) -> int:
    alt_idx, az_idx = np.nonzero(grid.amount > 0.03)
    alt = grid.alt_min_deg + (alt_idx + 0.5) * (
        (grid.alt_max_deg - grid.alt_min_deg) / grid.amount.shape[0]
    )
    az = grid.az_min_deg + (az_idx + 0.5) * (
        (grid.az_max_deg - grid.az_min_deg) / grid.amount.shape[1]
    )
    x_px, y_px = altaz_to_screen_coords(
        alt,
        az,
        width=size,
        height=size,
        center_alt_deg=45.0,
        center_az_deg=180.0,
        edge_fov_deg=117.0,
        mask_fov_deg=117.0,
        observer_lat_deg=grid.observer_lat,
        observer_lon_deg=grid.observer_lon,
    )
    return int(np.count_nonzero(np.isfinite(x_px) & np.isfinite(y_px)))


def _timings(renderer, grid: CloudAltAzGrid, size: int, repeats: int) -> list[float]:
    samples = []
    for _ in range(repeats):
        started = time.perf_counter()
        _render(renderer, grid, size)
        samples.append(time.perf_counter() - started)
    return samples


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=600, help="square output size")
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()

    print(
        f"Output: {args.size}x{args.size}; warm repetitions: {args.repeats}; "
        "FOV: 117 deg; grid: 90x720"
    )
    print(
        "density  active/visible  reference median/p95  numba median/p95  "
        "speedup  cold numba  max alpha delta  identical"
    )

    for case_index, density in enumerate((0.02, 0.10, 0.35)):
        grid = _make_grid(density, seed=20261002 + case_index)
        visible_count = _visible_active_count(grid, args.size)
        reference = _render(render_altaz_grid_circles_reference, grid, args.size)

        cold_started = time.perf_counter()
        accelerated = _render(render_altaz_grid_circles, grid, args.size)
        cold_numba = time.perf_counter() - cold_started
        alpha_delta = np.abs(
            reference[..., 3].astype(np.int16) - accelerated[..., 3].astype(np.int16)
        )

        reference_times = _timings(
            render_altaz_grid_circles_reference, grid, args.size, args.repeats
        )
        numba_times = _timings(render_altaz_grid_circles, grid, args.size, args.repeats)
        reference_median = float(np.median(reference_times))
        numba_median = float(np.median(numba_times))
        speedup = reference_median / numba_median if numba_median > 0 else float("inf")
        reference_p95 = float(np.percentile(reference_times, 95))
        numba_p95 = float(np.percentile(numba_times, 95))
        print(
            f"{density:5.0%}  "
            f"{int(np.count_nonzero(grid.amount > 0.03)):5d}/{visible_count:5d}  "
            f"{reference_median:7.4f}/{reference_p95:7.4f}s  "
            f"{numba_median:7.4f}/{numba_p95:7.4f}s  "
            f"{speedup:6.2f}x  {cold_numba:7.4f}s  {alpha_delta.max():d}  "
            f"{np.array_equal(reference, accelerated)}"
        )


if __name__ == "__main__":
    main()
