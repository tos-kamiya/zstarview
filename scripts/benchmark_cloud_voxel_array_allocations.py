"""Compare the production scalarized kernel with its array-based baseline.

Run with:
    uv run -p .venv/bin/python scripts/benchmark_cloud_voxel_array_allocations.py

The baseline keeps the current segment traversal and recreates the former
per-cell RGB array arithmetic. It is intentionally kept outside the package.
"""

from __future__ import annotations

import argparse
import time

import numpy as np
from numba import njit

from zstarview.render.cloud_voxels import (
    CLOUD_WHITENESS,
    SUNLIGHT_LEVELS,
    _light_with_environment,
    _render as render_production,
    _segments,
)


@njit(cache=True)
def render_array_based(
    density,
    origin,
    rays,
    sun,
    sunlight_mix,
    night_color_rgb,
    base,
    opacity,
    pixel_origin,
    pixel_from_enu,
    layer_centers,
    layer_bases,
):
    """Baseline _render variant using per-cell RGB arrays."""
    shape = np.array(density.shape)
    light_fractions = np.full(density.shape, -1.0)
    night_amounts = np.minimum(density.sum(axis=2), 1.0)
    output = base.copy()
    transmissions = np.ones(len(rays), dtype=np.float32)
    for i in range(len(rays)):
        transmission = 1.0
        color = np.zeros(3)
        for x, y, z, distance in _segments(
            origin,
            rays[i],
            shape,
            pixel_origin,
            pixel_from_enu,
            layer_centers,
        ):
            amount = density[x, y, z]
            if amount <= 0:
                continue
            if sunlight_mix > 0.0 and light_fractions[x, y, z] < 0:
                pixel_delta = np.array(
                    [
                        x + 0.5 - pixel_origin[0],
                        y + 0.5 - pixel_origin[1],
                    ]
                )
                center_enu = np.empty(2)
                for axis in range(2):
                    center_enu[axis] = layer_centers[z, axis] + (
                        layer_bases[z, axis, 0] * pixel_delta[0]
                        + layer_bases[z, axis, 1] * pixel_delta[1]
                    )
                center = np.array([center_enu[0], center_enu[1], z + 0.5])
                tau = 0.0
                for sx, sy, sz, length in _segments(
                    center,
                    sun,
                    shape,
                    pixel_origin,
                    pixel_from_enu,
                    layer_centers,
                ):
                    tau += density[sx, sy, sz] * length * 1.8
                light_fractions[x, y, z] = _light_with_environment(1.0, tau)
            solar_color = np.zeros(3)
            if sunlight_mix > 0.0:
                brightness = np.interp(
                    abs(sun[2]) * light_fractions[x, y, z],
                    SUNLIGHT_LEVELS,
                    CLOUD_WHITENESS,
                )
                solar_color = brightness * np.array([0.96, 0.975, 1.0])
            alpha = 1.0 - np.exp(-2.4 * amount * distance * opacity)
            color_rgb = (
                sunlight_mix * solar_color
                + (1.0 - sunlight_mix) * night_color_rgb * night_amounts[x, y]
            )
            color += transmission * alpha * color_rgb
            transmission *= 1.0 - alpha
            if transmission < 1e-5:
                break
        output[i] = color + transmission * base[i]
        transmissions[i] = transmission
    return output, transmissions


def _inputs(ray_count: int, seed: int):
    rng = np.random.default_rng(seed)
    density = np.zeros((64, 64, 9), dtype=np.float64)
    occupied = rng.random(density.shape) < 0.38
    density[occupied] = rng.uniform(0.08, 0.7, size=np.count_nonzero(occupied))

    altitude = rng.uniform(8.0, 88.0, size=ray_count)
    azimuth = rng.uniform(0.0, 360.0, size=ray_count)
    alt_rad = np.radians(altitude)
    az_rad = np.radians(azimuth)
    rays = np.column_stack(
        (
            np.cos(alt_rad) * np.sin(az_rad),
            np.cos(alt_rad) * np.cos(az_rad),
            np.sin(alt_rad),
        )
    )
    sun_alt = np.radians(28.0)
    sun_az = np.radians(245.0)
    sun = np.array(
        [
            np.cos(sun_alt) * np.sin(sun_az),
            np.cos(sun_alt) * np.cos(sun_az),
            np.sin(sun_alt),
        ]
    )
    basis = np.repeat((np.eye(2) * 4.0)[None, :, :], 9, axis=0)
    pixel_from_enu = np.repeat((np.eye(2) / 4.0)[None, :, :], 9, axis=0)
    layer_centers = np.zeros((9, 2), dtype=np.float64)
    base = rng.uniform(0.0, 0.4, size=(ray_count, 3)).astype(np.float32)
    return (
        density,
        np.array([0.0, 0.0, -2.5]),
        rays,
        sun,
        0.85,
        np.array([0.36, 0.385, 0.44]),
        base,
        0.85,
        np.array([32.0, 32.0]),
        pixel_from_enu,
        layer_centers,
        basis,
    )


def _time_calls(renderer, inputs, repeats: int) -> tuple[float, float, list[float]]:
    started = time.perf_counter()
    renderer(*inputs)
    first_call = time.perf_counter() - started

    started = time.perf_counter()
    renderer(*inputs)
    second_call = time.perf_counter() - started

    samples = []
    for _ in range(repeats):
        started = time.perf_counter()
        renderer(*inputs)
        samples.append(time.perf_counter() - started)
    return first_call, second_call, samples


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rays", type=int, default=2048)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    inputs = _inputs(args.rays, seed=20261007)

    array_first, array_second, array_samples = _time_calls(
        render_array_based, inputs, args.repeats
    )
    array_result = render_array_based(*inputs)
    production_first, production_second, production_samples = _time_calls(
        render_production, inputs, args.repeats
    )
    production_result = render_production(*inputs)
    max_rgb_diff = float(np.max(np.abs(array_result[0] - production_result[0])))
    max_transmission_diff = float(
        np.max(np.abs(array_result[1] - production_result[1]))
    )
    array_median = float(np.median(array_samples))
    production_median = float(np.median(production_samples))

    print(
        f"Synthetic workload: density=64x64x9; rays={args.rays}; "
        f"sunlight_mix=0.85; mode=voxel; overlay=off; repeats={args.repeats}"
    )
    print("kernel             first-call    second-call   warm-median    warm-p95")
    print(
        f"array baseline     {array_first:9.4f}s  {array_second:9.4f}s  "
        f"{array_median:9.4f}s  {np.percentile(array_samples, 95):9.4f}s"
    )
    print(
        f"production         {production_first:9.4f}s  {production_second:9.4f}s  "
        f"{production_median:9.4f}s  "
        f"{np.percentile(production_samples, 95):9.4f}s"
    )
    print(
        f"warm speed ratio: array/production = {array_median / production_median:.3f}x"
    )
    print(
        f"max output delta: rgb={max_rgb_diff:.3g}; "
        f"transmission={max_transmission_diff:.3g}"
    )


if __name__ == "__main__":
    main()
