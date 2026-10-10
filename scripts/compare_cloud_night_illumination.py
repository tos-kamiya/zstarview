#!/usr/bin/env python3
"""Export a deterministic synthetic nighttime cloud lighting comparison.

Run with the project's Python environment; no downloads are performed.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np

from zstarview.cloud_night_lighting import moon_cloud_strength
from zstarview.render.cloud_voxels import _altaz_to_directions, _render

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    x, y = np.meshgrid(np.arange(32) + 0.5, np.arange(32) + 0.5, indexing="ij")
    amount = np.clip(0.25 + 0.2 * np.sin(x * 0.6) * np.cos(y * 0.5), 0, 1)
    density = amount[..., None] * np.array([0.4, 0.2, 0.4])
    city = np.exp(-((x - 23) ** 2 + (y - 8) ** 2) / 30)
    ground = np.repeat(city[..., None], 3, axis=2)
    alt, az = np.meshgrid(np.linspace(65, 12, 120), np.linspace(-50, 50, 160), indexing="ij")
    rays = _altaz_to_directions(alt.ravel(), az.ravel())
    moon = _altaz_to_directions([45], [120])[0]
    basis = np.repeat(np.eye(2)[None], 3, axis=0)
    figure, axes = plt.subplots(1, 4, figsize=(12, 3), facecolor="#101018")
    panels = [
        ("Ambient only", 0, ground * 0),
        ("Full Moon", moon_cloud_strength(45, 180), ground * 0),
        ("Ground light", 0, ground),
        ("Moon + ground", moon_cloud_strength(45, 180), ground),
    ]
    for axis, (title, strength, surface) in zip(axes, panels):
        rgb, _ = _render(
            density, np.array([16.0, -8.0, -2.5]), rays, moon, 0.0,
            np.array([0.36, 0.385, 0.44]),
            np.tile([0.005, 0.009, 0.018], (len(rays), 1)), 1.0,
            np.zeros(2), basis, np.zeros((3, 2)), basis, moon, strength, surface,
        )
        axis.imshow(np.clip(rgb.reshape(120, 160, 3), 0, 1))
        axis.set_title(title, color="white")
        axis.axis("off")
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=150, facecolor=figure.get_facecolor())
    plt.close(figure)


if __name__ == "__main__":
    main()
