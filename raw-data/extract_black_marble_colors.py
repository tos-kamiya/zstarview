#!/usr/bin/env python3
"""Extract representative display colors from a 2012 Black Marble image.

This is a palette-extraction helper, not a radiometric analysis tool.  The
input JPEG is a display composite, so the reported colors are intended for
visual styling in zstarview.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image


def trimmed_mean(values: np.ndarray, trim: float = 0.05) -> np.ndarray:
    """Return a per-channel mean after dropping both tails."""
    ordered = np.sort(values.astype(np.float64), axis=0)
    count = len(ordered)
    margin = int(count * trim)
    if margin * 2 >= count:
        return np.mean(ordered, axis=0)
    return np.mean(ordered[margin : count - margin], axis=0)


def describe(name: str, pixels: np.ndarray) -> None:
    """Print robust RGB summaries for a selected pixel set."""
    median = np.median(pixels, axis=0)
    mean = trimmed_mean(pixels)
    representative = np.rint((median + mean) / 2).astype(int)
    normalized = np.rint(representative / representative.max() * 255).astype(int)
    print(f"{name} pixels: {len(pixels):,}")
    print(f"  median_rgb:         {np.rint(median).astype(int).tolist()}")
    print(f"  trimmed_mean_rgb:   {np.rint(mean).astype(int).tolist()}")
    print(f"  representative_rgb: {representative.tolist()}")
    print(f"  normalized_rgb:     {normalized.tolist()}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "image",
        nargs="?",
        type=Path,
        default=Path(__file__).with_name("BlackMarble_2012_01deg.jpg"),
    )
    args = parser.parse_args()

    with Image.open(args.image) as source:
        rgb = np.asarray(source.convert("RGB"), dtype=np.uint8)

    height, width, _ = rgb.shape
    latitude = 90.0 - (np.arange(height) + 0.5) * 180.0 / height
    lat = latitude[:, None]
    red = rgb[..., 0].astype(np.float32)
    green = rgb[..., 1].astype(np.float32)
    blue = rgb[..., 2].astype(np.float32)
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue

    # Exclude polar land from the city-light selection.  The warm/neutral
    # test rejects most of the blue-violet land background while retaining
    # yellow and white artificial lights.
    city_area = (lat > -55.0) & (lat < 55.0)
    city_color = (red >= blue * 1.03) & (green >= blue * 1.01)
    city_candidates = luminance[city_area & city_color]
    city_low_cutoff, city_high_cutoff = np.percentile(
        city_candidates, [97.0, 99.5]
    )
    city_color_mask = (
        city_area
        & city_color
        & (luminance >= city_low_cutoff)
        & (luminance < city_high_cutoff)
    )
    describe("city-light-color", rgb[city_color_mask])
    print(
        "  luminance_band:     "
        f"{city_low_cutoff:.2f} <= value < {city_high_cutoff:.2f}"
    )

    city_highlight_mask = city_area & city_color & (luminance >= city_high_cutoff)
    describe("city-light-highlight", rgb[city_highlight_mask])

    # Sample the interior ice sheet, avoiding the coastline and the very
    # bottom edge where the map transitions to the image/background boundary.
    antarctica_area = (lat > -82.0) & (lat < -72.0)
    antarctica_mask = antarctica_area & (luminance > 8.0)
    describe("antarctica-edge-glow", rgb[antarctica_mask])


if __name__ == "__main__":
    main()
