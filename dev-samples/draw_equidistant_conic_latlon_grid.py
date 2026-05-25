#!/usr/bin/env python3
"""Overlay a lat/lon grid on an image using an Equidistant Conic lon/lat map."""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image

from zstarview.utils.geostationary_grid_overlay import (
    draw_geostationary_latlon_grid,
    load_lonlat_grid_from_npz,
)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Overlay a 10-degree lat/lon grid using an Equidistant Conic lon/lat grid.",
    )
    parser.add_argument("image", type=Path, help="Input image path.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output image path. Defaults to <input>_grid.png.",
    )
    parser.add_argument(
        "--grid-npz",
        type=Path,
        required=True,
        help="Precomputed lon/lat grid saved as .npz, typically from fit_equidistant_conic_image_mapping.py.",
    )
    parser.add_argument(
        "--step-deg",
        type=int,
        default=10,
        help="Grid spacing in degrees for both latitude and longitude lines.",
    )
    parser.add_argument(
        "--major-step-deg",
        type=int,
        default=30,
        help="Draw multiples of this spacing in red instead of black.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    image_path = args.image.expanduser()
    if not image_path.exists():
        parser.error(f"Image file does not exist: {image_path}")

    grid_npz = args.grid_npz.expanduser()
    if not grid_npz.exists():
        parser.error(f"Grid file does not exist: {grid_npz}")

    output_path = (
        args.output.expanduser()
        if args.output is not None
        else image_path.with_name(f"{image_path.stem}_grid{image_path.suffix}")
    )

    grid = load_lonlat_grid_from_npz(grid_npz)
    with Image.open(image_path) as image:
        overlay = draw_geostationary_latlon_grid(
            image,
            grid.lon_deg,
            grid.lat_deg,
            step_deg=int(args.step_deg),
            major_step_deg=int(args.major_step_deg),
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output_path)
    print(f"Wrote: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
