#!/usr/bin/env python3
"""Overlay a geostationary lat/lon grid on top of an image."""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image

from zstarview.utils.geostationary_grid_overlay import (
    draw_geostationary_latlon_grid,
    load_lonlat_grid_for_image,
)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Overlay a 10-degree lat/lon grid on a geostationary image.",
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
        default=None,
        help="Use a precomputed lon/lat grid saved as .npz.",
    )
    parser.add_argument(
        "--control-map",
        type=Path,
        default=None,
        help="Use a latlonmap.txt-style control map to fit the image projection.",
    )
    parser.add_argument(
        "--lon0",
        type=float,
        default=0.0,
        help="Longitude of projection origin in degrees when fitting from control points.",
    )
    parser.add_argument(
        "--height-m",
        type=float,
        default=35785831.0,
        help="Perspective point height in meters when fitting from control points.",
    )
    parser.add_argument(
        "--sweep-axis",
        choices=("x", "y"),
        default="x",
        help="Geostationary sweep axis when fitting from control points.",
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


def _resolve_companion_map(image_path: Path) -> Path | None:
    candidates = [
        image_path.with_name(f"{image_path.stem}_lonlat.npz"),
        image_path.with_name(f"{image_path.stem}_grid.npz"),
        image_path.with_name(f"{image_path.stem}_latlonmap.txt"),
        image_path.with_name("latlonmap.txt"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    image_path = args.image.expanduser()
    if not image_path.exists():
        parser.error(f"Image file does not exist: {image_path}")

    output_path = (
        args.output.expanduser()
        if args.output is not None
        else image_path.with_name(f"{image_path.stem}_grid{image_path.suffix}")
    )

    companion = _resolve_companion_map(image_path)
    grid_npz = args.grid_npz.expanduser() if args.grid_npz is not None else (companion if companion is not None and companion.suffix == ".npz" else None)
    control_map = args.control_map.expanduser() if args.control_map is not None else (
        companion if companion is not None and companion.suffix == ".txt" else None
    )
    if grid_npz is None and control_map is None:
        parser.error(
            "No mapping source found. Provide --grid-npz or --control-map, or place a companion "
            "latlonmap.txt / *_lonlat.npz file next to the image."
        )

    grid = load_lonlat_grid_for_image(
        image_path,
        grid_npz=grid_npz,
        control_map=control_map,
        longitude_of_projection_origin=float(args.lon0),
        perspective_point_height=float(args.height_m),
        sweep_angle_axis=str(args.sweep_axis),
    )

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
