#!/usr/bin/env python3
"""Move or remove the old lower-left time-of-day marker in PNG images.

The marker geometry matches ``draw_instrument_time_of_day_marker``:
31 pixels on a side with a 4-pixel margin. The marker color is sampled from
each image, so images captured at different times can be processed together.

Examples:
    python scripts/edit_time_of_day_marker.py docs/images/scheduled
    python scripts/edit_time_of_day_marker.py docs/images/scheduled --remove
    python scripts/edit_time_of_day_marker.py docs/images/scheduled --in-place
    python scripts/edit_time_of_day_marker.py image.png --output-dir converted
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from PIL import Image

DEFAULT_MARKER_SIZE = 31
DEFAULT_MARGIN = 4
DEFAULT_COLOR_TOLERANCE = 2


def _triangle_pixels(
    *,
    height: int,
    size: int,
    margin: int,
) -> list[tuple[int, int]]:
    """Return pixels in the old lower-left right triangle."""
    left = margin
    top = height - margin - size
    bottom = height - margin - 1
    right = left + size
    pixels: list[tuple[int, int]] = []
    for y in range(max(0, top), min(height, bottom + 1)):
        diagonal_right = min(right, left + (y - top))
        pixels.extend((x, y) for x in range(left, diagonal_right + 1))
    return pixels


def _color_distance(first: tuple[int, ...], second: tuple[int, ...]) -> int:
    return max(abs(int(a) - int(b)) for a, b in zip(first[:3], second[:3]))


def edit_marker(
    image: Image.Image,
    *,
    remove: bool = False,
    marker_size: int = DEFAULT_MARKER_SIZE,
    margin: int = DEFAULT_MARGIN,
    color_tolerance: int = DEFAULT_COLOR_TOLERANCE,
) -> bool:
    """Move or remove the marker in ``image``; return whether it was found."""
    if image.mode != "RGBA":
        raise ValueError(f"expected RGBA image, got {image.mode}")
    width, height = image.size
    if width <= margin or height <= (margin * 2 + marker_size):
        raise ValueError(f"image is too small: {width}x{height}")

    pixels = image.load()
    left = margin
    top = height - margin - marker_size
    sample_x = left + max(1, marker_size // 3)
    sample_y = top + max(1, (marker_size * 2) // 3)
    sample = pixels[sample_x, sample_y]
    background = pixels[left - 1, sample_y]
    if _color_distance(sample, background) <= color_tolerance:
        return False

    source_pixels = _triangle_pixels(
        height=height,
        size=marker_size,
        margin=margin,
    )
    moved_pixels = {
        (x, height - 1 - y): pixels[x, y] for x, y in source_pixels
    }

    # The rounded sky viewport leaves the left margin outside the scene. Use
    # the adjacent pixel to restore the old marker area without hard-coding a
    # background color that may vary with the capture time.
    for x, y in source_pixels:
        pixels[x, y] = pixels[left - 1, y]
    if not remove:
        for (x, y), color in moved_pixels.items():
            pixels[x, y] = color
    return True


def _input_paths(values: list[Path]) -> list[Path]:
    paths: list[Path] = []
    for value in values:
        if value.is_dir():
            paths.extend(sorted(value.glob("*.png")))
        elif value.is_file():
            paths.append(value)
        else:
            raise FileNotFoundError(value)
    return paths


def _output_path(
    path: Path,
    output_dir: Path | None,
    in_place: bool,
    remove: bool,
) -> Path:
    if in_place:
        return path
    destination = output_dir or path.parent
    suffix = "-without-marker" if remove else "-top-left"
    return destination / f"{path.stem}{suffix}{path.suffix}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument(
        "--remove",
        action="store_true",
        help="remove the marker instead of moving it to the upper-left",
    )
    parser.add_argument("--in-place", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--marker-size", type=int, default=DEFAULT_MARKER_SIZE)
    parser.add_argument("--margin", type=int, default=DEFAULT_MARGIN)
    parser.add_argument("--color-tolerance", type=int, default=DEFAULT_COLOR_TOLERANCE)
    args = parser.parse_args()
    if args.in_place and args.output_dir is not None:
        parser.error("--in-place and --output-dir cannot be combined")

    paths = _input_paths(args.paths)
    if not paths:
        parser.error("no PNG files found")
    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)

    converted = 0
    skipped = 0
    for path in paths:
        source_stat = path.stat()
        with Image.open(path) as source:
            image = source.convert("RGBA")
        if edit_marker(
            image,
            remove=args.remove,
            marker_size=args.marker_size,
            margin=args.margin,
            color_tolerance=args.color_tolerance,
        ):
            destination = _output_path(
                path,
                args.output_dir,
                args.in_place,
                args.remove,
            )
            destination.parent.mkdir(parents=True, exist_ok=True)
            image.save(destination, format="PNG")
            os.utime(
                destination,
                ns=(source_stat.st_atime_ns, source_stat.st_mtime_ns),
            )
            converted += 1
            print(f"converted: {path} -> {destination}")
        else:
            skipped += 1
            print(f"skipped: {path} (marker not detected)")

    print(f"summary: converted={converted} skipped={skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
