#!/usr/bin/env python3
"""Export the prepared molecular-cloud grid as inspection images."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image

from zstarview.render.molecular_cloud_overlay import MOLECULAR_CLOUD_CACHE


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Export the prepared AKARI molecular-cloud grid without reprojection. "
            "The image axes remain Galactic longitude and latitude."
        ),
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=MOLECULAR_CLOUD_CACHE,
        help="Prepared NPZ path (default: the zstarview molecular-cloud cache).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("build/molecular-cloud-images"),
        help="Directory for exported images (default: build/molecular-cloud-images).",
    )
    parser.add_argument(
        "--format",
        choices=("png", "jpeg"),
        default="png",
        help="Image format (default: png).",
    )
    return parser


def _load_asset(path: Path) -> tuple[np.ndarray, tuple[int, ...]]:
    with np.load(path) as archive:
        data = np.asarray(archive["data"])
        bands = tuple(int(value) for value in np.asarray(archive["bands"]).tolist())
    if data.ndim != 3 or data.shape[0] != len(bands):
        raise ValueError(f"unexpected data shape or band list in {path}")
    if data.dtype != np.uint16:
        data = np.clip(data, 0, 65535).astype(np.uint16)
    return data, bands


def _to_u8(channel: np.ndarray) -> np.ndarray:
    return np.right_shift(channel.astype(np.uint16), 8).astype(np.uint8)


def _save_gray(channel: np.ndarray, path: Path, image_format: str) -> None:
    if image_format == "png":
        Image.fromarray(channel.astype(np.uint16)).save(path, format="PNG")
    else:
        Image.fromarray(_to_u8(channel), mode="L").save(path, format="JPEG", quality=95)


def _save_rgb(data: np.ndarray, bands: tuple[int, ...], path: Path, image_format: str) -> None:
    channels = {band: data[index] for index, band in enumerate(bands)}
    red = channels.get(160, np.zeros(data.shape[1:], dtype=np.uint16))
    green = channels.get(140, red)
    blue = channels.get(90, green)
    rgb = np.stack((_to_u8(red), _to_u8(green), _to_u8(blue)), axis=-1)
    image = Image.fromarray(rgb, mode="RGB")
    if image_format == "jpeg":
        image.save(path, format="JPEG", quality=95, subsampling=0)
    else:
        image.save(path, format="PNG")


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)
    input_path = args.input.expanduser()
    if not input_path.is_file():
        parser.error(f"Input file does not exist: {input_path}")

    try:
        data, bands = _load_asset(input_path)
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))

    output_dir = args.output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = str(args.format)
    _save_rgb(data, bands, output_dir / f"molecular-cloud-rgb.{suffix}", suffix)
    for index, band in enumerate(bands):
        _save_gray(data[index], output_dir / f"molecular-cloud-{band}um.{suffix}", suffix)

    height, width = data.shape[1:]
    print(f"Input: {input_path}")
    print(f"Grid: {width} x {height} (Galactic longitude x latitude)")
    print(f"Bands: {', '.join(f'{band} um' for band in bands)}")
    print(f"Output: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
