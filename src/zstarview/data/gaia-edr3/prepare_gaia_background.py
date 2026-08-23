"""Prepare equirectangular Gaia sky textures at several display sizes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from PIL import Image


DEFAULT_SOURCE = Path(
    "The_colour_of_the_sky_from_Gaia_s_Early_Data_Release_3-copy.png"
)
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent


def parse_sizes(values: list[str]) -> list[tuple[int, int]]:
    sizes = []
    for value in values:
        try:
            width_text, height_text = value.lower().split("x", 1)
            width, height = int(width_text), int(height_text)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                f"invalid texture size: {value!r}; use WIDTHxHEIGHT"
            ) from exc
        if width <= 0 or height <= 0:
            raise argparse.ArgumentTypeError("texture dimensions must be positive")
        sizes.append((width, height))
    return sizes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", nargs="?", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--size",
        dest="sizes",
        action="append",
        default=None,
        help="output size WIDTHxHEIGHT; may be repeated",
    )
    args = parser.parse_args()

    requested_sizes = parse_sizes(args.sizes or ["2048x1024"])
    with Image.open(args.source) as source:
        image = source.convert("RGB")
        source_size = image.size

        args.output_dir.mkdir(parents=True, exist_ok=True)
        outputs = []
        for width, height in requested_sizes:
            if (width, height) == source_size:
                prepared = image.copy()
            else:
                prepared = image.resize((width, height), Image.Resampling.LANCZOS)
            output = args.output_dir / f"gaia-edr3-colour-{width}x{height}.png"
            prepared.save(output, format="PNG", optimize=True)
            outputs.append({"path": str(output), "width": width, "height": height})
            print(f"saved: {output} ({width}x{height})")

    manifest = {
        "source": str(args.source),
        "source_size": list(source_size),
        "projection": "equirectangular",
        "coordinate_frame": "Galactic coordinates",
        "content": "Gaia EDR3 integrated stellar brightness and colour",
        "credit": "ESA/Gaia/DPAC; acknowledgement: A. Moitinho",
        "license": "CC BY-SA 3.0 IGO or ESA Standard Licence",
        "outputs": outputs,
    }
    manifest_path = args.output_dir / "gaia-edr3-colour-manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    print(f"saved: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
