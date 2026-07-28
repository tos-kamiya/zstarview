"""Console entry point for downloading and preparing AKARI dust maps."""

from __future__ import annotations

import argparse

from ..molecular_cloud_download import (
    AKARI_SOURCE_BASE_URL,
    DEFAULT_BANDS,
    DEFAULT_HEIGHT,
    DEFAULT_WIDTH,
    prepare_akari_data,
)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download AKARI far-infrared maps and prepare a display cache."
    )
    parser.add_argument(
        "--bands",
        default=",".join(DEFAULT_BANDS),
        help="comma-separated bands in micrometres (default: 90,140,160)",
    )
    parser.add_argument("--cache-dir", help="base directory for the molecular-cloud cache")
    parser.add_argument("--width", type=int, default=DEFAULT_WIDTH, help="output longitude samples")
    parser.add_argument("--height", type=int, default=DEFAULT_HEIGHT, help="output latitude samples")
    parser.add_argument("--source-base-url", default=AKARI_SOURCE_BASE_URL, help=argparse.SUPPRESS)
    parser.add_argument("--timeout", type=float, default=600.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)
    bands = tuple(part.strip() for part in args.bands.split(",") if part.strip())
    try:
        root = prepare_akari_data(
            bands=bands,
            cache_dir=args.cache_dir,
            source_base_url=args.source_base_url,
            width=args.width,
            height=args.height,
            timeout_s=args.timeout,
        )
    except Exception as exc:
        parser.exit(1, f"molecular-cloud data preparation failed: {exc}\n")
    print(f"Prepared molecular-cloud cache: {root}")
    return 0
