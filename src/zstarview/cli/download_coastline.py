"""Console entry point for downloading coastline vector data."""

from __future__ import annotations

import argparse

from ..coastline_download import download_coastline_data


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download coastline vector data for a longitude range."
    )
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--all", action="store_true", help="download all longitude columns")
    selection.add_argument("--lon-min", type=float, help="minimum longitude in degrees")
    parser.add_argument("--lon-max", type=float, help="maximum longitude in degrees")
    parser.add_argument("--cache-dir", help="coastline cache base directory")
    parser.add_argument("--base-url", help=argparse.SUPPRESS)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--download-timeout", type=float, default=600.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)
    if not args.all and args.lon_max is None:
        parser.error("--lon-max is required with --lon-min")
    try:
        columns = download_coastline_data(
            lon_min=args.lon_min,
            lon_max=args.lon_max,
            all_columns=args.all,
            cache_dir=args.cache_dir,
            base_url=args.base_url
            or "https://github.com/tos-kamiya/zstarview/releases/download/coastline-data-20260725",
            timeout_s=args.timeout,
            download_timeout_s=args.download_timeout,
            status_callback=print,
        )
    except Exception as exc:
        parser.exit(1, f"coastline download failed: {exc}\n")
    print("Downloaded coastline columns: " + ", ".join(f"x{column:02d}" for column in columns))
    return 0
