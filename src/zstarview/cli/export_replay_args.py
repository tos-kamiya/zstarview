"""Arguments for read-only replay cache coverage queries."""

from __future__ import annotations

import argparse
import math
import sys
from collections.abc import Sequence
from typing import NoReturn

from ..paths import WINDOW_HEIGHT, WINDOW_WIDTH
from ..replay.models import parse_utc
from .args import (
    _normalize_diffuse_sky_arguments,
    _normalize_location_arguments,
    _parse_image_size,
    _validate_fov_relationship,
    _validate_location_argument_combinations,
    _validate_urban_outline_argument_combinations,
    add_geo_satellite_argument,
    add_observing_arguments,
    add_overlay_arguments,
    add_sky_and_star_arguments,
)


class ReplayParser(argparse.ArgumentParser):
    def error(self, message: str) -> NoReturn:
        self.print_usage(sys.stderr)
        text = f"{self.prog}: error: {message}\n"
        self.exit(1, text.encode("ascii", "backslashreplace").decode("ascii"))


def parse_replay_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = ReplayParser(
        prog="zstarview-export-replay",
        description="List local cache coverage for offline replay. Image export is not implemented yet.",
    )
    add_observing_arguments(parser, include_time_arguments=False)
    add_sky_and_star_arguments(
        parser, include_sky_update_interval=False, include_twinkle_argument=False
    )
    add_overlay_arguments(parser)
    add_geo_satellite_argument(parser)
    parser.set_defaults(
        precipitation_opacity=0.0,
        tropical_cyclone_opacity=0.0,
        meteor_trails_opacity=0.0,
    )
    for action in parser._actions:
        if action.dest in {
            "precipitation_opacity",
            "tropical_cyclone_opacity",
            "meteor_trails_opacity",
        }:
            action.help = "Unsupported in replay; only 0 is accepted (default: 0)."
    parser.add_argument(
        "--list-cache-ranges",
        action="store_true",
        help="List cache timestamps and contiguous replay minute ranges.",
    )
    parser.add_argument(
        "--json", action="store_true", help="Write a versioned JSON report to stdout."
    )
    parser.add_argument(
        "--cache-dir", help="Input cache root (default: normal application cache)."
    )
    parser.add_argument("--start", help="Inclusive ISO 8601 minute with timezone.")
    parser.add_argument("--end", help="Exclusive ISO 8601 minute with timezone.")
    parser.add_argument(
        "--image-size", type=_parse_image_size, default=(WINDOW_WIDTH, WINDOW_HEIGHT)
    )
    args = parser.parse_args(argv)
    _normalize_diffuse_sky_arguments(
        parser, args, list(argv) if argv is not None else sys.argv[1:]
    )
    if not args.list_cache_ranges:
        parser.error(
            "--list-cache-ranges is required; image export is not implemented yet"
        )
    _normalize_location_arguments(parser, args)
    _validate_location_argument_combinations(parser, args)
    _validate_fov_relationship(parser, args)
    _validate_urban_outline_argument_combinations(parser, args)
    if not all(
        math.isfinite(value) for value in (args.view_center_alt, args.view_center_az)
    ):
        parser.error("Viewing angles must be finite")
    for name, value in vars(args).items():
        if (
            name.endswith("opacity")
            and value is not None
            and (not math.isfinite(value) or not 0 <= value <= 1)
        ):
            parser.error(f"--{name.replace('_', '-')} must be between 0 and 1")
    for name in (
        "precipitation_opacity",
        "tropical_cyclone_opacity",
        "meteor_trails_opacity",
    ):
        if getattr(args, name, 0) > 0:
            parser.error(f"--{name.replace('_', '-')} is unsupported in replay")
    for name in ("start", "end"):
        raw = getattr(args, name)
        try:
            parsed = parse_utc(raw) if raw else None
            if parsed and (parsed.second or parsed.microsecond):
                raise ValueError("Time must be on a minute boundary")
        except ValueError as exc:
            parser.error(str(exc))
        setattr(args, name, parsed)
    if args.start and args.end and args.start >= args.end:
        parser.error("--start must be earlier than --end")
    return args
