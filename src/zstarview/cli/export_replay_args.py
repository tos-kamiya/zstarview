"""Arguments for read-only replay cache coverage queries."""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections.abc import Sequence
from datetime import datetime, timedelta
from pathlib import Path
from typing import NoReturn

from ..paths import OVERLAY_FONT_SIZE_DEFAULT, WINDOW_HEIGHT, WINDOW_WIDTH
from ..replay.models import UTC, parse_utc
from .args import (
    _normalize_diffuse_sky_arguments,
    _normalize_location_arguments,
    _normalize_cloud_stripe_default,
    _normalize_vmag_limit,
    _validate_cloud_stripe_for_mode,
    _parse_image_size,
    _parse_theme,
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


_RELATIVE_TIME = re.compile(r"-(\d+):(\d{2})\Z")


def _parse_replay_time(raw: str | None, reference: datetime) -> datetime | None:
    if raw is None:
        return None
    match = _RELATIVE_TIME.fullmatch(raw)
    if match:
        hours, minutes = int(match.group(1)), int(match.group(2))
        if minutes > 59:
            raise ValueError("Relative time minutes must be between 00 and 59")
        return reference - timedelta(hours=hours, minutes=minutes)
    parsed = parse_utc(raw)
    if parsed.second or parsed.microsecond:
        raise ValueError("Time must be on a minute boundary")
    return parsed


def _parse_hud_items(raw: str) -> tuple[str, ...]:
    supported = {
        "datetime",
        "location",
        "coordinates",
        "height",
        "view-direction",
        "time-of-day",
    }
    if raw == "none":
        return ()
    if raw == "all":
        return tuple(sorted(supported))
    items = tuple(part.strip() for part in raw.split(","))
    if not items or any(not item for item in items):
        raise argparse.ArgumentTypeError("HUD items cannot be empty")
    unknown = sorted(set(items) - supported)
    if unknown:
        raise argparse.ArgumentTypeError(
            "Unsupported HUD item(s): " + ", ".join(unknown)
        )
    return tuple(dict.fromkeys(items))


def parse_replay_args(
    argv: Sequence[str] | None = None,
    *,
    reference_time: datetime | None = None,
) -> argparse.Namespace:
    if reference_time is None:
        reference_time = datetime.now(UTC)
    if reference_time.tzinfo is None:
        raise ValueError("Reference time must include a timezone")
    reference_time = reference_time.astimezone(UTC).replace(second=0, microsecond=0)
    parser = ReplayParser(
        prog="zstarview-export-replay",
        description="List local cache coverage or export a minute-stepped replay image sequence.",
    )
    add_observing_arguments(parser, include_time_arguments=False)
    add_sky_and_star_arguments(
        parser, include_sky_update_interval=False, include_twinkle_argument=False
    )
    add_overlay_arguments(parser)
    add_geo_satellite_argument(parser)
    parser.add_argument(
        "-t",
        "--theme",
        type=_parse_theme,
        default="night",
        help="Theme preset for background and star contrast (default: night).",
    )
    parser.add_argument("--visibility-boost", type=float, default=1.0)
    parser.add_argument(
        "--overlay-font-size", type=float, default=float(OVERLAY_FONT_SIZE_DEFAULT)
    )
    parser.set_defaults(
        precipitation_opacity=0.0,
        tropical_cyclone_opacity=0.0,
        meteor_trails_opacity=0.0,
        terrain_horizon_opacity=0.003,
        urban_outline_opacity=0.2,
        water_surface_opacity=0.0,
        night_light_opacity=0.22,
        road_light_opacity=0.0,
        ridge_glow_opacity=0.08,
        diffuse_sky_opacity=0.0,
        akari_ir_bands_opacity=0.0,
    )
    for action in parser._actions:
        if action.dest in {
            "precipitation_opacity",
            "tropical_cyclone_opacity",
            "meteor_trails_opacity",
            "water_surface_opacity",
            "road_light_opacity",
            "diffuse_sky_opacity",
            "akari_ir_bands_opacity",
        }:
            action.default = 0.0
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
    parser.add_argument(
        "--start",
        help="Inclusive zoned ISO 8601 minute or relative offset such as -08:30.",
    )
    parser.add_argument(
        "--end",
        help="Exclusive zoned ISO 8601 minute or relative offset such as -00:00.",
    )
    parser.add_argument(
        "--image-size", type=_parse_image_size, default=(WINDOW_WIDTH, WINDOW_HEIGHT)
    )
    parser.add_argument(
        "--jobs",
        type=int,
        choices=range(1, 9),
        default=2,
        metavar="N",
        help="Maximum parallel cloud-render subprocesses (1-8, default: 2).",
    )
    parser.add_argument(
        "-o", "--output-dir", type=Path, help="Directory for the replay image job."
    )
    parser.add_argument(
        "--plan-only",
        action="store_true",
        help="Validate and report without rendering.",
    )
    parser.add_argument(
        "--allow-partial-data",
        action="store_true",
        help="Export frames with unavailable layers marked in metadata.",
    )
    parser.add_argument(
        "--resume", action="store_true", help="Resume a matching existing image job."
    )
    parser.add_argument(
        "--include-direction-grid",
        action="store_true",
        help="Include the direction grid in each exported frame.",
    )
    parser.add_argument(
        "--hud",
        type=_parse_hud_items,
        default=(),
        metavar="ITEMS",
        help="HUD fields: datetime,location,coordinates,height,view-direction,time-of-day; default: none.",
    )
    parser.add_argument(
        "--observation-info",
        choices=("top", "bottom"),
        default="bottom",
        help="Position for selected HUD fields (default: bottom).",
    )
    args = parser.parse_args(argv)
    _normalize_cloud_stripe_default(
        parser, args, list(argv) if argv is not None else sys.argv[1:]
    )
    args.cloud_voxel = args.cloud_mode == "voxel"
    _validate_cloud_stripe_for_mode(parser, args)
    _normalize_diffuse_sky_arguments(
        parser, args, list(argv) if argv is not None else sys.argv[1:]
    )
    _normalize_vmag_limit(args)
    if args.list_cache_ranges and args.output_dir:
        parser.error("--output-dir cannot be used with --list-cache-ranges")
    if not args.list_cache_ranges and not args.plan_only and not args.output_dir:
        parser.error("--output-dir is required for image export")
    if args.list_cache_ranges and (
        args.plan_only or args.allow_partial_data or args.resume
    ):
        parser.error(
            "--plan-only, --allow-partial-data, and --resume require image export"
        )
    if args.json and not args.list_cache_ranges and not args.plan_only:
        parser.error("--json is only supported with --list-cache-ranges or --plan-only")
    if args.observation_info not in {"top", "bottom"}:
        parser.error("--observation-info must be top or bottom in replay")
    if args.resume and args.plan_only:
        parser.error("--resume cannot be combined with --plan-only")
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
            parsed = _parse_replay_time(raw, reference_time)
        except ValueError as exc:
            parser.error(str(exc))
        setattr(args, name, parsed)
    if args.start and args.end and args.start >= args.end:
        parser.error("--start must be earlier than --end")
    if not args.list_cache_ranges and (args.start is None or args.end is None):
        parser.error("--start and --end are required for image export")
    args.time_reference_utc = reference_time
    return args
