"""Read-only cache range queries for future replay image export."""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Any

from ..location_resolver.resolve import LocationResolveError
from ..paths import CACHE_PATH
from ..replay.dependencies import check_dependencies
from ..replay.inventory import Inventory
from ..replay.location import location_payload, resolve_location
from ..replay.report import build_report, format_report
from ..replay.static_sources import static_source_checks
from .export_replay_args import parse_replay_args


def query_cache_ranges(args: Any) -> dict[str, Any]:
    # Import services only when doing native work, never create a Qt app.
    from ..gui.application_services import ApplicationServices

    root = Path(args.cache_dir or CACHE_PATH).expanduser().resolve()
    services = ApplicationServices(max_workers=1)
    try:
        with services.native_work_lock:
            location, issues = resolve_location(args, root)
            height = (
                None
                if any(issue.startswith("building-") for issue in issues)
                else location.observer_height_m
            )
            inventory = Inventory(root, location.lat, location.lon, height)
            inventory.scan()
            enabled = set()
            if args.cloud_opacity > 0:
                enabled.add("cloud-geo" if args.geo_satellite else "cloud-native")
            if args.aircraft_opacity > 0:
                enabled.add("aircraft")
            if args.satellite_opacity > 0:
                enabled.update({"satellite-iss", "satellite-horizons"})
            checks, ephemeris = check_dependencies(root, args, issues)
            checks.extend(static_source_checks(root, args, location))
            report = build_report(
                inventory,
                enabled,
                start=args.start,
                end=args.end,
                ephemeris_ranges=ephemeris,
                dependencies=checks,
            )
            report["location"] = location_payload(location)
            if "ground-elevation-unverified" in issues:
                report["location"]["ground_elevation_m"] = None
            if "building-top-unverified" in issues:
                report["location"]["location_height_m"] = None
                report["location"]["observer_height_m"] = None
            report["settings"] = {
                key: value
                for key, value in vars(args).items()
                if key.endswith("opacity")
                or key
                in {
                    "view_center_alt",
                    "view_center_az",
                    "image_size",
                    "edge_fov_deg",
                    "content_fov_deg",
                    "geo_satellite",
                }
            }
            return report
    finally:
        services.shutdown()


def main() -> None:
    args = parse_replay_args()

    # Place names can be Unicode; keep console/log output ASCII on all systems.
    class AsciiFormatter(logging.Formatter):
        def format(self, record: logging.LogRecord) -> str:
            return (
                super()
                .format(record)
                .encode("ascii", "backslashreplace")
                .decode("ascii")
            )

    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(AsciiFormatter("%(levelname)s: %(message)s"))
    logger = logging.getLogger("zstarview")
    logger.addHandler(handler)
    try:
        report = query_cache_ranges(args)
        print(
            json.dumps(report, ensure_ascii=True, indent=2)
            if args.json
            else format_report(report)
        )
    except (OSError, ValueError, RuntimeError, LocationResolveError) as exc:
        print(
            (f"Cache range query failed: {exc}")
            .encode("ascii", "backslashreplace")
            .decode("ascii"),
            file=sys.stderr,
        )
        raise SystemExit(1) from exc
    finally:
        logger.removeHandler(handler)


if __name__ == "__main__":
    main()
