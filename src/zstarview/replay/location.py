"""Resolve observing inputs without downloading elevation/building data."""

from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

from ..location_resolver.resolve import (
    LocationDataPolicy,
    ResolvedLocation,
    resolve_launch_location,
)


class LocalLocationData:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.issues: set[str] = set()

    def ground(self, lat: float, lon: float) -> float:
        # Coordinate resolution needs a placeholder, but listing must not read DEMs.
        self.issues.add("ground-elevation-unverified")
        return 0.0

    def building(self, lat: float, lon: float) -> float | None:
        self.issues.add("building-top-unverified")
        return None


def resolve_location(args: Any, root: Path) -> tuple[ResolvedLocation, list[str]]:
    local = LocalLocationData(root)
    city, place = args.city, args.place
    if city == "auto":
        raise ValueError(
            "Automatic IP location is unavailable in replay; specify a location"
        )
    location = resolve_launch_location(
        city,
        place_query=place,
        place_countrycode=args.place_countrycode,
        place_lang=args.place_lang,
        use_building_top=args.use_building_top,
        save_last_city_func=lambda _: None,
        data_policy=LocationDataPolicy(
            local.ground, local.building, root / "nominatim" / "schema-1"
        ),
    )
    if args.observer_height_m is not None:
        location = replace(
            location,
            height_add_m=args.observer_height_m,
            observer_height_m=location.location_height_m + args.observer_height_m,
        )
    return location, sorted(local.issues)


def location_payload(location: ResolvedLocation) -> dict[str, Any]:
    return asdict(location)
