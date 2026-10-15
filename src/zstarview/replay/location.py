"""Resolve observing inputs without downloading elevation/building data."""

from __future__ import annotations

import math
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

from ..location_resolver.place_cache import build_place_cache_key, load_place_cache
from ..location_resolver.resolve import (
    LocationDataPolicy,
    ResolvedLocation,
    resolve_launch_location,
)
from .inventory import json_paths
from .models import parse_utc


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


def _saved_place(
    root: Path, query: str, country: str | None, language: str
) -> tuple[float, float] | None:
    import json

    key = build_place_cache_key(query, country, language)
    cached = load_place_cache(key, cache_root=root / "nominatim" / "schema-1")
    if cached and cached.results:
        item = cached.results[0]
        return float(str(item["lat"])), float(str(item["lon"]))
    matches: set[tuple[float, float]] = set()
    for path in json_paths(root / "replay" / "v1" / "launches"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            if payload.get("schema_version") != 1:
                continue
            parse_utc(payload["started_at_utc"])
            raw = payload["location"]["input"]
            if (
                build_place_cache_key(
                    raw["query"], raw.get("countrycode"), raw["language"]
                )
                != key
            ):
                continue
            resolved = payload["location"]["resolved"]
            lat, lon = float(resolved["latitude_deg"]), float(resolved["longitude_deg"])
            if (
                math.isfinite(lat)
                and math.isfinite(lon)
                and -90 <= lat <= 90
                and -180 <= lon <= 180
            ):
                matches.add((lat, lon))
        except (ValueError, KeyError, TypeError, FileNotFoundError):
            continue
    return next(iter(matches)) if len(matches) == 1 else None


def resolve_location(args: Any, root: Path) -> tuple[ResolvedLocation, list[str]]:
    local = LocalLocationData(root)
    city, place = args.city, args.place
    if city == "auto":
        raise ValueError(
            "Automatic IP location is unavailable in replay; specify a location"
        )
    saved = (
        _saved_place(root, place, args.place_countrycode, args.place_lang)
        if place
        else None
    )
    if saved:
        city, place = f"{saved[0]},{saved[1]}", None
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
    if saved:
        location = replace(location, display_name=args.place, kind="place-cache")
    if args.observer_height_m is not None:
        location = replace(
            location,
            height_add_m=args.observer_height_m,
            observer_height_m=location.location_height_m + args.observer_height_m,
        )
    return location, sorted(local.issues)


def location_payload(location: ResolvedLocation) -> dict[str, Any]:
    return asdict(location)
