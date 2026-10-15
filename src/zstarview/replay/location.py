"""Resolve observing inputs without downloading elevation/building data."""

from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

from ..location_resolver.resolve import (
    BUILDING_TOP_FETCH_RADIUS_KM,
    LocationDataPolicy,
    ResolvedLocation,
    _find_building_top_height_m,
    resolve_launch_location,
)
from ..terrain import GeoTiffDem, build_download_bbox, sample_ground_elevation
from ..terrain.dem import collect_copernicus_tile_keys
from .launch import load_latest_matching_place_resolution


class LocalLocationData:
    def __init__(
        self,
        root: Path,
        *,
        inspect_local_sources: bool = False,
        source_paths: set[Path] | None = None,
    ) -> None:
        self.root = root
        self.inspect_local_sources = inspect_local_sources
        self.source_paths = source_paths if source_paths is not None else set()
        self.issues: set[str] = set()

    def ground(self, lat: float, lon: float) -> float:
        if not self.inspect_local_sources:
            self.issues.add("ground-elevation-unverified")
            return 0.0
        bbox = build_download_bbox(lat_deg=lat, lon_deg=lon, radius_km=0.25)
        dem_root = self.root / "copernicus-dem"
        paths = tuple(
            dem_root / key
            for key in collect_copernicus_tile_keys(bbox)
            if (dem_root / key).is_file()
        )
        if paths:
            dem = None
            try:
                dem = GeoTiffDem(paths, default_elevation_m=None)
                grid = dem.build_grid(bbox)
                ground = sample_ground_elevation(
                    grid,
                    latitude_deg=lat,
                    longitude_deg=lon,
                    dem_resampling="bilinear",
                )
                self.source_paths.update(paths)
                return ground
            except Exception:
                pass
            finally:
                if dem is not None:
                    dem.close()
        self.issues.add("ground-elevation-unverified")
        return 0.0

    def building(self, lat: float, lon: float) -> float | None:
        if not self.inspect_local_sources:
            self.issues.add("building-top-unverified")
            return None
        try:
            from ..data.derived_tile_cache import (
                parse_derived_tile_buildings,
                select_derived_tile_envelopes,
            )
            from ..data.building_source import select_prepared_building_source

            source = select_prepared_building_source(
                observer_lat_deg=lat,
                observer_lon_deg=lon,
                radius_km=BUILDING_TOP_FETCH_RADIUS_KM,
                plateau_root_dir=self.root / "plateau_buildings",
            )
            if source.source == "plateau":
                buildings: list[Any] = []
                for derived_dir in source.derived_dirs:
                    envelopes = select_derived_tile_envelopes(
                        derived_dir,
                        observer_lat_deg=lat,
                        observer_lon_deg=lon,
                        radius_km=BUILDING_TOP_FETCH_RADIUS_KM,
                    )
                    for envelope in envelopes:
                        self.source_paths.add(envelope.path)
                        buildings.extend(parse_derived_tile_buildings(envelope.path))
                if buildings:
                    height = _find_building_top_height_m(
                        tuple(buildings), lon_deg=lon, lat_deg=lat
                    )
                    if height is not None:
                        return height
        except (OSError, ValueError, RuntimeError):
            pass
        self.issues.add("building-top-unverified")
        return None


def resolve_location(
    args: Any,
    root: Path,
    *,
    inspect_local_sources: bool = False,
    source_paths: set[Path] | None = None,
) -> tuple[ResolvedLocation, list[str]]:
    local = LocalLocationData(
        root,
        inspect_local_sources=inspect_local_sources,
        source_paths=source_paths,
    )
    city, place = args.city, args.place
    if city == "auto":
        raise ValueError(
            "Automatic IP location is unavailable in replay; specify a location"
        )
    previous_place_result = (
        load_latest_matching_place_resolution(
            query=place,
            countrycode=args.place_countrycode,
            language=args.place_lang,
            cache_root=root,
        )
        if isinstance(place, str) and place.strip()
        else None
    )
    location = resolve_launch_location(
        city,
        place_query=place,
        place_countrycode=args.place_countrycode,
        place_lang=args.place_lang,
        use_building_top=args.use_building_top,
        save_last_city_func=lambda _: None,
        previous_place_result=previous_place_result,
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
