#!/usr/bin/env python3
"""Find strong straight-line brightness boundaries in a night-light GeoTIFF.

The script scans a roughly 128 km square window over one or more EPSG:4326
night-light tiles. Each window is split by a line through its center. The line
direction is tested at 0, 45, 90, and 135 degrees, and the difference between
the two half-window brightness sums is used as the score.

The angle is measured counter-clockwise from east: 0 degrees is an
east-west boundary and 90 degrees is a north-south boundary. Coordinates are
reported as latitude and longitude in decimal degrees, matching zstarview's
location argument convention.
"""

from __future__ import annotations

import argparse
import heapq
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from affine import Affine
from rasterio.enums import Resampling


EARTH_KM_PER_DEGREE = 111.195
BOUNDARY_ANGLES_DEG = (0.0, 45.0, 90.0, 135.0)
DEFAULT_DARK_SIDE_OFFSET_KM = 60.0


@dataclass(frozen=True, slots=True)
class Candidate:
    score: float
    longitude: float
    latitude: float
    angle_deg: float
    bright_sum: float
    dark_sum: float
    ratio: float
    bright_side_azimuth_deg: float
    dark_side_azimuth_deg: float
    recommended_latitude: float
    recommended_longitude: float
    nonzero_pixels: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find strong straight-line boundaries in night-light GeoTIFFs."
    )
    parser.add_argument(
        "source",
        type=Path,
        nargs="+",
        help="One or more EPSG:4326 single-band night-light GeoTIFFs.",
    )
    parser.add_argument(
        "--window-km",
        type=float,
        default=128.0,
        help="Side length of each square scan window in km (default: 128).",
    )
    parser.add_argument(
        "--step-km",
        type=float,
        default=16.0,
        help="Distance between scan-window centers in km (default: 16).",
    )
    parser.add_argument(
        "--sample-km",
        type=float,
        default=4.0,
        help="Approximate raster cell size used during scanning (default: 4).",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Number of highest-scoring candidates to print (default: 10).",
    )
    parser.add_argument(
        "--bound-lat",
        metavar="MIN,MAX",
        help="Latitude scan range, for example 24,50.",
    )
    parser.add_argument(
        "--bound-lon",
        metavar="MIN,MAX",
        help="Longitude scan range, for example -125,-66.",
    )
    parser.add_argument(
        "--min-separation-km",
        type=float,
        default=128.0,
        help=(
            "Minimum distance between reported candidate centers in km "
            "(default: 128)."
        ),
    )
    parser.add_argument(
        "--dark-side-offset-km",
        type=float,
        default=DEFAULT_DARK_SIDE_OFFSET_KM,
        help=(
            "Distance to move from the candidate center toward the dark side "
            f"(default: {DEFAULT_DARK_SIDE_OFFSET_KM:g})."
        ),
    )
    # argparse treats a comma-separated value beginning with '-' as another
    # option. Normalize the friendly `--bound-lon "-125,-66"` spelling to an
    # equals form before handing the arguments to argparse.
    argv = list(sys.argv[1:])
    for option in ("--bound-lat", "--bound-lon"):
        try:
            option_index = argv.index(option)
        except ValueError:
            continue
        if option_index + 1 < len(argv):
            argv[option_index] = f"{option}={argv[option_index + 1]}"
            del argv[option_index + 1]
    return parser.parse_args(argv)


def parse_bound(value: str, name: str, minimum: float, maximum: float) -> tuple[float, float]:
    parts = value.split(",")
    if len(parts) != 2:
        raise ValueError(f"{name} must contain two comma-separated numbers")
    try:
        lower, upper = (float(part.strip()) for part in parts)
    except ValueError as exc:
        raise ValueError(f"{name} must contain two comma-separated numbers") from exc
    if not minimum <= lower < upper <= maximum:
        raise ValueError(
            f"{name} must satisfy {minimum} <= min < max <= {maximum}"
        )
    return lower, upper


def load_raster(path: Path, sample_km: float) -> tuple[np.ndarray, rasterio.Affine]:
    with rasterio.open(path) as dataset:
        if dataset.crs is None or dataset.crs.to_epsg() != 4326:
            raise ValueError(f"{path}: source CRS must be EPSG:4326, got {dataset.crs!s}")
        if dataset.count < 1:
            raise ValueError(f"{path}: source has no raster bands")
        if sample_km <= 0:
            raise ValueError("sample-km must be positive")

        native_km = abs(dataset.res[0]) * EARTH_KM_PER_DEGREE
        reduction = max(1, int(round(sample_km / native_km)))
        output_width = max(1, math.ceil(dataset.width / reduction))
        output_height = max(1, math.ceil(dataset.height / reduction))
        values = dataset.read(
            1,
            out_shape=(output_height, output_width),
            resampling=Resampling.average,
            masked=True,
        ).filled(0.0)
        transform = dataset.transform * dataset.transform.scale(
            dataset.width / output_width,
            dataset.height / output_height,
        )

    return np.asarray(values, dtype=np.float64), transform


def crop_to_bounds(
    values: np.ndarray,
    transform: rasterio.Affine,
    *,
    latitude_bounds: tuple[float, float] | None,
    longitude_bounds: tuple[float, float] | None,
) -> tuple[np.ndarray, rasterio.Affine]:
    if latitude_bounds is None and longitude_bounds is None:
        return values, transform
    if latitude_bounds is None or longitude_bounds is None:
        raise ValueError("bound-lat and bound-lon must be specified together")
    if transform.a <= 0 or transform.e >= 0:
        raise ValueError("bound cropping requires a north-up EPSG:4326 raster")

    west, east = longitude_bounds
    south, north = latitude_bounds
    left = transform.c
    top = transform.f
    pixel_width = transform.a
    pixel_height = abs(transform.e)
    col_start = max(0, math.floor((west - left) / pixel_width))
    col_stop = min(values.shape[1], math.ceil((east - left) / pixel_width))
    row_start = max(0, math.floor((top - north) / pixel_height))
    row_stop = min(values.shape[0], math.ceil((top - south) / pixel_height))
    if col_start >= col_stop or row_start >= row_stop:
        return values[:0, :0], transform

    cropped_transform = transform * Affine.translation(col_start, row_start)
    return values[row_start:row_stop, col_start:col_stop], cropped_transform


def coordinate(transform: rasterio.Affine, row: int, col: int) -> tuple[float, float]:
    longitude, latitude = transform * (col + 0.5, row + 0.5)
    return float(longitude), float(latitude)


def offset_location(
    latitude: float, longitude: float, bearing_deg: float, distance_km: float
) -> tuple[float, float]:
    """Move a location by a short distance using a local spherical approximation."""
    bearing_rad = math.radians(bearing_deg)
    latitude_scale = max(0.05, math.cos(math.radians(latitude)))
    offset_latitude = (
        latitude + distance_km * math.cos(bearing_rad) / EARTH_KM_PER_DEGREE
    )
    offset_longitude = (
        longitude
        + distance_km
        * math.sin(bearing_rad)
        / (EARTH_KM_PER_DEGREE * latitude_scale)
    )
    return offset_latitude, offset_longitude


def retain_candidate(
    candidates: list[tuple[float, int, Candidate]], candidate: Candidate, limit: int
) -> None:
    item = (candidate.score, id(candidate), candidate)
    if len(candidates) < limit:
        heapq.heappush(candidates, item)
    elif candidate.score > candidates[0][0]:
        heapq.heapreplace(candidates, item)


def distance_km(first: Candidate, second: Candidate) -> float:
    """Return an approximate great-circle distance between two candidates."""
    latitude_rad = math.radians((first.latitude + second.latitude) * 0.5)
    longitude_delta_rad = math.radians(second.longitude - first.longitude)
    latitude_delta_rad = math.radians(second.latitude - first.latitude)
    east_km = (
        math.degrees(longitude_delta_rad)
        * math.cos(latitude_rad)
        * EARTH_KM_PER_DEGREE
    )
    north_km = math.degrees(latitude_delta_rad) * EARTH_KM_PER_DEGREE
    return math.hypot(east_km, north_km)


def select_diverse_candidates(
    candidates: list[Candidate], *, top: int, min_separation_km: float
) -> list[Candidate]:
    if min_separation_km < 0:
        raise ValueError("min-separation-km must not be negative")
    selected: list[Candidate] = []
    for candidate in candidates:
        if all(
            distance_km(candidate, previous) >= min_separation_km
            for previous in selected
        ):
            selected.append(candidate)
            if len(selected) >= top:
                break
    return selected


def scan_raster(
    values: np.ndarray,
    transform: rasterio.Affine,
    *,
    window_km: float,
    step_km: float,
    top: int,
    min_separation_km: float,
    dark_side_offset_km: float,
) -> list[Candidate]:
    if window_km <= 0 or step_km <= 0:
        raise ValueError("window-km and step-km must be positive")
    if top <= 0:
        raise ValueError("top must be positive")
    if dark_side_offset_km < 0:
        raise ValueError("dark-side-offset-km must not be negative")

    height, width = values.shape
    cell_lat_km = abs(transform.e) * EARTH_KM_PER_DEGREE
    cell_lon_deg = abs(transform.a)
    # Keep the best candidate in each geographic bucket before ranking. This
    # prevents one bright city from filling the candidate pool with overlapping
    # windows before spatial suppression is applied.
    bucket_best: dict[tuple[int, int], Candidate] = {}
    candidates: list[tuple[float, int, Candidate]] = []

    # The window width in raster columns changes with latitude because the
    # GeoTIFF is in EPSG:4326. The local east/north approximation is accurate
    # enough for a 128 km exploratory window.
    row_step = max(1, round(step_km / cell_lat_km))
    for center_row in range(0, height, row_step):
        _, latitude = coordinate(transform, center_row, 0)
        longitude_scale = max(0.05, math.cos(math.radians(latitude)))
        cell_lon_km = cell_lon_deg * EARTH_KM_PER_DEGREE * longitude_scale
        half_height = max(1, math.ceil(window_km / (2.0 * cell_lat_km)))
        half_width = max(1, math.ceil(window_km / (2.0 * cell_lon_km)))

        if center_row - half_height < 0 or center_row + half_height >= height:
            continue

        col_step = max(1, round(step_km / cell_lon_km))
        for center_col in range(0, width, col_step):
            if center_col - half_width < 0 or center_col + half_width >= width:
                continue

            row_start = center_row - half_height
            row_stop = center_row + half_height + 1
            col_start = center_col - half_width
            col_stop = center_col + half_width + 1
            window = values[row_start:row_stop, col_start:col_stop]

            rows, cols = np.indices(window.shape, dtype=np.float64)
            x_km = (cols - half_width) * cell_lon_km
            y_km = (half_height - rows) * cell_lat_km
            nonzero_pixels = int(np.count_nonzero(window > 0.0))
            if nonzero_pixels == 0:
                continue

            longitude, latitude = coordinate(transform, center_row, center_col)
            for angle_deg in BOUNDARY_ANGLES_DEG:
                angle_rad = math.radians(angle_deg)
                # The sign of the cross product with the boundary direction
                # assigns pixels to the two sides of the center line.
                signed_side = math.cos(angle_rad) * y_km - math.sin(angle_rad) * x_km
                first_sum = float(window[signed_side >= 0.0].sum())
                second_sum = float(window[signed_side < 0.0].sum())
                positive_side_azimuth = (-angle_deg) % 360.0
                if first_sum >= second_sum:
                    bright_side_azimuth = positive_side_azimuth
                else:
                    bright_side_azimuth = (positive_side_azimuth + 180.0) % 360.0
                dark_side_azimuth = (bright_side_azimuth + 180.0) % 360.0
                recommended_latitude, recommended_longitude = offset_location(
                    latitude,
                    longitude,
                    dark_side_azimuth,
                    dark_side_offset_km,
                )
                bright_sum = max(first_sum, second_sum)
                dark_sum = min(first_sum, second_sum)
                candidate = Candidate(
                    score=bright_sum - dark_sum,
                    longitude=longitude,
                    latitude=latitude,
                    angle_deg=angle_deg,
                    bright_sum=bright_sum,
                    dark_sum=dark_sum,
                    ratio=(bright_sum / dark_sum if dark_sum > 0.0 else math.inf),
                    bright_side_azimuth_deg=bright_side_azimuth,
                    dark_side_azimuth_deg=dark_side_azimuth,
                    recommended_latitude=recommended_latitude,
                    recommended_longitude=recommended_longitude,
                    nonzero_pixels=nonzero_pixels,
                )
                if min_separation_km > 0.0:
                    latitude_bin = math.floor(
                        candidate.latitude
                        / (min_separation_km / EARTH_KM_PER_DEGREE)
                    )
                    longitude_scale = max(
                        0.05, math.cos(math.radians(candidate.latitude))
                    )
                    longitude_bin = math.floor(
                        candidate.longitude
                        / (min_separation_km / (EARTH_KM_PER_DEGREE * longitude_scale))
                    )
                    bucket = (latitude_bin, longitude_bin)
                    previous = bucket_best.get(bucket)
                    if previous is None or candidate.score > previous.score:
                        bucket_best[bucket] = candidate
                else:
                    retain_candidate(candidates, candidate, max(top * 1000, top))

    if min_separation_km > 0.0:
        ranked = sorted(bucket_best.values(), key=lambda item: item.score, reverse=True)
    else:
        ranked = [
            item[2]
            for item in sorted(candidates, key=lambda item: item[0], reverse=True)
        ]
    return select_diverse_candidates(
        ranked,
        top=top,
        min_separation_km=min_separation_km,
    )


def main() -> int:
    args = parse_args()
    if (args.bound_lat is None) != (args.bound_lon is None):
        raise SystemExit("--bound-lat and --bound-lon must be specified together")
    latitude_bounds = (
        parse_bound(args.bound_lat, "bound-lat", -90.0, 90.0)
        if args.bound_lat is not None
        else None
    )
    longitude_bounds = (
        parse_bound(args.bound_lon, "bound-lon", -180.0, 180.0)
        if args.bound_lon is not None
        else None
    )
    for path in args.source:
        values, transform = load_raster(path, args.sample_km)
        values, transform = crop_to_bounds(
            values,
            transform,
            latitude_bounds=latitude_bounds,
            longitude_bounds=longitude_bounds,
        )
        if values.size == 0:
            print(f"source={path}")
            print("No raster pixels overlap the requested bounds.")
            continue
        candidates = scan_raster(
            values,
            transform,
            window_km=args.window_km,
            step_km=args.step_km,
            top=args.top,
            min_separation_km=args.min_separation_km,
            dark_side_offset_km=args.dark_side_offset_km,
        )
        print(f"source={path}")
        print(
            "rank latitude longitude recommended_latitude recommended_longitude "
            "angle_deg bright_side_azimuth_deg dark_side_azimuth_deg bright_sum "
            "dark_sum difference ratio nonzero_pixels"
        )
        for rank, candidate in enumerate(candidates, start=1):
            print(
                f"{rank} {candidate.latitude:.6f} {candidate.longitude:.6f} "
                f"{candidate.recommended_latitude:.6f} "
                f"{candidate.recommended_longitude:.6f} "
                f"{candidate.angle_deg:.0f} "
                f"{candidate.bright_side_azimuth_deg:.0f} "
                f"{candidate.dark_side_azimuth_deg:.0f} "
                f"{candidate.bright_sum:.6g} {candidate.dark_sum:.6g} "
                f"{candidate.score:.6g} {candidate.ratio:.6g} "
                f"{candidate.nonzero_pixels}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
