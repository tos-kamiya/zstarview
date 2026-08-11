from __future__ import annotations

import json
import math
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Any

import numpy as np
from pyproj import Geod

from .location_resolver.place_projection import project_place_targets_to_altaz
from .road_night_lights import build_road_night_light_ground_sampler
from .types import ViewerData
from .user_agent import build_user_agent

OPEN_METEO_ENDPOINT = "https://api.open-meteo.com/v1/forecast"
PRECIPITATION_REFRESH_SECONDS = 10 * 60
PRECIPITATION_SAMPLE_COUNT = 24
PRECIPITATION_MIN_DISTANCE_KM = 8.0
PRECIPITATION_MAX_DISTANCE_KM = 32.0
PRECIPITATION_GOLDEN_ANGLE_DEG = 137.507764
PRECIPITATION_MIN_RATE_MM_H = 0.1
PRECIPITATION_MAX_DISPLAY_HEIGHT_DEG = 16.0
PRECIPITATION_NEAR_WIDTH_PX = 4.0
PRECIPITATION_FAR_WIDTH_PX = 1.5
PRECIPITATION_COLUMN_COLOR_RGB = (70, 150, 255)
PRECIPITATION_CACHE_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class PrecipitationSampleLocation:
    latitude_deg: float
    longitude_deg: float
    azimuth_deg: float
    distance_km: float


@dataclass(frozen=True)
class PrecipitationForecastValue:
    amount_mm: float | None
    rain_mm: float | None
    showers_mm: float | None
    rate_mm_h: float | None
    forecast_time_utc: datetime
    interval_seconds: int


@dataclass(frozen=True)
class PrecipitationSnapshot:
    samples: tuple[PrecipitationSampleLocation, ...]
    values: tuple[PrecipitationForecastValue, ...]
    fetched_at_utc: datetime


@dataclass(frozen=True)
class ProjectedPrecipitationColumn:
    base_alt_deg: float
    base_az_deg: float
    top_alt_deg: float
    top_az_deg: float
    distance_km: float
    rate_mm_h: float


def generate_precipitation_samples(
    observer_latitude_deg: float,
    observer_longitude_deg: float,
) -> tuple[PrecipitationSampleLocation, ...]:
    geod = Geod(ellps="WGS84")
    result: list[PrecipitationSampleLocation] = []
    for index in range(PRECIPITATION_SAMPLE_COUNT):
        azimuth_deg = (index * PRECIPITATION_GOLDEN_ANGLE_DEG) % 360.0
        t = (index + 0.5) / PRECIPITATION_SAMPLE_COUNT
        distance_km = math.sqrt(
            PRECIPITATION_MIN_DISTANCE_KM**2
            + (
                PRECIPITATION_MAX_DISTANCE_KM**2
                - PRECIPITATION_MIN_DISTANCE_KM**2
            )
            * t
        )
        longitude_deg, latitude_deg, _ = geod.fwd(
            float(observer_longitude_deg),
            float(observer_latitude_deg),
            azimuth_deg,
            distance_km * 1000.0,
        )
        result.append(
            PrecipitationSampleLocation(
                latitude_deg=float(latitude_deg),
                longitude_deg=float(longitude_deg),
                azimuth_deg=float(azimuth_deg),
                distance_km=float(distance_km),
            )
        )
    return tuple(result)


def precipitation_rate_mm_h(
    amount_mm: float | None, interval_seconds: int
) -> float | None:
    if interval_seconds <= 0:
        raise ValueError("Open-Meteo precipitation interval must be positive")
    if amount_mm is None:
        return None
    amount = float(amount_mm)
    if not math.isfinite(amount) or amount < 0.0:
        raise ValueError("Open-Meteo precipitation amount is invalid")
    return amount * 3600.0 / float(interval_seconds)


def precipitation_column_display_height_deg(rate_mm_h: float) -> float:
    rate = max(0.0, float(rate_mm_h))
    return min(
        PRECIPITATION_MAX_DISPLAY_HEIGHT_DEG,
        3.0 * math.log2(1.0 + rate),
    )


def precipitation_column_width_px(distance_km: float) -> float:
    distance_span_km = (
        PRECIPITATION_MAX_DISTANCE_KM - PRECIPITATION_MIN_DISTANCE_KM
    )
    t = (float(distance_km) - PRECIPITATION_MIN_DISTANCE_KM) / distance_span_km
    t = min(1.0, max(0.0, t))
    return PRECIPITATION_NEAR_WIDTH_PX + (
        (PRECIPITATION_FAR_WIDTH_PX - PRECIPITATION_NEAR_WIDTH_PX) * t
    )


def precipitation_cache_key(
    samples: tuple[PrecipitationSampleLocation, ...],
) -> tuple[object, ...]:
    coordinates = tuple(
        (round(sample.latitude_deg, 5), round(sample.longitude_deg, 5))
        for sample in samples
    )
    return (
        PRECIPITATION_CACHE_SCHEMA_VERSION,
        "open-meteo-noncommercial",
        coordinates,
        ("precipitation", "rain", "showers"),
        "mm",
        "nearest",
    )


def precipitation_snapshot_is_fresh(
    snapshot: PrecipitationSnapshot,
    *,
    now_utc: datetime,
) -> bool:
    age = now_utc.astimezone(timezone.utc) - snapshot.fetched_at_utc.astimezone(
        timezone.utc
    )
    return timedelta(0) <= age < timedelta(seconds=PRECIPITATION_REFRESH_SECONDS)


def _optional_amount(value: object) -> float | None:
    if value is None:
        return None
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ValueError("Open-Meteo precipitation component is invalid")
    amount = float(value)
    if not math.isfinite(amount) or amount < 0.0:
        raise ValueError("Open-Meteo precipitation component is invalid")
    return amount


def parse_open_meteo_response(
    payload: object,
    samples: tuple[PrecipitationSampleLocation, ...],
    *,
    fetched_at_utc: datetime,
) -> PrecipitationSnapshot:
    if not isinstance(payload, list) or len(payload) != len(samples):
        raise ValueError("Open-Meteo returned an unexpected location count")
    values: list[PrecipitationForecastValue] = []
    for item in payload:
        if not isinstance(item, dict):
            raise ValueError("Open-Meteo location response is invalid")
        current = item.get("current")
        units = item.get("current_units")
        if not isinstance(current, dict) or not isinstance(units, dict):
            raise ValueError("Open-Meteo current precipitation is missing")
        if units.get("precipitation") != "mm":
            raise ValueError("Open-Meteo precipitation unit is not mm")
        try:
            interval_seconds = int(current["interval"])
            forecast_time = datetime.fromisoformat(
                str(current["time"]).replace("Z", "+00:00")
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("Open-Meteo current time or interval is invalid") from exc
        if forecast_time.tzinfo is None:
            forecast_time = forecast_time.replace(tzinfo=timezone.utc)
        forecast_time = forecast_time.astimezone(timezone.utc)
        amount = _optional_amount(current.get("precipitation"))
        values.append(
            PrecipitationForecastValue(
                amount_mm=amount,
                rain_mm=_optional_amount(current.get("rain")),
                showers_mm=_optional_amount(current.get("showers")),
                rate_mm_h=precipitation_rate_mm_h(amount, interval_seconds),
                forecast_time_utc=forecast_time,
                interval_seconds=interval_seconds,
            )
        )
    return PrecipitationSnapshot(samples, tuple(values), fetched_at_utc)


def fetch_open_meteo_precipitation(
    samples: tuple[PrecipitationSampleLocation, ...],
    *,
    timeout_seconds: float = 20.0,
    opener: Any = urllib.request.urlopen,
) -> PrecipitationSnapshot:
    query = urllib.parse.urlencode(
        {
            "latitude": ",".join(f"{sample.latitude_deg:.5f}" for sample in samples),
            "longitude": ",".join(f"{sample.longitude_deg:.5f}" for sample in samples),
            "current": "precipitation,rain,showers",
            "precipitation_unit": "mm",
            "timezone": "GMT",
            "cell_selection": "nearest",
        }
    )
    request = urllib.request.Request(
        f"{OPEN_METEO_ENDPOINT}?{query}",
        headers={"User-Agent": build_user_agent("open-meteo-precipitation")},
    )
    with opener(request, timeout=float(timeout_seconds)) as response:
        payload = json.loads(response.read().decode("utf-8"))
    return parse_open_meteo_response(
        payload,
        samples,
        fetched_at_utc=datetime.now(timezone.utc),
    )


def project_precipitation_columns(
    snapshot: PrecipitationSnapshot,
    viewer_data: ViewerData,
) -> tuple[ProjectedPrecipitationColumn, ...]:
    latitudes = [sample.latitude_deg for sample in snapshot.samples]
    longitudes = [sample.longitude_deg for sample in snapshot.samples]
    try:
        sampler = build_road_night_light_ground_sampler(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            radius_km=PRECIPITATION_MAX_DISTANCE_KM,
        )
        ground = np.asarray(sampler(latitudes, longitudes), dtype=np.float64)
        observer_ground = float(
            np.asarray(
                sampler([viewer_data.lat_deg], [viewer_data.lon_deg]),
                dtype=np.float64,
            )[0]
        )
    except (OSError, RuntimeError, TypeError, ValueError):
        ground = np.zeros(len(snapshot.samples), dtype=np.float64)
        observer_ground = float(getattr(viewer_data, "ground_elevation_m", 0.0))
    observer_height = observer_ground + float(viewer_data.observer_height_m)
    bases = project_place_targets_to_altaz(
        observer_latitude_deg=float(viewer_data.lat_deg),
        observer_longitude_deg=float(viewer_data.lon_deg),
        observer_height_m=observer_height,
        target_latitude_deg=latitudes,
        target_longitude_deg=longitudes,
        target_height_m=ground,
    )
    result: list[ProjectedPrecipitationColumn] = []
    for value, base in zip(snapshot.values, bases, strict=True):
        rate = value.rate_mm_h
        if rate is None or rate < PRECIPITATION_MIN_RATE_MM_H:
            continue
        display_height_deg = precipitation_column_display_height_deg(rate)
        result.append(
            ProjectedPrecipitationColumn(
                base_alt_deg=base.alt_deg,
                base_az_deg=base.az_deg,
                top_alt_deg=min(90.0, base.alt_deg + display_height_deg),
                top_az_deg=base.az_deg,
                distance_km=base.distance_km,
                rate_mm_h=rate,
            )
        )
    return tuple(result)
