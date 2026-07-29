from __future__ import annotations

import datetime as dt
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Iterable
from typing import Any

from ..user_agent import build_user_agent
from .models import (
    TropicalCyclonePoint,
    TropicalCyclonePolygon,
    TropicalCycloneSnapshot,
    TropicalCycloneSnapshotCollection,
)

DEFAULT_SERVICE_URL = (
    "https://services9.arcgis.com/RHVPKKiFTONKtxq3/arcgis/rest/services/"
    "Active_Hurricanes_v1/FeatureServer"
)
DEFAULT_USER_AGENT = build_user_agent("tropical-cyclone")
DEFAULT_TIMEOUT_S = 30.0
CURRENT_OBSERVED_LAYER_ID = 1
FORECAST_LAYER_ID = 0
WIND_LAYER_IDS = (7, 8, 9, 11)
DEFAULT_FEATURE_LIMIT = 250


class TropicalCycloneFetchError(RuntimeError):
    pass


def fetch_json(url: str, *, timeout_s: float, user_agent: str) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"User-Agent": user_agent})
    try:
        with urllib.request.urlopen(request, timeout=timeout_s) as response:
            payload = response.read()
    except urllib.error.HTTPError as exc:
        raise TropicalCycloneFetchError(
            f"HTTP error for {url}: {exc.code} {exc.reason}"
        ) from exc
    except urllib.error.URLError as exc:
        raise TropicalCycloneFetchError(
            f"Network error for {url}: {exc.reason}"
        ) from exc
    try:
        loaded = json.loads(payload.decode("utf-8"))
    except Exception as exc:
        raise TropicalCycloneFetchError(f"Invalid JSON from {url}") from exc
    if not isinstance(loaded, dict):
        raise TropicalCycloneFetchError(f"Expected JSON object from {url}")
    return loaded


def query_json(
    base_url: str,
    path: str,
    *,
    params: dict[str, str | int | float | bool],
    timeout_s: float,
    user_agent: str,
) -> dict[str, Any]:
    query = urllib.parse.urlencode(params)
    url = f"{base_url.rstrip('/')}/{path.lstrip('/')}?{query}"
    return fetch_json(url, timeout_s=timeout_s, user_agent=user_agent)


def service_label(service_url: str, service: dict[str, Any]) -> str:
    for key in ("name", "serviceDescription", "description", "documentInfo"):
        value = service.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
        if isinstance(value, dict):
            title = value.get("Title") or value.get("title")
            if isinstance(title, str) and title.strip():
                return title.strip()
    return service_url.rstrip("/").rsplit("/", 1)[-1]


def iter_layer_ids(service: dict[str, Any]) -> Iterable[int]:
    layers = service.get("layers")
    if not isinstance(layers, list):
        return ()
    result: list[int] = []
    for layer in layers:
        if not isinstance(layer, dict):
            continue
        layer_id = layer.get("id")
        if isinstance(layer_id, int):
            result.append(layer_id)
    return result


def _feature_attrs(feature: dict[str, Any]) -> dict[str, Any]:
    attrs = feature.get("attributes")
    if isinstance(attrs, dict):
        return attrs
    return {}


def _geometry(feature: dict[str, Any]) -> dict[str, Any] | None:
    geometry = feature.get("geometry")
    if isinstance(geometry, dict):
        return geometry
    return None


def _parse_epoch_ms(value: object) -> dt.datetime | None:
    if not isinstance(value, (int, float)):
        return None
    stamp = dt.datetime.fromtimestamp(float(value) / 1000.0, tz=dt.timezone.utc)
    return stamp


def _parse_label_time_utc(label: object) -> dt.datetime | None:
    if not isinstance(label, str) or not label.strip():
        return None
    text = label.strip()
    match = re.fullmatch(
        r"(?P<date>\d{4}-\d{2}-\d{2}) "
        r"(?P<hour>\d{1,2}):(?P<minute>\d{2})"
        r"(?: (?P<ampm>AM|PM))?"
        r"(?: [A-Za-z]{3})? UTC",
        text,
    )
    if match is not None:
        hour = int(match.group("hour"))
        minute = int(match.group("minute"))
        ampm = match.group("ampm")
        if ampm is not None:
            hour %= 12
            if ampm == "PM":
                hour += 12
        year_s, month_s, day_s = match.group("date").split("-")
        return dt.datetime(
            int(year_s),
            int(month_s),
            int(day_s),
            hour,
            minute,
            tzinfo=dt.timezone.utc,
        )
    return None


def _parse_valid_time_utc(attrs: dict[str, Any], *, prefer_label: bool = False) -> dt.datetime | None:
    if prefer_label:
        parsed_label = _parse_label_time_utc(attrs.get("FLDATELBL"))
        if parsed_label is not None:
            return parsed_label

    valid_time = attrs.get("VALIDTIME")
    parsed = _parse_epoch_ms(valid_time)
    if parsed is not None:
        return parsed

    if not prefer_label:
        parsed_label = _parse_label_time_utc(attrs.get("FLDATELBL"))
        if parsed_label is not None:
            return parsed_label

    return _parse_epoch_ms(attrs.get("DTG"))


def _parse_point(feature: dict[str, Any], *, label_field: str | None = None) -> TropicalCyclonePoint | None:
    attrs = _feature_attrs(feature)
    geometry = _geometry(feature)
    if geometry is None:
        return None
    x = geometry.get("x")
    y = geometry.get("y")
    if not isinstance(x, (int, float)) or not isinstance(y, (int, float)):
        return None
    tau = attrs.get("TAU")
    if not isinstance(tau, (int, float)):
        tau = attrs.get("FCSTPRD")
    label = None
    if label_field is not None:
        value = attrs.get(label_field)
        if isinstance(value, str) and value.strip():
            label = value.strip()
    valid_time_utc = _parse_valid_time_utc(attrs, prefer_label=label_field == "FLDATELBL")
    if valid_time_utc is None:
        advdate = _parse_epoch_ms(attrs.get("ADVDATE"))
        if advdate is not None and isinstance(tau, (int, float)):
            valid_time_utc = advdate + dt.timedelta(hours=float(tau))
        elif advdate is not None:
            valid_time_utc = advdate
    return TropicalCyclonePoint(
        lat_deg=float(y),
        lon_deg=float(x),
        valid_time_utc=valid_time_utc,
        label=label,
        tau_hr=int(tau) if isinstance(tau, (int, float)) else None,
        maxwind_kt=float(attrs["MAXWIND"]) if isinstance(attrs.get("MAXWIND"), (int, float)) else None,
        gust_kt=float(attrs["GUST"]) if isinstance(attrs.get("GUST"), (int, float)) else None,
        mslp_hpa=float(attrs["MSLP"]) if isinstance(attrs.get("MSLP"), (int, float)) else None,
        dvlbl=attrs.get("DVLBL") if isinstance(attrs.get("DVLBL"), str) else None,
    )


def _polygon_from_feature(feature: dict[str, Any], *, layer_id: int, name: str) -> TropicalCyclonePolygon | None:
    geometry = _geometry(feature)
    if geometry is None:
        return None
    rings_data = geometry.get("rings")
    if not isinstance(rings_data, list):
        return None
    rings: list[tuple[tuple[float, float], ...]] = []
    for ring in rings_data:
        if not isinstance(ring, list):
            continue
        coords: list[tuple[float, float]] = []
        for point in ring:
            if (
                isinstance(point, list)
                and len(point) >= 2
                and isinstance(point[0], (int, float))
                and isinstance(point[1], (int, float))
            ):
                coords.append((float(point[1]), float(point[0])))
        if coords:
            rings.append(tuple(coords))
    if not rings:
        return None
    return TropicalCyclonePolygon(layer_id=layer_id, name=name, rings=tuple(rings))


def _storm_identity(attrs: dict[str, Any]) -> str | None:
    storm_id = attrs.get("ATCFID")
    if isinstance(storm_id, str) and storm_id.strip():
        return storm_id.strip()
    storm_num = attrs.get("STORMNUM")
    if storm_num is not None:
        storm_num_text = str(storm_num).strip()
        if storm_num_text:
            basin = attrs.get("BASIN")
            basin_text = basin.strip() if isinstance(basin, str) and basin.strip() else ""
            if basin_text:
                return f"{basin_text}:{storm_num_text}"
            return storm_num_text
    storm_name = attrs.get("STORMNAME")
    if isinstance(storm_name, str) and storm_name.strip():
        basin = attrs.get("BASIN")
        basin_text = basin.strip() if isinstance(basin, str) and basin.strip() else ""
        if basin_text:
            return f"{basin_text}:{storm_name.strip()}"
        return storm_name.strip()
    return None


def _storm_name_for_attrs(attrs: dict[str, Any]) -> str | None:
    storm_name = attrs.get("STORMNAME")
    if isinstance(storm_name, str) and storm_name.strip():
        return storm_name.strip()
    return None


def _storm_basin_for_attrs(attrs: dict[str, Any]) -> str | None:
    basin = attrs.get("BASIN")
    if isinstance(basin, str) and basin.strip():
        return basin.strip()
    return None


def _query_layer(
    service_url: str,
    layer_id: int,
    *,
    where: str,
    limit: int,
    timeout_s: float,
    user_agent: str,
    order_by_fields: str | None = None,
) -> dict[str, Any]:
    params: dict[str, str | int | float | bool] = {
        "where": where,
        "outFields": "*",
        "returnGeometry": "true",
        "resultRecordCount": max(1, int(limit)),
        "f": "json",
    }
    if order_by_fields is not None:
        params["orderByFields"] = order_by_fields
    return query_json(
        service_url,
        f"{layer_id}/query",
        params=params,
        timeout_s=timeout_s,
        user_agent=user_agent,
    )


def _matches_storm(feature: dict[str, Any], storm_name: str, basin: str | None) -> bool:
    attrs = _feature_attrs(feature)
    name = attrs.get("STORMNAME")
    if not isinstance(name, str) or name != storm_name:
        return False
    if basin is None:
        return True
    value = attrs.get("BASIN")
    return isinstance(value, str) and value == basin


def _latest_forecast_advdate(features: list[dict[str, Any]], storm_name: str, basin: str | None) -> int | None:
    for feature in features:
        if not isinstance(feature, dict) or not _matches_storm(feature, storm_name, basin):
            continue
        advdate = _feature_attrs(feature).get("ADVDATE")
        if isinstance(advdate, (int, float)):
            return int(advdate)
    return None


def _sort_forecast_features(features: list[dict[str, Any]]) -> list[dict[str, Any]]:
    def sort_key(feature: dict[str, Any]) -> tuple[int, int, int]:
        attrs = _feature_attrs(feature)
        valid_time = attrs.get("VALIDTIME")
        tau = attrs.get("TAU") or attrs.get("FCSTPRD") or 0
        object_id = attrs.get("OBJECTID") or 0
        return (
            int(valid_time) if isinstance(valid_time, (int, float)) else 0,
            int(tau) if isinstance(tau, (int, float)) else 0,
            int(object_id) if isinstance(object_id, (int, float)) else 0,
        )

    return sorted(features, key=sort_key)


def fetch_latest_observed_feature(
    *,
    service_url: str = DEFAULT_SERVICE_URL,
    timeout_s: float = DEFAULT_TIMEOUT_S,
    user_agent: str = DEFAULT_USER_AGENT,
) -> dict[str, Any] | None:
    query = _query_layer(
        service_url,
        CURRENT_OBSERVED_LAYER_ID,
        where="1=1",
        limit=1,
        timeout_s=timeout_s,
        user_agent=user_agent,
        order_by_fields="DTG DESC, OBJECTID DESC",
    )
    error = query.get("error")
    if isinstance(error, dict):
        raise TropicalCycloneFetchError(str(error.get("message", "observed query failed")))
    features = query.get("features")
    if not isinstance(features, list):
        raise TropicalCycloneFetchError("Observed query returned no features array")
    if not features:
        return None
    first = features[0]
    if not isinstance(first, dict):
        raise TropicalCycloneFetchError("Unexpected observed position payload")
    return first


def fetch_active_hurricanes_snapshot(
    *,
    service_url: str = DEFAULT_SERVICE_URL,
    timeout_s: float = DEFAULT_TIMEOUT_S,
    user_agent: str = DEFAULT_USER_AGENT,
) -> TropicalCycloneSnapshotCollection:
    service = fetch_json(
        f"{service_url}?f=json",
        timeout_s=timeout_s,
        user_agent=user_agent,
    )
    service_name = service_label(service_url, service)
    observed_query = _query_layer(
        service_url,
        CURRENT_OBSERVED_LAYER_ID,
        where="1=1",
        limit=DEFAULT_FEATURE_LIMIT,
        timeout_s=timeout_s,
        user_agent=user_agent,
        order_by_fields="DTG DESC, OBJECTID DESC",
    )
    observed_error = observed_query.get("error")
    if isinstance(observed_error, dict):
        raise TropicalCycloneFetchError(str(observed_error.get("message", "observed query failed")))
    observed_features = observed_query.get("features")
    if not isinstance(observed_features, list):
        raise TropicalCycloneFetchError("Observed query returned no features array")
    observed_by_storm_key: dict[str, dict[str, Any]] = {}
    observed_order: list[str] = []
    for feature in observed_features:
        if not isinstance(feature, dict):
            continue
        attrs = _feature_attrs(feature)
        storm_key = _storm_identity(attrs)
        if storm_key is None or storm_key in observed_by_storm_key:
            continue
        observed_by_storm_key[storm_key] = feature
        observed_order.append(storm_key)

    forecast_limit = DEFAULT_FEATURE_LIMIT
    forecast_query = _query_layer(
        service_url,
        FORECAST_LAYER_ID,
        where="1=1",
        limit=forecast_limit,
        timeout_s=timeout_s,
        user_agent=user_agent,
        order_by_fields="ADVDATE DESC, FCSTPRD ASC, OBJECTID ASC",
    )
    forecast_error = forecast_query.get("error")
    if isinstance(forecast_error, dict):
        raise TropicalCycloneFetchError(str(forecast_error.get("message", "forecast query failed")))
    forecast_features = forecast_query.get("features")
    if not isinstance(forecast_features, list):
        raise TropicalCycloneFetchError("Forecast query returned no features array")
    active_storm_keys = set(observed_by_storm_key.keys())
    forecast_features_by_storm_key: dict[str, list[dict[str, Any]]] = {key: [] for key in active_storm_keys}
    for feature in forecast_features:
        if not isinstance(feature, dict):
            continue
        attrs = _feature_attrs(feature)
        storm_key = _storm_identity(attrs)
        if storm_key is None or storm_key not in active_storm_keys:
            continue
        forecast_features_by_storm_key.setdefault(storm_key, []).append(feature)

    wind_features_by_storm_key: dict[str, list[tuple[int, str, dict[str, Any]]]] = {
        key: [] for key in active_storm_keys
    }
    for layer_id in WIND_LAYER_IDS:
        wind_query = _query_layer(
            service_url,
            layer_id,
            where="1=1",
            limit=forecast_limit,
            timeout_s=timeout_s,
            user_agent=user_agent,
        )
        wind_error = wind_query.get("error")
        if isinstance(wind_error, dict):
            continue
        wind_features = wind_query.get("features")
        if not isinstance(wind_features, list):
            continue
        layer_meta = service.get("layers")
        layer_name = f"layer-{layer_id}"
        if isinstance(layer_meta, list):
            for layer in layer_meta:
                if isinstance(layer, dict) and layer.get("id") == layer_id:
                    candidate_name = layer.get("name")
                    if isinstance(candidate_name, str) and candidate_name.strip():
                        layer_name = candidate_name.strip()
                    break
        for feature in wind_features:
            if not isinstance(feature, dict):
                continue
            attrs = _feature_attrs(feature)
            storm_key = _storm_identity(attrs)
            if storm_key is None or storm_key not in active_storm_keys:
                continue
            wind_features_by_storm_key.setdefault(storm_key, []).append((layer_id, layer_name, feature))

    snapshots: list[TropicalCycloneSnapshot] = []
    refreshed_at_utc = dt.datetime.now(dt.timezone.utc)
    for storm_key in observed_order:
        observed_feature = observed_by_storm_key.get(storm_key)
        if observed_feature is None:
            continue
        observed_attrs = _feature_attrs(observed_feature)
        storm_name = _storm_name_for_attrs(observed_attrs)
        if storm_name is None:
            continue
        basin_text = _storm_basin_for_attrs(observed_attrs)
        storm_id = observed_attrs.get("ATCFID")
        if not isinstance(storm_id, str) or not storm_id.strip():
            storm_id = observed_attrs.get("STORMNUM")
            storm_id = str(storm_id).strip() if storm_id is not None else None
        advdate = observed_attrs.get("ADVDATE")
        advdate_utc = _parse_epoch_ms(advdate)
        observed_point = _parse_point(observed_feature, label_field="DATELBL")
        if observed_point is None:
            continue

        storm_forecast_features = forecast_features_by_storm_key.get(storm_key, [])
        forecast_advdate = _latest_forecast_advdate(
            storm_forecast_features,
            storm_name,
            basin_text,
        )
        filtered_forecast = [
            feature
            for feature in storm_forecast_features
            if _matches_storm(feature, storm_name, basin_text)
            and (
                forecast_advdate is None
                or _feature_attrs(feature).get("ADVDATE") == forecast_advdate
            )
        ]
        forecast_points: list[TropicalCyclonePoint] = []
        for feature in _sort_forecast_features(filtered_forecast):
            point = _parse_point(feature, label_field="FLDATELBL")
            if point is not None:
                forecast_points.append(point)

        wind_polygons: list[TropicalCyclonePolygon] = []
        for layer_id, layer_name, feature in wind_features_by_storm_key.get(storm_key, []):
            polygon = _polygon_from_feature(feature, layer_id=layer_id, name=layer_name)
            if polygon is not None:
                wind_polygons.append(polygon)

        snapshots.append(
            TropicalCycloneSnapshot(
                storm_name=storm_name,
                basin=basin_text,
                advdate_utc=advdate_utc,
                observed_position=observed_point,
                forecast_positions=tuple(forecast_points),
                wind_polygons=tuple(wind_polygons),
                source_url=service_url,
                service_name=service_name,
                refreshed_at_utc=refreshed_at_utc,
                current_storm_id=storm_id,
            )
        )

    snapshots.sort(
        key=lambda snapshot: (
            snapshot.advdate_utc or dt.datetime.min.replace(tzinfo=dt.timezone.utc),
            snapshot.storm_name,
            snapshot.current_storm_id or "",
        ),
        reverse=True,
    )
    return TropicalCycloneSnapshotCollection(
        snapshots=tuple(snapshots),
        source_url=service_url,
        service_name=service_name,
        refreshed_at_utc=refreshed_at_utc,
    )
