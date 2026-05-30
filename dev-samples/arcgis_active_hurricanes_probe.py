#!/usr/bin/env python3
"""Fetch active tropical cyclone data from a public ArcGIS FeatureServer.

This probe demonstrates that some ArcGIS-hosted hurricane/cyclone layers can be
queried without an API key or token when the item is publicly accessible.
It resolves the active storm from the latest observed position and then prints
the current position, forecast position, and wind polygon layers for that same
storm.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Iterable
from typing import Any

DEFAULT_SERVICE_URL = (
    "https://services9.arcgis.com/RHVPKKiFTONKtxq3/arcgis/rest/services/"
    "Active_Hurricanes_v1/FeatureServer"
)
DEFAULT_USER_AGENT = "zstarview-arcgis-active-hurricanes-probe/0.1"
DEFAULT_TIMEOUT_S = 30.0
DEFAULT_LIMIT = 3
WIND_LAYER_IDS = (7, 8, 9, 11)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Fetch public ArcGIS tropical cyclone data without an API key or token. "
            "The default service is NHC/Esri Active_Hurricanes_v1."
        )
    )
    parser.add_argument(
        "--service-url",
        default=DEFAULT_SERVICE_URL,
        help=f"ArcGIS FeatureServer URL (default: {DEFAULT_SERVICE_URL}).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=DEFAULT_LIMIT,
        help=f"Maximum number of features to print per layer (default: {DEFAULT_LIMIT}).",
    )
    parser.add_argument(
        "--layers",
        type=int,
        nargs="*",
        help="Optional layer IDs to query. If omitted, all layers in the service are queried.",
    )
    parser.add_argument(
        "--timeout-s",
        type=float,
        default=DEFAULT_TIMEOUT_S,
        help=f"HTTP timeout in seconds (default: {DEFAULT_TIMEOUT_S:.0f}).",
    )
    parser.add_argument(
        "--user-agent",
        default=DEFAULT_USER_AGENT,
        help=f"HTTP User-Agent string (default: {DEFAULT_USER_AGENT}).",
    )
    return parser


def fetch_json(url: str, *, timeout_s: float, user_agent: str) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"User-Agent": user_agent})
    try:
        with urllib.request.urlopen(request, timeout=timeout_s) as response:
            payload = response.read()
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"HTTP error for {url}: {exc.code} {exc.reason}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Network error for {url}: {exc.reason}") from exc
    try:
        loaded = json.loads(payload.decode("utf-8"))
    except Exception as exc:
        raise RuntimeError(f"Invalid JSON from {url}") from exc
    if not isinstance(loaded, dict):
        raise RuntimeError(f"Expected JSON object from {url}")
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


def format_coord(value: object) -> str:
    if isinstance(value, (int, float)):
        return f"{value:.4f}"
    return "?"


def format_epoch_ms(value: object) -> str | None:
    if not isinstance(value, (int, float)):
        return None
    seconds = float(value) / 1000.0
    stamp = dt.datetime.fromtimestamp(seconds, tz=dt.timezone.utc)
    return stamp.strftime("%Y-%m-%d %H:%M:%S UTC")


def format_observed_date(attributes: dict[str, Any]) -> str | None:
    dtg = attributes.get("DTG")
    formatted = format_epoch_ms(dtg)
    if formatted is not None:
        return f"DTG={formatted}"
    year = attributes.get("YEAR")
    month = attributes.get("MONTH")
    day = attributes.get("DAY")
    hhmm = attributes.get("HHMM")
    if all(value is not None for value in (year, month, day, hhmm)):
        return f"DTG={year}-{month}-{day} {hhmm}"
    return None


def format_valid_time(attributes: dict[str, Any]) -> str | None:
    valid_time = attributes.get("VALIDTIME")
    formatted = format_epoch_ms(valid_time)
    if formatted is not None:
        return f"VALIDTIME={formatted}"
    return None


def summarize_point_geometry(geometry: dict[str, Any] | None) -> str:
    if not geometry:
        return "geometry=none"
    x = geometry.get("x")
    y = geometry.get("y")
    return f"point=({format_coord(x)}, {format_coord(y)})"


def summarize_paths_geometry(geometry: dict[str, Any] | None) -> str:
    if not geometry:
        return "geometry=none"
    paths = geometry.get("paths")
    if isinstance(paths, list):
        path_count = len(paths)
        point_count = sum(len(path) for path in paths if isinstance(path, list))
        return f"polyline=paths:{path_count} points:{point_count}"
    rings = geometry.get("rings")
    if isinstance(rings, list):
        ring_count = len(rings)
        point_count = sum(len(ring) for ring in rings if isinstance(ring, list))
        return f"polygon=rings:{ring_count} points:{point_count}"
    return "geometry=unknown"


def pick_summary_fields(attributes: dict[str, Any]) -> list[str]:
    preferred = [
        "STORMNAME",
        "STORMTYPE",
        "DTG",
        "ADVDATE",
        "ADVISNUM",
        "FCSTPRD",
        "TAU",
        "VALIDTIME",
        "DATELBL",
        "FLDATELBL",
        "LAT",
        "LON",
        "MAXWIND",
        "GUST",
        "MSLP",
        "TCDIR",
        "TCSPD",
        "DVLBL",
        "SSNUM",
        "STORMNUM",
        "STORMSRC",
        "RADII",
        "TCWW",
        "STATUS",
        "BASIN",
        "ATCFID",
        "NAME",
    ]
    parts: list[str] = []
    for key in preferred:
        value = attributes.get(key)
        if value is not None:
            if key == "ADVDATE":
                formatted = format_epoch_ms(value)
                if formatted is not None:
                    parts.append(f"{key}={formatted}")
                    continue
            parts.append(f"{key}={value}")
    observed_date = format_observed_date(attributes)
    if observed_date is not None:
        parts.insert(0, observed_date)
    valid_time = format_valid_time(attributes)
    if valid_time is not None:
        parts.insert(1 if observed_date is not None else 0, valid_time)
    return parts


def print_feature(feature: dict[str, Any], geometry_type: str | None, *, label: str | None = None) -> None:
    attributes = feature.get("attributes")
    if not isinstance(attributes, dict):
        attributes = {}
    geometry = feature.get("geometry")
    if geometry_type == "esriGeometryPoint":
        geom_summary = summarize_point_geometry(geometry if isinstance(geometry, dict) else None)
    else:
        geom_summary = summarize_paths_geometry(geometry if isinstance(geometry, dict) else None)
    fields = pick_summary_fields(attributes)
    field_text = ", ".join(fields) if fields else "no key attributes"
    prefix = f"{label}: " if label else ""
    print(f"    - {prefix}{field_text}; {geom_summary}")


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


def escape_sql_literal(value: str) -> str:
    return value.replace("'", "''")


def where_for_storm(storm_name: str, basin: str | None = None, *, advdate: int | None = None) -> str:
    parts = [f"STORMNAME = '{escape_sql_literal(storm_name)}'"]
    if basin:
        parts.append(f"BASIN = '{escape_sql_literal(basin)}'")
    if advdate is not None:
        parts.append(f"ADVDATE = {advdate}")
    return " AND ".join(parts)


def query_params_for_layer(layer_id: int, limit: int, where: str) -> dict[str, str | int | float | bool]:
    params: dict[str, str | int | float | bool] = {
        "where": where,
        "outFields": "*",
        "returnGeometry": "true",
        "resultRecordCount": limit,
        "f": "json",
    }
    if layer_id == 1:
        params["orderByFields"] = "DTG DESC, OBJECTID DESC"
        params["resultRecordCount"] = 1
    elif layer_id == 0:
        params["orderByFields"] = "ADVDATE DESC, FCSTPRD ASC, OBJECTID ASC"
    return params


def query_layer(
    args: argparse.Namespace,
    layer_id: int,
    *,
    where: str,
    limit: int,
) -> dict[str, Any]:
    return query_json(
        args.service_url,
        f"{layer_id}/query",
        params=query_params_for_layer(layer_id, limit, where),
        timeout_s=args.timeout_s,
        user_agent=args.user_agent,
    )


def feature_attrs(feature: dict[str, Any]) -> dict[str, Any]:
    attrs = feature.get("attributes")
    if isinstance(attrs, dict):
        return attrs
    return {}


def matches_storm(feature: dict[str, Any], storm_name: str, basin: str | None = None) -> bool:
    attrs = feature_attrs(feature)
    name = attrs.get("STORMNAME")
    if not isinstance(name, str) or name != storm_name:
        return False
    if basin is None:
        return True
    value = attrs.get("BASIN")
    return isinstance(value, str) and value == basin


def wind_layer_filter(
    features: list[dict[str, Any]],
    storm_name: str,
    basin: str | None,
) -> list[dict[str, Any]]:
    return [feature for feature in features if matches_storm(feature, storm_name, basin)]


def main() -> int:
    args = build_arg_parser().parse_args()

    service = fetch_json(
        f"{args.service_url}?f=json",
        timeout_s=args.timeout_s,
        user_agent=args.user_agent,
    )

    service_name = service_label(args.service_url, service)
    service_description = service.get("serviceDescription") or service.get("description") or ""
    print(f"Service: {service_name}")
    if isinstance(service_description, str) and service_description.strip():
        print(f"Description: {service_description.strip()}")

    available_layer_ids = list(iter_layer_ids(service))
    if not available_layer_ids:
        print("No layers found in service metadata.", file=sys.stderr)
        return 1

    if args.layers:
        layer_ids = [layer_id for layer_id in args.layers if layer_id in available_layer_ids]
    else:
        layer_ids = [layer_id for layer_id in (1, 0, *WIND_LAYER_IDS) if layer_id in available_layer_ids]

    layer_meta_cache: dict[int, dict[str, Any]] = {}

    def layer_meta(layer_id: int) -> dict[str, Any]:
        cached = layer_meta_cache.get(layer_id)
        if cached is not None:
            return cached
        meta = fetch_json(
            f"{args.service_url}/{layer_id}?f=json",
            timeout_s=args.timeout_s,
            user_agent=args.user_agent,
        )
        layer_meta_cache[layer_id] = meta
        return meta

    print("Layers:")
    current_attrs: dict[str, Any] = {}
    if 1 in layer_ids:
        meta = layer_meta(1)
        print(f"- [1] {meta.get('name', 'Observed Position')} ({meta.get('geometryType') or 'unknown geometry'})")
        query = query_layer(args, 1, where="1=1", limit=1)
        error = query.get("error")
        if isinstance(error, dict):
            print(f"    ! query error: {error.get('message', 'query failed')}")
            return 1
        features = query.get("features")
        if not isinstance(features, list) or not features:
            print("    ! no current position returned")
            return 1
        first = features[0]
        if not isinstance(first, dict):
            print("    ! current position returned unexpected payload")
            return 1
        attrs = first.get("attributes")
        if isinstance(attrs, dict):
            current_attrs = attrs
        print("    current_position:")
        print_feature(first, meta.get("geometryType") if isinstance(meta.get("geometryType"), str) else None, label="latest")

    storm_name = current_attrs.get("STORMNAME")
    basin = current_attrs.get("BASIN")
    if not isinstance(storm_name, str) or not storm_name.strip():
        print("No storm name found in current position record.", file=sys.stderr)
        return 1
    basin_text = basin if isinstance(basin, str) and basin.strip() else None
    print(f"Storm: {storm_name}" + (f" ({basin_text})" if basin_text else ""))

    forecast_advdate: int | None = None
    if 0 in layer_ids:
        meta = layer_meta(0)
        print(f"- [0] {meta.get('name', 'Forecast Position')} ({meta.get('geometryType') or 'unknown geometry'})")
        latest_query = query_layer(
            args,
            0,
            where="1=1",
            limit=1,
        )
        error = latest_query.get("error")
        if isinstance(error, dict):
            print(f"    ! query error: {error.get('message', 'query failed')}")
        else:
            features = latest_query.get("features")
            if isinstance(features, list) and features and isinstance(features[0], dict):
                match = next(
                    (
                        feature
                        for feature in features
                        if matches_storm(feature, storm_name, basin_text)
                    ),
                    None,
                )
                if match is not None:
                    attrs = feature_attrs(match)
                    advdate = attrs.get("ADVDATE")
                    if isinstance(advdate, (int, float)):
                        forecast_advdate = int(advdate)
        query = query_layer(
            args,
            0,
            where="1=1",
            limit=max(1, args.limit * 10),
        )
        error = query.get("error")
        if isinstance(error, dict):
            print(f"    ! query error: {error.get('message', 'query failed')}")
        else:
            features = query.get("features")
            if isinstance(features, list):
                filtered = [
                    feature
                    for feature in features
                    if isinstance(feature, dict)
                    and matches_storm(feature, storm_name, basin_text)
                ]
                if forecast_advdate is not None:
                    filtered = [
                        feature
                        for feature in filtered
                        if feature_attrs(feature).get("ADVDATE") == forecast_advdate
                    ]
                print(f"    forecast_position_count: {len(filtered)}")
                for feature in filtered:
                    print_feature(feature, meta.get("geometryType") if isinstance(meta.get("geometryType"), str) else None)
            else:
                print("    ! query returned no features array")

    for layer_id in WIND_LAYER_IDS:
        if layer_id not in layer_ids:
            continue
        meta = layer_meta(layer_id)
        print(f"- [{layer_id}] {meta.get('name', f'layer-{layer_id}')} ({meta.get('geometryType') or 'unknown geometry'})")
        query = query_layer(
            args,
            layer_id,
            where="1=1",
            limit=max(1, args.limit * 10),
        )
        error = query.get("error")
        if isinstance(error, dict):
            print(f"    ! query error: {error.get('message', 'query failed')}")
            continue
        features = query.get("features")
        if not isinstance(features, list):
            print("    ! query returned no features array")
            continue
        filtered_features = wind_layer_filter(
            [feature for feature in features if isinstance(feature, dict)],
            storm_name,
            basin_text,
        )
        print(f"    wind_polygon_count: {len(filtered_features)}")
        for feature in filtered_features:
            print_feature(feature, meta.get("geometryType") if isinstance(meta.get("geometryType"), str) else None)
        if not filtered_features:
            print("    (no features returned)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
