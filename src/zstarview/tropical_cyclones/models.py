from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from typing import Any


def _dt_to_iso(value: datetime | None) -> str | None:
    if value is None:
        return None
    return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def _dt_from_iso(value: object) -> datetime | None:
    if not isinstance(value, str) or not value.strip():
        return None
    text = value.strip().replace("Z", "+00:00")
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _wrap_lon_deg(value: float) -> float:
    wrapped = (float(value) + 180.0) % 360.0 - 180.0
    if wrapped == -180.0:
        return 180.0
    return wrapped


def _wrap_lon_delta_deg(delta_deg: float) -> float:
    return _wrap_lon_deg(float(delta_deg))


def _lerp(value0: float, value1: float, ratio: float) -> float:
    return float(value0) + (float(value1) - float(value0)) * float(ratio)


@dataclass(frozen=True, slots=True)
class TropicalCyclonePoint:
    lat_deg: float
    lon_deg: float
    valid_time_utc: datetime | None = None
    label: str | None = None
    tau_hr: int | None = None
    maxwind_kt: float | None = None
    gust_kt: float | None = None
    mslp_hpa: float | None = None
    dvlbl: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "lat_deg": float(self.lat_deg),
            "lon_deg": float(self.lon_deg),
            "valid_time_utc": _dt_to_iso(self.valid_time_utc),
            "label": self.label,
            "tau_hr": self.tau_hr,
            "maxwind_kt": self.maxwind_kt,
            "gust_kt": self.gust_kt,
            "mslp_hpa": self.mslp_hpa,
            "dvlbl": self.dvlbl,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> TropicalCyclonePoint | None:
        try:
            lat_deg = float(data["lat_deg"])
            lon_deg = float(data["lon_deg"])
        except Exception:
            return None
        tau_hr = data.get("tau_hr")
        maxwind_kt = data.get("maxwind_kt")
        gust_kt = data.get("gust_kt")
        mslp_hpa = data.get("mslp_hpa")
        return cls(
            lat_deg=lat_deg,
            lon_deg=lon_deg,
            valid_time_utc=_dt_from_iso(data.get("valid_time_utc")),
            label=data.get("label") if isinstance(data.get("label"), str) else None,
            tau_hr=int(tau_hr) if isinstance(tau_hr, (int, float)) else None,
            maxwind_kt=float(maxwind_kt) if isinstance(maxwind_kt, (int, float)) else None,
            gust_kt=float(gust_kt) if isinstance(gust_kt, (int, float)) else None,
            mslp_hpa=float(mslp_hpa) if isinstance(mslp_hpa, (int, float)) else None,
            dvlbl=data.get("dvlbl") if isinstance(data.get("dvlbl"), str) else None,
        )


@dataclass(frozen=True, slots=True)
class TropicalCyclonePolygon:
    layer_id: int
    name: str
    rings: tuple[tuple[tuple[float, float], ...], ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "layer_id": int(self.layer_id),
            "name": self.name,
            "rings": [
                [[float(lat), float(lon)] for lat, lon in ring]
                for ring in self.rings
            ],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> TropicalCyclonePolygon | None:
        try:
            layer_id = int(data["layer_id"])
            name = str(data["name"])
            rings_data = data["rings"]
        except Exception:
            return None
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
                    coords.append((float(point[0]), float(point[1])))
            if coords:
                rings.append(tuple(coords))
        return cls(layer_id=layer_id, name=name, rings=tuple(rings))


@dataclass(frozen=True, slots=True)
class TropicalCycloneSnapshot:
    storm_name: str
    basin: str | None
    advdate_utc: datetime | None
    observed_position: TropicalCyclonePoint
    forecast_positions: tuple[TropicalCyclonePoint, ...] = field(default_factory=tuple)
    wind_polygons: tuple[TropicalCyclonePolygon, ...] = field(default_factory=tuple)
    source_url: str = ""
    service_name: str = ""
    refreshed_at_utc: datetime | None = None
    current_storm_id: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "storm_name": self.storm_name,
            "basin": self.basin,
            "advdate_utc": _dt_to_iso(self.advdate_utc),
            "observed_position": self.observed_position.to_dict(),
            "forecast_positions": [point.to_dict() for point in self.forecast_positions],
            "wind_polygons": [polygon.to_dict() for polygon in self.wind_polygons],
            "source_url": self.source_url,
            "service_name": self.service_name,
            "refreshed_at_utc": _dt_to_iso(self.refreshed_at_utc),
            "current_storm_id": self.current_storm_id,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> TropicalCycloneSnapshot | None:
        storm_name = data.get("storm_name")
        if not isinstance(storm_name, str) or not storm_name.strip():
            return None
        observed = TropicalCyclonePoint.from_dict(data.get("observed_position") or {})
        if observed is None:
            return None
        forecast_positions_data = data.get("forecast_positions")
        forecast_positions: list[TropicalCyclonePoint] = []
        if isinstance(forecast_positions_data, list):
            for item in forecast_positions_data:
                if isinstance(item, dict):
                    point = TropicalCyclonePoint.from_dict(item)
                    if point is not None:
                        forecast_positions.append(point)
        wind_polygons_data = data.get("wind_polygons")
        wind_polygons: list[TropicalCyclonePolygon] = []
        if isinstance(wind_polygons_data, list):
            for item in wind_polygons_data:
                if isinstance(item, dict):
                    polygon = TropicalCyclonePolygon.from_dict(item)
                    if polygon is not None:
                        wind_polygons.append(polygon)
        basin = data.get("basin")
        service_name = data.get("service_name")
        source_url = data.get("source_url")
        current_storm_id = data.get("current_storm_id")
        return cls(
            storm_name=storm_name,
            basin=basin if isinstance(basin, str) and basin.strip() else None,
            advdate_utc=_dt_from_iso(data.get("advdate_utc")),
            observed_position=observed,
            forecast_positions=tuple(forecast_positions),
            wind_polygons=tuple(wind_polygons),
            source_url=source_url if isinstance(source_url, str) else "",
            service_name=service_name if isinstance(service_name, str) else "",
            refreshed_at_utc=_dt_from_iso(data.get("refreshed_at_utc")),
            current_storm_id=current_storm_id if isinstance(current_storm_id, str) else None,
        )

    def has_projectable_timeline(self) -> bool:
        if self.advdate_utc is None and self.refreshed_at_utc is None:
            return False
        if self.forecast_positions:
            for point in self.forecast_positions:
                if point.valid_time_utc is not None:
                    return True
                if point.tau_hr is not None:
                    return True
        return False


def _point_time_utc(
    snapshot: TropicalCycloneSnapshot,
    point: TropicalCyclonePoint,
) -> datetime | None:
    if point.valid_time_utc is not None:
        return point.valid_time_utc
    if point.tau_hr is not None and snapshot.advdate_utc is not None:
        return snapshot.advdate_utc + timedelta(hours=float(point.tau_hr))
    return snapshot.advdate_utc or snapshot.refreshed_at_utc


def _interpolate_point(
    start: TropicalCyclonePoint,
    end: TropicalCyclonePoint,
    ratio: float,
    *,
    valid_time_utc: datetime | None = None,
) -> TropicalCyclonePoint:
    lon_delta = _wrap_lon_delta_deg(end.lon_deg - start.lon_deg)
    lon_deg = _wrap_lon_deg(start.lon_deg + lon_delta * float(ratio))
    return TropicalCyclonePoint(
        lat_deg=_lerp(start.lat_deg, end.lat_deg, ratio),
        lon_deg=lon_deg,
        valid_time_utc=valid_time_utc or end.valid_time_utc or start.valid_time_utc,
        label=start.label or end.label,
        tau_hr=start.tau_hr if ratio <= 0.5 else end.tau_hr,
        maxwind_kt=_lerp(start.maxwind_kt, end.maxwind_kt, ratio)
        if start.maxwind_kt is not None and end.maxwind_kt is not None
        else start.maxwind_kt or end.maxwind_kt,
        gust_kt=_lerp(start.gust_kt, end.gust_kt, ratio)
        if start.gust_kt is not None and end.gust_kt is not None
        else start.gust_kt or end.gust_kt,
        mslp_hpa=_lerp(start.mslp_hpa, end.mslp_hpa, ratio)
        if start.mslp_hpa is not None and end.mslp_hpa is not None
        else start.mslp_hpa or end.mslp_hpa,
        dvlbl=start.dvlbl or end.dvlbl,
    )


def _translate_polygon(
    polygon: TropicalCyclonePolygon,
    *,
    lat_delta_deg: float,
    lon_delta_deg: float,
) -> TropicalCyclonePolygon:
    rings: list[tuple[tuple[float, float], ...]] = []
    for ring in polygon.rings:
        translated_ring: list[tuple[float, float]] = []
        for lat_deg, lon_deg in ring:
            translated_ring.append(
                (
                    float(lat_deg) + float(lat_delta_deg),
                    _wrap_lon_deg(float(lon_deg) + float(lon_delta_deg)),
                )
            )
        if translated_ring:
            rings.append(tuple(translated_ring))
    return TropicalCyclonePolygon(layer_id=polygon.layer_id, name=polygon.name, rings=tuple(rings))


def project_tropical_cyclone_snapshot(
    snapshot: TropicalCycloneSnapshot,
    when_utc: datetime,
) -> TropicalCycloneSnapshot:
    if when_utc.tzinfo is None:
        when_utc = when_utc.replace(tzinfo=timezone.utc)
    when_utc = when_utc.astimezone(timezone.utc)

    timeline: list[tuple[datetime, TropicalCyclonePoint]] = []
    for point in (snapshot.observed_position, *snapshot.forecast_positions):
        point_time = _point_time_utc(snapshot, point)
        if point_time is not None:
            timeline.append((point_time, point))
    if not timeline:
        return snapshot
    timeline.sort(key=lambda item: item[0])
    if len(timeline) == 1:
        projected_position = timeline[0][1]
    else:
        projected_position = timeline[0][1]
        for (start_time, start_point), (end_time, end_point) in zip(timeline, timeline[1:]):
            if when_utc <= start_time:
                projected_position = start_point
                break
            if start_time <= when_utc <= end_time:
                span_seconds = max(1.0, (end_time - start_time).total_seconds())
                ratio = (when_utc - start_time).total_seconds() / span_seconds
                projected_position = _interpolate_point(
                    start_point,
                    end_point,
                    ratio,
                    valid_time_utc=when_utc,
                )
                break
        else:
            projected_position = timeline[-1][1]

    observed = snapshot.observed_position
    lat_delta = float(projected_position.lat_deg) - float(observed.lat_deg)
    lon_delta = _wrap_lon_delta_deg(float(projected_position.lon_deg) - float(observed.lon_deg))
    translated_polygons = tuple(
        _translate_polygon(
            polygon,
            lat_delta_deg=lat_delta,
            lon_delta_deg=lon_delta,
        )
        for polygon in snapshot.wind_polygons
    )
    return TropicalCycloneSnapshot(
        storm_name=snapshot.storm_name,
        basin=snapshot.basin,
        advdate_utc=snapshot.advdate_utc,
        observed_position=projected_position,
        forecast_positions=snapshot.forecast_positions,
        wind_polygons=translated_polygons,
        source_url=snapshot.source_url,
        service_name=snapshot.service_name,
        refreshed_at_utc=when_utc,
        current_storm_id=snapshot.current_storm_id,
    )
