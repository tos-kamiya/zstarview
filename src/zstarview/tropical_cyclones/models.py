from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
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
