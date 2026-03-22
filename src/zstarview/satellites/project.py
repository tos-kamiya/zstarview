from __future__ import annotations

from collections.abc import Mapping, Sequence

import astropy.time
from skyfield.api import Topos
import skyfield.api

from .fetch import build_earth_satellites
from .types import SatelliteOmmRecord, SatelliteOverlayPoint

_DEFAULT_GROUP_ORDER = ("station", "starlink")
_MAX_MARKERS_BY_GROUP = {
    "station": 8,
    "starlink": 20,
}
_MARKER_SCALE_BY_GROUP = {
    "station": 0.3,
    "starlink": 0.156,
}


def project_satellite_records(
    records_by_group: Mapping[str, Sequence[SatelliteOmmRecord]],
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    time_obj: astropy.time.Time,
) -> list[SatelliteOverlayPoint]:
    ts = skyfield.api.load.timescale()
    t = ts.from_astropy(time_obj)
    observer = Topos(
        latitude_degrees=float(observer_lat),
        longitude_degrees=float(observer_lon),
        elevation_m=float(observer_height_m),
    )

    overlay_points: list[SatelliteOverlayPoint] = []
    for group_key in _iter_group_order(records_by_group):
        records = records_by_group.get(group_key) or ()
        satellites = build_earth_satellites(records, ts=ts)
        visible_points: list[SatelliteOverlayPoint] = []
        for satellite in satellites:
            topocentric = (satellite - observer).at(t)
            alt, az, _distance = topocentric.altaz()
            alt_deg = float(alt.degrees)
            if alt_deg <= 0.0:
                continue
            visible_points.append(
                SatelliteOverlayPoint(
                    group_key=str(group_key),
                    satellite_name=str(getattr(satellite, "name", "") or group_key.upper()),
                    alt_deg=alt_deg,
                    az_deg=float(az.degrees),
                    marker_scale=float(_MARKER_SCALE_BY_GROUP.get(group_key, 0.13)),
                    show_label=group_key == "station",
                )
            )
        visible_points.sort(key=lambda point: float(point.alt_deg), reverse=True)
        overlay_points.extend(visible_points[: _MAX_MARKERS_BY_GROUP.get(group_key, len(visible_points))])
    return overlay_points


def _iter_group_order(records_by_group: Mapping[str, Sequence[SatelliteOmmRecord]]) -> list[str]:
    ordered = [group_key for group_key in _DEFAULT_GROUP_ORDER if group_key in records_by_group]
    for group_key in records_by_group:
        if group_key not in ordered:
            ordered.append(str(group_key))
    return ordered
