from __future__ import annotations

import re

from .latlon_format import format_lat_lon_display

_LAT_LON_DISPLAY_PATTERN = re.compile(
    r"Lat:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)),\s*"
    r"Lon:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+))\Z"
)


def _is_same_lat_lon_display(text: str, lat_deg: float, lon_deg: float) -> bool:
    if text == format_lat_lon_display(lat_deg, lon_deg):
        return True
    match = _LAT_LON_DISPLAY_PATTERN.fullmatch(text)
    if match is None:
        return False
    return float(match.group(1)) == float(lat_deg) and float(match.group(2)) == float(lon_deg)


def format_height_m(value_m: float) -> str:
    rounded = round(float(value_m))
    if abs(float(value_m) - rounded) < 0.05:
        return f"{int(rounded)} m"
    return f"{float(value_m):.1f} m"


def build_location_info_lines(
    display_name: str,
    lat_deg: float,
    lon_deg: float,
    *,
    ground_elevation_m: float = 0.0,
    location_height_m: float = 0.0,
    height_add_m: float = 1.7,
) -> list[str]:
    lines = []
    lat_lon_text = format_lat_lon_display(lat_deg, lon_deg)
    ground = max(0.0, float(ground_elevation_m))
    structure = max(0.0, float(location_height_m))
    height_add = max(0.0, float(height_add_m))
    display_text = str(display_name).strip()
    if display_text and not _is_same_lat_lon_display(display_text, lat_deg, lon_deg):
        lines.append(display_text)
    height_parts = [f"ground {format_height_m(ground)}"]
    if structure > 0.0:
        height_parts.append(f"building {format_height_m(structure)}")
    height_parts.append(f"add {format_height_m(height_add)}")
    lines.append(lat_lon_text)
    lines.append(f"Height: {', '.join(height_parts)}")
    return lines


def format_location_summary(
    display_name: str,
    lat_deg: float,
    lon_deg: float,
    *,
    ground_elevation_m: float = 0.0,
    location_height_m: float = 0.0,
    height_add_m: float = 1.7,
) -> str:
    return " | ".join(
        build_location_info_lines(
            display_name,
            lat_deg,
            lon_deg,
            ground_elevation_m=ground_elevation_m,
            location_height_m=location_height_m,
            height_add_m=height_add_m,
        )
    )
