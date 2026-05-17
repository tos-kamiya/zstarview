from __future__ import annotations

from .latlon_format import format_lat_lon_display


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
    if display_text and display_text != lat_lon_text:
        lines.append(display_text)
    lines.append(
        f"{lat_lon_text} | Ground: {format_height_m(ground)}, Building: {format_height_m(structure)} | Height add: {format_height_m(height_add)}"
    )
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
