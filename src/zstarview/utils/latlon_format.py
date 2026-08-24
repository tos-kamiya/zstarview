from __future__ import annotations

LAT_LON_DISPLAY_DECIMALS = None
LAT_LON_CACHE_DECIMALS = 4


def format_lat_lon_value(value: float, *, decimals: int | None = LAT_LON_DISPLAY_DECIMALS) -> str:
    if decimals is None:
        return str(float(value))
    return f"{float(value):.{int(decimals)}f}"


def format_lat_lon_display(
    lat_deg: float,
    lon_deg: float,
    *,
    decimals: int | None = LAT_LON_DISPLAY_DECIMALS,
) -> str:
    return (
        f"Lat: {format_lat_lon_value(lat_deg, decimals=decimals)}, "
        f"Lon: {format_lat_lon_value(lon_deg, decimals=decimals)}"
    )


def format_lat_lon_cache_segment(value: float, *, decimals: int = LAT_LON_CACHE_DECIMALS) -> str:
    return format_lat_lon_value(value, decimals=decimals).replace("-", "m").replace(".", "p")
