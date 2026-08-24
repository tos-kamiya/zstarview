from __future__ import annotations

from zstarview.utils.latlon_format import (
    LAT_LON_CACHE_DECIMALS,
    LAT_LON_DISPLAY_DECIMALS,
    format_lat_lon_cache_segment,
    format_lat_lon_display,
    format_lat_lon_value,
)


def test_latlon_display_format_uses_shared_precision() -> None:
    assert LAT_LON_DISPLAY_DECIMALS is None
    assert LAT_LON_CACHE_DECIMALS == 4
    assert format_lat_lon_value(35.4824704) == "35.4824704"
    assert format_lat_lon_value(133.0683567) == "133.0683567"
    assert format_lat_lon_display(35.4824704, 133.0683567) == "Lat: 35.4824704, Lon: 133.0683567"


def test_latlon_display_format_can_request_fixed_decimals() -> None:
    assert format_lat_lon_value(35.4824704, decimals=5) == "35.48247"


def test_latlon_cache_segment_uses_same_rounding() -> None:
    assert format_lat_lon_cache_segment(35.4824704) == "35p4825"
    assert format_lat_lon_cache_segment(-3.0666666666667) == "m3p0667"
