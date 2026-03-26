from zstarview.aircraft_constants import (
    AIRCRAFT_BBOX_AREA_SQ_DEG,
    AIRCRAFT_BBOX_DELTA_DEG,
    AIRCRAFT_REFRESH_INTERVAL_SECONDS,
    AIRCRAFT_TRAIL_HALF_SPAN_SECONDS,
)


def test_aircraft_refresh_interval_is_five_minutes() -> None:
    assert AIRCRAFT_REFRESH_INTERVAL_SECONDS == 300


def test_aircraft_bbox_is_one_degree_each_side() -> None:
    assert AIRCRAFT_BBOX_DELTA_DEG == 1.0


def test_aircraft_bbox_area_stays_within_one_credit_band() -> None:
    assert AIRCRAFT_BBOX_AREA_SQ_DEG == 4.0
    assert AIRCRAFT_BBOX_AREA_SQ_DEG <= 25.0


def test_aircraft_trail_half_span_is_four_seconds() -> None:
    assert AIRCRAFT_TRAIL_HALF_SPAN_SECONDS == 4.0
