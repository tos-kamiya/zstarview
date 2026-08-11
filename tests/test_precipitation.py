from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace
from typing import Any, cast

import pytest

from zstarview.cli.args import parse_args, parse_export_image_args
from zstarview.gui.window_updates import SkyWindowUpdatesMixin
from zstarview.precipitation import (
    PRECIPITATION_FAR_OPACITY_FACTOR,
    PRECIPITATION_MAX_DISTANCE_KM,
    PRECIPITATION_MIN_DISTANCE_KM,
    ProjectedPrecipitationColumn,
    fetch_open_meteo_precipitation,
    generate_precipitation_samples,
    parse_open_meteo_response,
    precipitation_cache_key,
    precipitation_distance_opacity_factor,
    precipitation_rate_mm_h,
    precipitation_streak_count,
    precipitation_streak_height_deg,
    precipitation_snapshot_is_fresh,
)
from zstarview.render import precipitation as render_precipitation
from zstarview.types import ViewerData


def test_precipitation_samples_are_deterministic_and_inside_annulus() -> None:
    samples = generate_precipitation_samples(35.0, 139.0)
    assert samples == generate_precipitation_samples(35.0, 139.0)
    assert len(samples) == 48
    assert all(
        PRECIPITATION_MIN_DISTANCE_KM < sample.distance_km
        < PRECIPITATION_MAX_DISTANCE_KM
        for sample in samples
    )
    assert samples[0].azimuth_deg == pytest.approx(0.0)


def test_precipitation_rate_normalizes_interval_amount() -> None:
    assert precipitation_rate_mm_h(1.0, 900) == pytest.approx(4.0)
    assert precipitation_rate_mm_h(1.0, 3600) == pytest.approx(1.0)
    assert precipitation_rate_mm_h(None, 900) is None
    with pytest.raises(ValueError, match="interval"):
        precipitation_rate_mm_h(1.0, 0)


def test_precipitation_streak_count_encodes_rate() -> None:
    assert precipitation_streak_count(0.0) == 0
    assert precipitation_streak_count(0.1) == 1
    assert precipitation_streak_count(1.0) == 1
    assert precipitation_streak_count(3.0) == 2
    assert precipitation_streak_count(1.0e9) == 6


def test_precipitation_streak_height_and_opacity_encode_distance() -> None:
    assert precipitation_streak_height_deg(8.0) == pytest.approx(16.0 / 3.0)
    assert precipitation_streak_height_deg(20.0) == pytest.approx(11.0 / 3.0)
    assert precipitation_streak_height_deg(32.0) == pytest.approx(2.0)
    assert precipitation_distance_opacity_factor(8.0) == pytest.approx(1.0)
    assert precipitation_distance_opacity_factor(20.0) == pytest.approx(0.675)
    assert precipitation_distance_opacity_factor(32.0) == pytest.approx(
        PRECIPITATION_FAR_OPACITY_FACTOR
    )


def _response_item(*, amount: float | None = 1.0, interval: int = 900) -> dict:
    return {
        "current_units": {
            "time": "iso8601",
            "interval": "seconds",
            "precipitation": "mm",
            "rain": "mm",
            "showers": "mm",
        },
        "current": {
            "time": "2026-08-11T12:15",
            "interval": interval,
            "precipitation": amount,
            "rain": amount,
            "showers": 0.0 if amount is not None else None,
        },
    }


def test_parse_open_meteo_response_keeps_missing_distinct_from_zero() -> None:
    samples = generate_precipitation_samples(35.0, 139.0)
    payload = [_response_item(amount=0.0) for _sample in samples]
    payload[1] = _response_item(amount=None)
    snapshot = parse_open_meteo_response(
        payload,
        samples,
        fetched_at_utc=datetime(2026, 8, 11, tzinfo=timezone.utc),
    )
    assert snapshot.values[0].rate_mm_h == 0.0
    assert snapshot.values[1].rate_mm_h is None


def test_parse_open_meteo_response_rejects_wrong_count_and_unit() -> None:
    samples = generate_precipitation_samples(35.0, 139.0)
    with pytest.raises(ValueError, match="location count"):
        parse_open_meteo_response([], samples, fetched_at_utc=datetime.now(timezone.utc))
    payload = [_response_item() for _sample in samples]
    payload[0]["current_units"]["precipitation"] = "inch"
    with pytest.raises(ValueError, match="unit"):
        parse_open_meteo_response(
            payload, samples, fetched_at_utc=datetime.now(timezone.utc)
        )


def test_fetch_open_meteo_uses_one_multiple_coordinate_request() -> None:
    samples = generate_precipitation_samples(35.0, 139.0)
    seen = {}

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return None

        def read(self) -> bytes:
            return json.dumps([_response_item() for _sample in samples]).encode()

    def opener(request, *, timeout):
        seen["url"] = request.full_url
        seen["timeout"] = timeout
        return Response()

    snapshot = fetch_open_meteo_precipitation(samples, opener=opener)
    assert len(snapshot.values) == 48
    assert seen["url"].count("latitude=") == 1
    assert "cell_selection=nearest" in seen["url"]
    assert "apikey" not in seen["url"]


def test_precipitation_cache_key_tracks_order_and_has_no_secret() -> None:
    samples = generate_precipitation_samples(35.0, 139.0)
    key = precipitation_cache_key(samples)
    assert key != precipitation_cache_key(tuple(reversed(samples)))
    assert "apikey" not in repr(key).casefold()


def test_precipitation_cache_freshness_uses_utc_and_rejects_future_time() -> None:
    samples = generate_precipitation_samples(35.0, 139.0)
    now = datetime(2026, 8, 11, 12, 0, tzinfo=timezone.utc)
    snapshot = parse_open_meteo_response(
        [_response_item() for _sample in samples],
        samples,
        fetched_at_utc=now,
    )
    assert precipitation_snapshot_is_fresh(
        snapshot, now_utc=now + timedelta(minutes=9, seconds=59)
    )
    assert not precipitation_snapshot_is_fresh(
        snapshot, now_utc=now + timedelta(minutes=10)
    )
    assert not precipitation_snapshot_is_fresh(
        snapshot, now_utc=now - timedelta(seconds=1)
    )


def test_precipitation_option_is_viewer_only() -> None:
    assert parse_args(["--precipitation-opacity", "0.6"]).precipitation_opacity == 0.6
    with pytest.raises(SystemExit):
        parse_export_image_args(["--precipitation-opacity", "0.6"])


def test_precipitation_renderer_draws_blue_rain_streaks(monkeypatch) -> None:
    lines = []

    class Painter:
        def save(self):
            pass

        def restore(self):
            pass

        def setPen(self, pen):
            assert pen.color().blue() == 255
            assert pen.color().alpha() == 86
            assert pen.widthF() == pytest.approx(1.8)

        def drawLine(self, start, end):
            lines.append((start, end))

    monkeypatch.setattr(render_precipitation, "is_in_fov", lambda *args, **kwargs: True)
    monkeypatch.setattr(
        render_precipitation,
        "altaz_to_normalized_xy",
        lambda alt, az, *args, **kwargs: (alt / 10.0, az / 10.0),
    )
    monkeypatch.setattr(
        render_precipitation,
        "normalized_to_screen_xy",
        lambda x, y, geometry: (x * 100.0, y * 100.0),
    )
    column = ProjectedPrecipitationColumn(1.0, 2.0, 3.0, 2.0, 20.0, 3.0)
    render_precipitation.draw_precipitation_columns(
        cast(Any, Painter()),
        object(),
        ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Test"),
        [column],
        opacity=0.5,
    )
    assert len(lines) == 2
    assert all(start.x() < end.x() for start, end in lines)


def test_precipitation_failure_removes_existing_columns() -> None:
    invalidations = []
    updates = []
    owner = SimpleNamespace(
        state=SimpleNamespace(
            precipitation_columns=[object()], precipitation_next_refresh_utc=None
        ),
        precipitation_status="ready",
        _compositor=SimpleNamespace(invalidate=lambda: invalidations.append(True)),
        request_client_update=lambda: updates.append(True),
    )

    SkyWindowUpdatesMixin._on_precipitation_failed(owner, {"error": "offline"})

    assert owner.state.precipitation_columns is None
    assert owner.precipitation_status == "unavailable"
    assert owner.state.precipitation_next_refresh_utc is not None
    assert invalidations == [True]
    assert updates == [True]
