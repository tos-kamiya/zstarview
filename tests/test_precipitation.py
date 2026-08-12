from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace
from typing import Any, cast

import pytest

from zstarview.cli import export_image as export_image_module
from zstarview.cli.args import parse_args, parse_export_image_args
from zstarview.gui.window_actions import SkyWindowActionsMixin
from zstarview.gui.window_updates import SkyWindowUpdatesMixin
from zstarview.precipitation import (
    OBSERVER_PRECIPITATION_MARKER_SCALE,
    PRECIPITATION_FAR_OPACITY_FACTOR,
    PRECIPITATION_MAX_DISTANCE_KM,
    PRECIPITATION_MIN_DISTANCE_KM,
    ObserverPrecipitationMarker,
    ProjectedPrecipitationColumn,
    fetch_open_meteo_precipitation,
    generate_precipitation_request_samples,
    generate_precipitation_samples,
    parse_open_meteo_response,
    precipitation_cache_key,
    precipitation_distance_opacity_factor,
    precipitation_rate_mm_h,
    precipitation_streak_count,
    precipitation_streak_height_deg,
    precipitation_snapshot_is_fresh,
    project_precipitation_columns,
)
from zstarview.render import precipitation as render_precipitation
from zstarview.types import ScreenGeometry, ViewerData


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


def test_precipitation_request_samples_put_observer_before_surroundings() -> None:
    samples = generate_precipitation_request_samples(35.0, 139.0)
    assert len(samples) == 49
    assert samples[0].latitude_deg == pytest.approx(35.0)
    assert samples[0].longitude_deg == pytest.approx(139.0)
    assert samples[0].distance_km == 0.0
    assert samples[1:] == generate_precipitation_samples(35.0, 139.0)


def test_precipitation_rate_normalizes_interval_amount() -> None:
    assert precipitation_rate_mm_h(1.0, 900) == pytest.approx(4.0)
    assert precipitation_rate_mm_h(1.0, 3600) == pytest.approx(1.0)
    assert precipitation_rate_mm_h(None, 900) is None
    with pytest.raises(ValueError, match="interval"):
        precipitation_rate_mm_h(1.0, 0)


def test_precipitation_streak_count_encodes_rate() -> None:
    assert precipitation_streak_count(0.0) == 0
    assert precipitation_streak_count(0.1) == 1
    assert precipitation_streak_count(2.0) == 1
    assert precipitation_streak_count(2.001) == 2
    assert precipitation_streak_count(6.0) == 2
    assert precipitation_streak_count(6.001) == 3
    assert precipitation_streak_count(14.0) == 3
    assert precipitation_streak_count(14.001) == 4
    assert precipitation_streak_count(30.0) == 4
    assert precipitation_streak_count(30.001) == 5
    assert precipitation_streak_count(62.0) == 5
    assert precipitation_streak_count(62.001) == 6
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
    samples = generate_precipitation_request_samples(35.0, 139.0)
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
    assert len(snapshot.values) == 49
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


def test_precipitation_option_is_available_to_viewer_and_export() -> None:
    assert parse_args(["--precipitation-opacity", "0.6"]).precipitation_opacity == 0.6
    assert parse_args(["-P", "0.6"]).precipitation_opacity == 0.6
    assert (
        parse_export_image_args(
            ["-P", "0.6", "--output", "test.png"]
        ).precipitation_opacity
        == 0.6
    )


def test_export_precipitation_requires_prior_terms_acceptance(
    monkeypatch, capsys
) -> None:
    monkeypatch.setattr(
        export_image_module,
        "open_meteo_noncommercial_terms_accepted",
        lambda: False,
    )

    export_image_module._require_open_meteo_consent_for_export(0.0)
    with pytest.raises(SystemExit) as exc_info:
        export_image_module._require_open_meteo_consent_for_export(0.5)

    assert exc_info.value.code == 1
    assert "Start zstarview or zstarview-gui" in capsys.readouterr().err


def test_export_precipitation_accepts_saved_terms(monkeypatch) -> None:
    monkeypatch.setattr(
        export_image_module,
        "open_meteo_noncommercial_terms_accepted",
        lambda: True,
    )

    export_image_module._require_open_meteo_consent_for_export(0.5)


def test_precipitation_menu_toggle_preserves_configured_opacity() -> None:
    updates: list[str] = []
    invalidations: list[bool] = []
    action = SimpleNamespace(
        checked=True,
        isChecked=lambda: action.checked,
        setChecked=lambda value: setattr(action, "checked", bool(value)),
    )
    owner = SimpleNamespace(
        _precipitation_toggle_supported=True,
        _precipitation_opacity_when_enabled=0.6,
        precipitation_opacity=0.6,
        _action_toggle_precipitation=action,
        _compositor=SimpleNamespace(
            invalidate=lambda: invalidations.append(True)
        ),
        start_background_precipitation_update=lambda **kwargs: (
            updates.append(str(kwargs["reason"])) or True
        ),
        request_client_update=lambda: updates.append("repaint"),
    )

    SkyWindowActionsMixin.toggle_precipitation(owner)
    assert owner.precipitation_opacity == 0.0
    assert action.checked is False

    SkyWindowActionsMixin.toggle_precipitation(owner)
    assert owner.precipitation_opacity == pytest.approx(0.6)
    assert action.checked is True
    assert updates == ["repaint", "toggle-on", "repaint"]
    assert invalidations == [True, True]


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


def test_precipitation_projection_keeps_observer_out_of_altaz_projection(
    monkeypatch,
) -> None:
    samples = generate_precipitation_request_samples(35.0, 139.0)
    payload = [_response_item(amount=0.5) for _sample in samples]
    snapshot = parse_open_meteo_response(
        payload,
        samples,
        fetched_at_utc=datetime(2026, 8, 11, tzinfo=timezone.utc),
    )
    sampled_coordinate_counts = []

    def ground_sampler(*args, **kwargs):
        def sample(latitudes, longitudes):
            sampled_coordinate_counts.append(len(latitudes))
            return [0.0] * len(latitudes)

        return sample

    monkeypatch.setattr(
        "zstarview.precipitation.build_road_night_light_ground_sampler",
        ground_sampler,
    )
    items = project_precipitation_columns(
        snapshot,
        ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Test"),
    )

    assert isinstance(items[-1], ObserverPrecipitationMarker)
    assert items[-1].rate_mm_h == pytest.approx(2.0)
    assert len(items) == 49
    assert sampled_coordinate_counts == [48, 1]


def test_observer_precipitation_marker_is_centered_and_enlarged() -> None:
    lines = []
    pen_widths = []

    class Painter:
        def save(self):
            pass

        def restore(self):
            pass

        def setPen(self, pen):
            pen_widths.append(pen.widthF())

        def drawLine(self, start, end):
            lines.append((start, end))

    geometry = ScreenGeometry(center=(120, 90), radius=90)
    render_precipitation.draw_precipitation_columns(
        cast(Any, Painter()),
        geometry,
        ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Test"),
        [ObserverPrecipitationMarker(rate_mm_h=3.0)],
        opacity=0.5,
    )

    assert len(lines) == 2
    assert pen_widths == [pytest.approx(1.8 * OBSERVER_PRECIPITATION_MARKER_SCALE)]
    for start, end in lines:
        assert (start.y() + end.y()) / 2.0 == pytest.approx(90.0)
        assert abs(start.y() - end.y()) == pytest.approx(
            90.0
            * (16.0 / 3.0)
            / 90.0
            * OBSERVER_PRECIPITATION_MARKER_SCALE
        )
    marker_center_x = sum(
        (start.x() + end.x()) / 2.0 for start, end in lines
    ) / len(lines)
    assert marker_center_x == pytest.approx(120.0)


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
