from __future__ import annotations

import json
import math
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace
from typing import Any, cast

import pytest
from PySide6.QtCore import Qt

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
    OBSERVER_PRECIPITATION_MARKER_DISTANCE_M,
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
        PRECIPITATION_MIN_DISTANCE_KM <= sample.distance_km
        < PRECIPITATION_MAX_DISTANCE_KM
        for sample in samples
    )
    assert samples[0].azimuth_deg == pytest.approx(0.0)


def test_precipitation_request_samples_use_surrounding_sunflower_points() -> None:
    samples = generate_precipitation_request_samples(35.0, 139.0)
    assert len(samples) == 48
    assert samples == generate_precipitation_samples(35.0, 139.0)
    assert samples[0].distance_km == pytest.approx(32.0 * math.sqrt(0.5 / 48.0))


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
    assert precipitation_streak_height_deg(0.0) == pytest.approx(7.0)
    assert precipitation_streak_height_deg(16.0) == pytest.approx(4.3)
    assert precipitation_streak_height_deg(32.0) == pytest.approx(1.6)
    assert precipitation_distance_opacity_factor(0.0) == pytest.approx(1.0)
    assert precipitation_distance_opacity_factor(16.0) == pytest.approx(0.675)
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


def test_precipitation_renderer_draws_clipped_solid_tile_lines(monkeypatch) -> None:
    lines = []
    pens = []
    clip_rects = []

    class Painter:
        def save(self):
            pass

        def restore(self):
            pass

        def setClipRect(self, rect):
            clip_rects.append(rect)

        def setPen(self, pen):
            pens.append(pen)

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
        ScreenGeometry(center=(50, 50), radius=100),
        ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Test"),
        [column],
        opacity=0.5,
    )
    assert len(lines) == 2
    assert len(clip_rects) == 1
    assert clip_rects[0].bottom() == pytest.approx(20.0)
    assert all(point.y() <= 20.0 for line in lines for point in line)
    assert all(start.x() < end.x() and start.y() > end.y() for start, end in lines)
    assert len(pens) == 1
    pen = pens[0]
    assert pen.color().getRgb()[:3] == render_precipitation.PRECIPITATION_COLUMN_COLOR_RGB
    assert pen.color().alpha() == 106
    assert pen.widthF() == pytest.approx(
        render_precipitation._precipitation_tile_line_width(
            100.0 * precipitation_streak_height_deg(20.0) / 90.0, 1.0
        )
    )
    assert pen.capStyle() == Qt.PenCapStyle.FlatCap
    assert pen.style() == Qt.PenStyle.SolidLine


@pytest.mark.parametrize(
    "line_count, expected_offsets",
    [(1, (0.0,)), (2, (-0.5, 0.5)), (3, (-1.0, 0.0, 1.0))],
)
def test_precipitation_line_offsets_keep_even_counts_open(
    line_count: int, expected_offsets: tuple[float, ...]
) -> None:
    assert render_precipitation._precipitation_line_offsets(line_count) == expected_offsets


def test_precipitation_tile_lines_are_clipped_to_square() -> None:
    lines = render_precipitation._precipitation_tile_lines(
        50.0, 40.0, 20.0, 6, 1.0
    )
    assert len(lines) == 6
    for start, end in lines:
        for point in (start, end):
            assert 40.0 <= point.x() <= 60.0
            assert 30.0 <= point.y() <= 50.0


def test_precipitation_line_spacing_scales_with_tile_size() -> None:
    small = render_precipitation._precipitation_tile_line_spacing(8.0, 1.0)
    large = render_precipitation._precipitation_tile_line_spacing(20.0, 1.0)
    assert small == pytest.approx(2.0)
    assert large == pytest.approx(5.0)
    assert large / small == pytest.approx(20.0 / 8.0)


def test_precipitation_line_width_scales_with_tile_size() -> None:
    far = render_precipitation._precipitation_tile_line_width(5.0, 1.0)
    middle = render_precipitation._precipitation_tile_line_width(12.0, 1.0)
    near = render_precipitation._precipitation_tile_line_width(24.0, 1.0)
    assert far == pytest.approx(0.45)
    assert middle == pytest.approx(12.0 * 0.08)
    assert near == pytest.approx(24.0 * 0.08)
    assert far < middle < near


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

    assert not any(isinstance(item, ObserverPrecipitationMarker) for item in items)
    assert len(items) == 48
    assert sampled_coordinate_counts == [48]


def test_observer_precipitation_marker_is_centered_and_enlarged() -> None:
    lines = []
    pen_widths = []

    class Painter:
        def save(self):
            pass

        def restore(self):
            pass

        def setClipRect(self, rect):
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
    assert OBSERVER_PRECIPITATION_MARKER_DISTANCE_M == pytest.approx(50.0)
    assert pen_widths == [
        pytest.approx(
            render_precipitation._precipitation_tile_line_width(
                90.0 * 7.0 / 90.0 * OBSERVER_PRECIPITATION_MARKER_SCALE,
                1.0,
            )
        ),
    ]
    for start, end in lines:
        assert start.x() < end.x()
        assert start.y() > end.y()
        assert 120.0 - 5.0 <= start.x() <= 120.0 + 5.0
        assert 120.0 - 5.0 <= end.x() <= 120.0 + 5.0
        assert 90.0 - 5.0 <= start.y() <= 90.0 + 5.0
        assert 90.0 - 5.0 <= end.y() <= 90.0 + 5.0


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
