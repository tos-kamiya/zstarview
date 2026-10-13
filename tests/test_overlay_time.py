from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest

from zstarview.overlay_time import (
    classify_delta_t,
    classify_target_time,
    overlay_availability_for_delta,
    overlay_availability_for_time_mode,
    overlay_source_time_utc_from_delta,
    target_time_utc_from_delta,
)


def test_target_time_utc_from_delta_uses_utc_now() -> None:
    now_utc = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    assert target_time_utc_from_delta(timedelta(hours=-2), now_utc=now_utc) == datetime(
        2026, 3, 22, 10, 0, tzinfo=timezone.utc
    )


def test_classify_target_time_returns_present_within_tolerance() -> None:
    now_utc = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    assert (
        classify_target_time(
            now_utc + timedelta(minutes=4),
            now_utc=now_utc,
            present_tolerance_seconds=5 * 60,
        )
        == "present"
    )


def test_classify_delta_t_distinguishes_past_and_future() -> None:
    now_utc = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    assert classify_delta_t(timedelta(hours=-1), now_utc=now_utc) == "past"
    assert classify_delta_t(timedelta(hours=1), now_utc=now_utc) == "future"


def test_overlay_availability_for_time_mode_matches_policy() -> None:
    assert overlay_availability_for_time_mode("present").cloud is True
    assert overlay_availability_for_time_mode("present").aircraft is True
    assert overlay_availability_for_time_mode("present").satellite is True
    assert overlay_availability_for_time_mode("present").tropical_cyclone is True
    assert overlay_availability_for_time_mode("present").precipitation is True
    assert overlay_availability_for_time_mode("past").satellite is False
    assert overlay_availability_for_time_mode("past").cloud is False
    assert overlay_availability_for_time_mode("future").satellite is False
    assert overlay_availability_for_time_mode("future").tropical_cyclone is False
    assert overlay_availability_for_time_mode("past").precipitation is False
    assert overlay_availability_for_time_mode("future").precipitation is False


def test_overlay_availability_for_delta_matches_past_policy() -> None:
    now_utc = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    availability = overlay_availability_for_delta(timedelta(days=-1), now_utc=now_utc)
    assert availability.cloud is False
    assert availability.aircraft is False
    assert availability.satellite is False
    assert availability.tropical_cyclone is False
    assert availability.precipitation is False


@pytest.mark.parametrize("history_hours", [1, 24, 720])
@pytest.mark.parametrize("offset_seconds", [-1, 0, 1])
def test_past_live_layers_and_source_time_agree_at_retention_boundary(
    monkeypatch, history_hours, offset_seconds
) -> None:
    def unexpected_config_read():
        pytest.fail("Explicit retention must not read the user's settings")

    monkeypatch.setattr(
        "zstarview.overlay_time.overlay_history_seconds", unexpected_config_read
    )
    now = datetime(2026, 10, 13, 12, tzinfo=timezone.utc)
    retention = history_hours * 3600
    delta = timedelta(seconds=-(retention + offset_seconds))
    eligible = offset_seconds <= 0
    source_time = overlay_source_time_utc_from_delta(
        delta, now_utc=now, history_seconds=retention
    )
    assert source_time == (now if eligible else now + delta)
    availability = overlay_availability_for_delta(
        delta, now_utc=now, allow_cached_history=True, history_seconds=retention
    )
    assert availability.cloud is eligible
    assert availability.aircraft is eligible
    assert availability.satellite is eligible
    assert availability.tropical_cyclone is False
    assert availability.precipitation is False
