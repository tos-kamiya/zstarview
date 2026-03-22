from __future__ import annotations

from datetime import datetime, timezone

from zstarview.startup import _startup_parse_time_arguments


class _FixedDateTime(datetime):
    @classmethod
    def now(cls, tz=None):
        return cls(2026, 3, 22, 0, 0, 0, tzinfo=tz or timezone.utc)


def test_startup_parse_time_arguments_uses_location_timezone_when_tz_omitted(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.startup.datetime", _FixedDateTime)

    delta = _startup_parse_time_arguments(
        "2026-03-22 09:00",
        0,
        0,
        timezone_name="Asia/Tokyo",
    )

    assert delta.total_seconds() == 0


def test_startup_parse_time_arguments_prefers_timezone_override(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.startup.datetime", _FixedDateTime)

    delta = _startup_parse_time_arguments(
        "2026-03-22 09:00",
        0,
        0,
        timezone_name="Asia/Tokyo",
        timezone_override="UTC",
    )

    assert delta.total_seconds() == 9 * 3600
