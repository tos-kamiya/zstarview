from datetime import datetime, timedelta, timezone

import pytest

from zstarview.render.meteors import meteor_age_label, meteor_age_opacity


@pytest.mark.parametrize(
    ("age", "expected"),
    [
        (timedelta(seconds=-1), 0.0),
        (timedelta(0), 1.0),
        (timedelta(hours=18), 1.0),
        (timedelta(hours=21), 0.5),
        (timedelta(hours=24), 0.0),
        (timedelta(hours=25), 0.0),
    ],
)
def test_meteor_age_opacity(age: timedelta, expected: float) -> None:
    display_time = datetime(2026, 8, 12, 12, tzinfo=timezone.utc)
    assert meteor_age_opacity(display_time - age, display_time) == pytest.approx(expected)


def test_meteor_age_opacity_accepts_naive_utc() -> None:
    display_time = datetime(2026, 8, 12, 12)
    assert meteor_age_opacity(display_time - timedelta(hours=21), display_time) == pytest.approx(0.5)


def test_meteor_age_opacity_can_use_delayed_window_end() -> None:
    display_time = datetime(2026, 8, 12, 12, tzinfo=timezone.utc)
    window_end = display_time - timedelta(hours=30)
    observation_time = window_end - timedelta(hours=21)

    assert meteor_age_opacity(observation_time, window_end) == pytest.approx(0.5)
    assert meteor_age_opacity(observation_time, display_time) == 0.0


def test_meteor_age_label_uses_signed_display_time_hours() -> None:
    display_time = datetime(2026, 8, 12, 12, tzinfo=timezone.utc)
    assert meteor_age_label(display_time - timedelta(hours=32), display_time) == "-32h"
    assert meteor_age_label(display_time + timedelta(hours=2), display_time) == "+2h"
