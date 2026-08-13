from datetime import datetime, timedelta, timezone

import pytest

from zstarview.render.meteors import METEOR_TRAIL_COLOR, meteor_age_label, meteor_age_opacity


def test_meteor_color_is_slightly_whiter_green() -> None:
    assert METEOR_TRAIL_COLOR == (230, 245, 205)


@pytest.mark.parametrize(
    ("age", "expected"),
    [
        (timedelta(seconds=-1), 0.0),
        (timedelta(0), 1.0),
        (timedelta(hours=24), 1.0 - 0.7 * (24 / 72)),
        (timedelta(hours=48), 1.0 - 0.7 * (48 / 72)),
        (timedelta(hours=72), 0.3),
        (timedelta(hours=96), 0.3),
        (timedelta(hours=97), 0.3),
    ],
)
def test_meteor_age_opacity(age: timedelta, expected: float) -> None:
    display_time = datetime(2026, 8, 12, 12, tzinfo=timezone.utc)
    assert meteor_age_opacity(display_time - age, display_time) == pytest.approx(expected)


def test_meteor_age_opacity_accepts_naive_utc() -> None:
    display_time = datetime(2026, 8, 12, 12)
    assert meteor_age_opacity(display_time - timedelta(hours=72), display_time) == pytest.approx(0.3)


def test_meteor_age_opacity_uses_display_time_not_window_end() -> None:
    display_time = datetime(2026, 8, 12, 12, tzinfo=timezone.utc)
    window_end = display_time - timedelta(hours=30)
    observation_time = window_end - timedelta(hours=18)

    assert meteor_age_opacity(observation_time, display_time) == pytest.approx(
        1.0 - 0.7 * (48 / 72)
    )
    assert meteor_age_opacity(observation_time, window_end) == pytest.approx(
        1.0 - 0.7 * (18 / 72)
    )


def test_meteor_age_label_uses_signed_display_time_hours() -> None:
    display_time = datetime(2026, 8, 12, 12, tzinfo=timezone.utc)
    assert meteor_age_label(display_time - timedelta(hours=32), display_time) == "-32h"
    assert meteor_age_label(display_time + timedelta(hours=2), display_time) == "+2h"
