from datetime import datetime, timedelta, timezone

import pytest
from PySide6.QtCore import QPointF
from PySide6.QtGui import QColor

from zstarview.render.meteors import (
    METEOR_CORE_COLOR,
    METEOR_GLOW_COLOR,
    METEOR_LABEL_COLOR,
    _draw_meteor_trail_shape,
    meteor_age_label,
    meteor_age_opacity,
)


def test_meteor_uses_white_core_yellow_glow_and_intermediate_label() -> None:
    assert METEOR_CORE_COLOR == (255, 255, 255)
    assert METEOR_GLOW_COLOR == (255, 220, 120)
    assert METEOR_LABEL_COLOR == (255, 238, 188)


def test_meteor_trail_shape_peaks_at_four_fifths() -> None:
    class PainterProbe:
        def __init__(self) -> None:
            self.polygon = None

        def setPen(self, _pen: object) -> None:
            pass

        def setBrush(self, _brush: object) -> None:
            pass

        def drawPolygon(self, polygon: object) -> None:
            self.polygon = polygon

    painter = PainterProbe()
    _draw_meteor_trail_shape(
        painter,  # type: ignore[arg-type]
        QPointF(0.0, 0.0),
        QPointF(10.0, 0.0),
        color=QColor(255, 255, 255),
        start_half_width=0.4,
        peak_half_width=1.25,
        end_half_width=0.4,
    )

    assert painter.polygon is not None
    assert painter.polygon[1] == QPointF(8.0, 1.25)
    assert painter.polygon[3] == QPointF(8.0, -1.25)


@pytest.mark.parametrize(
    ("age", "expected"),
    [
        (timedelta(seconds=-1), 0.0),
        (timedelta(0), 1.0),
        (timedelta(hours=24), 1.0),
        (timedelta(hours=48), 1.0 - 0.7 * (24 / 72)),
        (timedelta(hours=72), 1.0 - 0.7 * (48 / 72)),
        (timedelta(hours=96), 0.3),
        (timedelta(hours=97), 0.3),
    ],
)
def test_meteor_age_opacity(age: timedelta, expected: float) -> None:
    display_time = datetime(2026, 8, 12, 12, tzinfo=timezone.utc)
    assert meteor_age_opacity(display_time - age, display_time) == pytest.approx(expected)


def test_meteor_age_opacity_accepts_naive_utc() -> None:
    display_time = datetime(2026, 8, 12, 12)
    assert meteor_age_opacity(display_time - timedelta(hours=96), display_time) == pytest.approx(0.3)


def test_meteor_age_opacity_uses_display_time_not_window_end() -> None:
    display_time = datetime(2026, 8, 12, 12, tzinfo=timezone.utc)
    window_end = display_time - timedelta(hours=30)
    observation_time = window_end - timedelta(hours=18)

    assert meteor_age_opacity(observation_time, display_time) == pytest.approx(
        1.0 - 0.7 * (24 / 72)
    )
    assert meteor_age_opacity(observation_time, window_end) == pytest.approx(1.0)


def test_meteor_age_label_uses_signed_display_time_hours() -> None:
    display_time = datetime(2026, 8, 12, 12, tzinfo=timezone.utc)
    assert meteor_age_label(display_time - timedelta(hours=32), display_time) == "-32h"
    assert meteor_age_label(display_time + timedelta(hours=2), display_time) == "+2h"
