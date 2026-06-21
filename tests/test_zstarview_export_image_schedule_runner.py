from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path
import sys

from zoneinfo import ZoneInfo

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))

import zstarview_export_image_schedule_runner as scheduler  # noqa: E402


def test_parse_job_line_accepts_bare_utc_offset_and_repeat() -> None:
    job = scheduler._parse_job_line(
        '05:30:00 +08 x3 zstarview-export-image "@1.28815;103.85905" -o screenshot-%t.png',
        1,
    )
    assert job is not None
    assert job.hour == 5
    assert job.minute == 30
    assert job.second == 0
    assert job.repeat_count == 3
    assert job.command[0] == "zstarview-export-image"
    assert job.command[-1] == "screenshot-%t.png"


def test_nonexistent_dst_time_is_skipped_with_warning(caplog) -> None:
    job = scheduler._parse_job_line("02:30:00 America/New_York echo run", 1)
    assert job is not None

    now_utc = datetime(2026, 3, 8, 0, 0, tzinfo=timezone.utc)
    with caplog.at_level("WARNING"):
        occurrence = scheduler._first_valid_occurrence(0, job, now_utc)

    assert occurrence.scheduled_utc == datetime(2026, 3, 9, 6, 30, tzinfo=timezone.utc)
    assert "skipping nonexistent local time 2026-03-08 02:30:00" in caplog.text


def test_ambiguous_dst_time_uses_first_occurrence() -> None:
    job = scheduler._parse_job_line("01:30:00 America/New_York echo run", 1)
    assert job is not None

    now_utc = datetime(2026, 11, 1, 4, 0, tzinfo=timezone.utc)
    occurrence = scheduler._first_valid_occurrence(0, job, now_utc)

    assert occurrence.scheduled_utc == datetime(2026, 11, 1, 5, 30, tzinfo=timezone.utc)
    assert occurrence.scheduled_local.tzinfo == ZoneInfo("America/New_York")
    assert occurrence.scheduled_local.fold == 0
    assert occurrence.scheduled_local.utcoffset() == timedelta(hours=-4)


def test_expand_placeholders_uses_actual_start_time() -> None:
    command = ("zstarview-export-image", "-o", "shot-%t.png")
    expanded = scheduler._expand_placeholders(
        command,
        start_time_utc=datetime(2026, 6, 21, 12, 34, 56, tzinfo=timezone.utc),
    )
    assert expanded[-1] == "shot-20260621T123456Z.png"


def test_next_occurrence_advances_by_five_minutes() -> None:
    job = scheduler._parse_job_line("05:30:00 UTC x3 echo run", 1)
    assert job is not None

    first = scheduler._first_valid_occurrence(0, job, datetime(2026, 6, 21, 0, 0, tzinfo=timezone.utc))
    second = scheduler._next_occurrence(first)
    third = scheduler._next_occurrence(second)

    assert first.scheduled_utc == datetime(2026, 6, 21, 5, 30, tzinfo=timezone.utc)
    assert second.scheduled_utc == datetime(2026, 6, 21, 5, 35, tzinfo=timezone.utc)
    assert third.scheduled_utc == datetime(2026, 6, 21, 5, 40, tzinfo=timezone.utc)
