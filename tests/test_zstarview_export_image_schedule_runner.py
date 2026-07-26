from __future__ import annotations

from datetime import datetime, timedelta, timezone
from io import StringIO
from pathlib import Path
import sys
import threading
import time

import pytest
from zoneinfo import ZoneInfo

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))

import zstarview_export_image_schedule_runner as scheduler  # noqa: E402


class _TtyBuffer:
    def __init__(self) -> None:
        self.parts: list[str] = []

    def write(self, text: str) -> int:
        self.parts.append(text)
        return len(text)

    def flush(self) -> None:
        pass

    def isatty(self) -> bool:
        return True

    def getvalue(self) -> str:
        return "".join(self.parts)


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


def test_schedule_listing_orders_tasks_and_lists_repeats_once() -> None:
    jobs = [
        scheduler._parse_job_line(
            "06:00:00 UTC x8 zstarview-export-image -o late-%t.png", 1
        ),
        scheduler._parse_job_line(
            "05:30:00 UTC x4 zstarview-export-image -p 'Circular Quay' -A5 -Z90 "
            "-o screenshot-sydney-%t.png",
            2,
        ),
    ]
    assert all(job is not None for job in jobs)
    stream = StringIO()

    scheduler._print_schedule(
        jobs, datetime(2026, 6, 21, 5, 0, tzinfo=timezone.utc), stream  # type: ignore[arg-type]
    )

    lines = stream.getvalue().splitlines()
    assert lines == [
        "+00:30:00 Task 2 -> zstarview-export-image -p 'Circular Quay' -A5 -Z90 -o screenshot-sydney-%t.png",
        "+01:00:00 Task 1 -> zstarview-export-image -o late-%t.png",
    ]
    assert "残り時間" not in stream.getvalue()
    assert "[1/4]" not in stream.getvalue()
    assert "[1/8]" not in stream.getvalue()


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


def test_next_occurrence_advances_by_six_minutes() -> None:
    job = scheduler._parse_job_line("05:30:00 UTC x3 echo run", 1)
    assert job is not None

    first = scheduler._first_valid_occurrence(0, job, datetime(2026, 6, 21, 0, 0, tzinfo=timezone.utc))
    second = scheduler._next_occurrence(first)
    third = scheduler._next_occurrence(second)

    assert first.scheduled_utc == datetime(2026, 6, 21, 5, 30, tzinfo=timezone.utc)
    assert second.scheduled_utc == datetime(2026, 6, 21, 5, 36, tzinfo=timezone.utc)
    assert third.scheduled_utc == datetime(2026, 6, 21, 5, 42, tzinfo=timezone.utc)


def test_parse_job_line_strips_inline_comment() -> None:
    job = scheduler._parse_job_line(
        "05:30:00 UTC echo run # this is an inline comment",
        1,
    )
    assert job is not None
    assert job.hour == 5
    assert job.minute == 30
    assert job.second == 0
    assert job.timezone_text == "UTC"
    assert job.command == ("echo", "run")


def test_parse_job_line_ignores_full_line_comment_and_blank_line() -> None:
    assert scheduler._parse_job_line("# full-line comment", 1) is None
    assert scheduler._parse_job_line("   # indented comment", 2) is None
    assert scheduler._parse_job_line("   ", 3) is None
    assert scheduler._parse_job_line("\t", 4) is None


def test_parse_job_line_inline_comment_without_preceding_space() -> None:
    job = scheduler._parse_job_line(
        '05:30:00 UTC echo run#attached-comment',
        1,
    )
    assert job is not None
    assert job.command == ("echo", "run")


def test_repeat_label_for_single_run_is_empty() -> None:
    job = scheduler._parse_job_line("05:30:00 UTC echo run", 1)
    assert job is not None
    occurrence = scheduler._first_valid_occurrence(
        0, job, datetime(2026, 6, 21, 0, 0, tzinfo=timezone.utc)
    )
    assert occurrence.repeat_label == ""


def test_repeat_label_for_repeated_run() -> None:
    job = scheduler._parse_job_line("05:30:00 UTC x3 echo run", 1)
    assert job is not None
    first = scheduler._first_valid_occurrence(
        0, job, datetime(2026, 6, 21, 0, 0, tzinfo=timezone.utc)
    )
    assert first.repeat_label == "[1/3]"
    second = scheduler._next_occurrence(first)
    assert second.repeat_label == "[2/3]"
    third = scheduler._next_occurrence(second)
    assert third.repeat_label == "[3/3]"


def test_sleep_target_text_shows_filename_when_command_contains_time_placeholder() -> None:
    job = scheduler._parse_job_line(
        "05:30:00 UTC x3 echo -o screenshot-%t.png", 1
    )
    assert job is not None
    target = datetime(2026, 6, 21, 5, 30, tzinfo=timezone.utc)
    target_text, label = scheduler._sleep_target_text(job.command, target)
    assert target_text == "screenshot-20260621T053000Z.png"
    assert label == "for"


def test_sleep_target_text_shows_utc_when_command_has_no_time_placeholder() -> None:
    job = scheduler._parse_job_line("05:30:00 UTC echo run", 1)
    assert job is not None
    target = datetime(2026, 6, 21, 5, 30, tzinfo=timezone.utc)
    target_text, label = scheduler._sleep_target_text(job.command, target)
    assert target_text == "20260621T053000Z"
    assert label == "until"


def test_sleep_until_includes_task_line_and_repeat_label() -> None:
    job = scheduler._parse_job_line(
        "05:30:00 UTC x3 echo -o screenshot-%t.png", 8
    )
    assert job is not None
    occurrence = scheduler._first_valid_occurrence(
        0, job, datetime(2026, 6, 21, 0, 0, tzinfo=timezone.utc)
    )
    target = datetime.now(timezone.utc) + timedelta(minutes=5)
    stop_event = scheduler.Event()
    stream = _TtyBuffer()
    renderer = scheduler.WaitLineRenderer(stream=stream)
    thread = threading.Thread(
        target=scheduler._sleep_until,
        args=(target, stop_event, occurrence, renderer),
    )
    thread.start()
    time.sleep(0.2)
    stop_event.set()
    thread.join(timeout=1.0)
    output = stream.getvalue()
    assert "Waiting" in output
    assert "task 8" in output
    assert " [1/3]" in output
    assert "screenshot-" in output
    assert "[INFO] zstarview-export-image-schedule-runner" not in output


def test_sleep_until_has_no_repeat_label_for_single_run() -> None:
    job = scheduler._parse_job_line("05:30:00 UTC echo run", 1)
    assert job is not None
    occurrence = scheduler._first_valid_occurrence(
        0, job, datetime(2026, 6, 21, 0, 0, tzinfo=timezone.utc)
    )
    target = datetime.now(timezone.utc) + timedelta(minutes=5)
    stop_event = scheduler.Event()
    stream = _TtyBuffer()
    renderer = scheduler.WaitLineRenderer(stream=stream)
    thread = threading.Thread(
        target=scheduler._sleep_until,
        args=(target, stop_event, occurrence, renderer),
    )
    thread.start()
    time.sleep(0.2)
    stop_event.set()
    thread.join(timeout=1.0)
    output = stream.getvalue()
    assert "Waiting" in output
    assert "task 1" in output
    assert " [1/1]" not in output
    assert "until 20" in output


def test_sleep_until_drops_seconds_when_wait_is_longer_than_five_minutes() -> None:
    job = scheduler._parse_job_line(
        "05:30:00 UTC x3 echo -o screenshot-%t.png", 8
    )
    assert job is not None
    occurrence = scheduler._first_valid_occurrence(
        0, job, datetime(2026, 6, 21, 0, 0, tzinfo=timezone.utc)
    )
    target = datetime.now(timezone.utc) + timedelta(minutes=6, seconds=13)
    stop_event = scheduler.Event()
    stream = _TtyBuffer()
    renderer = scheduler.WaitLineRenderer(stream=stream)
    thread = threading.Thread(
        target=scheduler._sleep_until,
        args=(target, stop_event, occurrence, renderer),
    )
    thread.start()
    time.sleep(0.2)
    stop_event.set()
    thread.join(timeout=1.0)
    output = stream.getvalue()
    assert "Waiting 6m" in output
    assert "13s" not in output
    assert "task 8" in output


def test_wait_line_renderer_overwrites_tty_output() -> None:
    stream = _TtyBuffer()
    renderer = scheduler.WaitLineRenderer(stream=stream)

    renderer.write("Waiting 46m for task 12")
    renderer.write("Waiting 45m for task 12")
    renderer.separate_block()
    renderer.write("Waiting 44m for task 12")

    assert (
        stream.getvalue()
        == "\x1b[34mWaiting 46m for task 12\x1b[0m"
        "\r\x1b[2K\x1b[34mWaiting 45m for task 12\x1b[0m\n"
        "\x1b[34mWaiting 44m for task 12\x1b[0m"
    )


@pytest.mark.parametrize(
    ("remaining_seconds", "expected"),
    [
        (601.0, 60.0),
        (301.0, 60.0),
        (300.0, 10.0),
        (31.0, 10.0),
        (30.0, 10.0),
        (29.9, 1.0),
        (1.0, 1.0),
    ],
)
def test_sleep_poll_interval_seconds_switches_by_remaining_time(
    remaining_seconds: float, expected: float
) -> None:
    assert scheduler._sleep_poll_interval_seconds(remaining_seconds) == expected


@pytest.mark.parametrize(
    ("remaining_seconds", "expected"),
    [
        (720.0, 60.0),
        (360.0, 60.0),
        (330.0, 30.0),
        (301.0, 1.0),
        (300.0, 10.0),
        (35.0, 5.0),
        (31.0, 1.0),
        (30.0, 1.0),
        (29.0, 1.0),
    ],
)
def test_sleep_timeout_seconds_hits_boundary_values(
    remaining_seconds: float, expected: float
) -> None:
    assert scheduler._sleep_timeout_seconds(remaining_seconds) == expected
