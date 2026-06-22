#!/usr/bin/env python3
"""Run zstarview-export-image jobs from a daily schedule file.

This script is independent from the main zstarview application. It reads a
configuration file and keeps running until interrupted with Ctrl+C.

Config format
- One job per line.
- Blank lines and lines starting with '#' are ignored.
- Each job has the form:

    HH:MM:SS TZ [xN] COMMAND...

  where:
  - HH:MM:SS is required.
  - TZ accepts the same timezone forms as zstarview, plus bare offsets such as
    +08 or -07:30.
  - xN is optional. It means run the command N times at 5-minute intervals.
  - COMMAND is parsed with shell-style quoting and executed directly.

Execution rules
- Jobs run every day at the specified local time.
- If a local time does not exist because of DST, skip that occurrence and emit
  a warning.
- If a local time occurs twice because of DST, execute only the first
  occurrence.
- If multiple runs would overlap, do not run them in parallel. Wait for the
  current command to finish before starting the next queued run.
- The %t placeholder in COMMAND arguments expands to the actual UTC start time
  of that specific command invocation.
- The UTC timestamp format is YYYYMMDDTHHMMSSZ.
- The script logs scheduling, waiting, start, finish, exit codes, and failures.
- The script keeps running until interrupted with Ctrl+C.

Example:
    05:30:00 +08 x3 zstarview-export-image "@1.28815;103.85905" -A5 -Z135 -o screenshot-marinabay-%t.png
    04:00:00 AEST x3 zstarview-export-image -p "Circular Quay" -A5 -Z90 -o screenshot-sydney-%t.png
"""

from __future__ import annotations

import argparse
import heapq
import logging
import os
import re
import shlex
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import date, datetime, time as time_of_day, timedelta, timezone, tzinfo
from pathlib import Path
from threading import Event
from typing import Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.utils.timezone_parser import parse_tz_string  # noqa: E402

LOGGER = logging.getLogger("zstarview-export-image-schedule-runner")
REPEAT_INTERVAL = timedelta(minutes=6)
MAX_SEQUENCE_LOOKAHEAD_DAYS = 370
UTC_TIMESTAMP_FORMAT = "%Y%m%dT%H%M%SZ"
TIME_RE = re.compile(r"^(?P<h>\d{2}):(?P<m>\d{2}):(?P<s>\d{2})$")
BARE_OFFSET_RE = re.compile(r"^[+-](?P<h>\d{1,2})(?::?(?P<m>\d{2}))?$")
REPEAT_RE = re.compile(r"^x(?P<count>\d+)$", re.IGNORECASE)


class ConfigError(ValueError):
    """Raised when the schedule file contains an invalid line."""


@dataclass(frozen=True)
class JobSpec:
    line_no: int
    raw_line: str
    hour: int
    minute: int
    second: int
    timezone_text: str
    timezone_info: tzinfo
    repeat_count: int
    command: tuple[str, ...]

    @property
    def wall_time(self) -> time_of_day:
        return time_of_day(self.hour, self.minute, self.second)


@dataclass(frozen=True)
class ScheduledOccurrence:
    job_index: int
    job: JobSpec
    scheduled_day: date
    repeat_index: int
    scheduled_local: datetime
    scheduled_utc: datetime


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run zstarview-export-image jobs from a daily schedule file."
    )
    parser.add_argument(
        "config_file",
        type=Path,
        help="Path to the schedule file.",
    )
    return parser


def _parse_wall_time(value: str) -> tuple[int, int, int]:
    match = TIME_RE.fullmatch(value.strip())
    if not match:
        raise ConfigError(
            f"Invalid time {value!r}: use zero-padded HH:MM:SS."
        )
    hour = int(match.group("h"))
    minute = int(match.group("m"))
    second = int(match.group("s"))
    if not (0 <= hour <= 23):
        raise ConfigError(f"Hour out of range in {value!r}: {hour}")
    if not (0 <= minute <= 59):
        raise ConfigError(f"Minute out of range in {value!r}: {minute}")
    if not (0 <= second <= 59):
        raise ConfigError(f"Second out of range in {value!r}: {second}")
    return hour, minute, second


def _parse_timezone(value: str):
    text = value.strip()
    match = BARE_OFFSET_RE.fullmatch(text)
    if match:
        hours = int(match.group("h"))
        minutes = int(match.group("m")) if match.group("m") is not None else 0
        if not (0 <= hours <= 14):
            raise ConfigError(f"Timezone offset hour out of range in {value!r}: {hours}")
        if not (0 <= minutes < 60):
            raise ConfigError(f"Timezone offset minute out of range in {value!r}: {minutes}")
        sign = -1 if text.startswith("-") else 1
        return timezone(timedelta(hours=sign * hours, minutes=sign * minutes))
    try:
        return parse_tz_string(text)
    except Exception as exc:  # pragma: no cover - defensive, parser already validates
        raise ConfigError(f"Invalid timezone {value!r}: {exc}") from exc


def _parse_job_line(line: str, line_no: int) -> JobSpec | None:
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return None

    try:
        parts = shlex.split(line, comments=True, posix=True)
    except ValueError as exc:
        raise ConfigError(f"Line {line_no}: failed to parse shell-style quoting: {exc}") from exc

    if len(parts) < 3:
        raise ConfigError(
            f"Line {line_no}: expected 'HH:MM:SS TZ [xN] COMMAND...'"
        )

    hour, minute, second = _parse_wall_time(parts[0])
    timezone_text = parts[1]
    repeat_count = 1
    command_start = 2

    repeat_match = REPEAT_RE.fullmatch(parts[2])
    if repeat_match:
        repeat_count = int(repeat_match.group("count"))
        if repeat_count <= 0:
            raise ConfigError(f"Line {line_no}: xN must be greater than zero.")
        command_start = 3

    command = tuple(parts[command_start:])
    if not command:
        raise ConfigError(f"Line {line_no}: missing command after schedule fields.")

    timezone_info = _parse_timezone(timezone_text)
    return JobSpec(
        line_no=line_no,
        raw_line=line.rstrip("\n"),
        hour=hour,
        minute=minute,
        second=second,
        timezone_text=timezone_text,
        timezone_info=timezone_info,
        repeat_count=repeat_count,
        command=command,
    )


def _load_jobs(config_path: Path) -> list[JobSpec]:
    try:
        lines = config_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ConfigError(f"Failed to read {config_path}: {exc}") from exc

    jobs: list[JobSpec] = []
    for line_no, line in enumerate(lines, start=1):
        job = _parse_job_line(line, line_no)
        if job is not None:
            jobs.append(job)
    if not jobs:
        raise ConfigError("No jobs found in the schedule file.")
    return jobs


def _format_aware(dt: datetime) -> str:
    return dt.strftime("%Y-%m-%d %H:%M:%S %Z%z")


def _format_utc(dt: datetime) -> str:
    return dt.astimezone(timezone.utc).strftime(UTC_TIMESTAMP_FORMAT)


def _localize_wall_time(naive_local: datetime, tzinfo) -> datetime | None:
    candidates: list[datetime] = []
    seen: set[datetime] = set()
    for fold in (0, 1):
        aware = naive_local.replace(tzinfo=tzinfo, fold=fold)
        round_trip = aware.astimezone(timezone.utc).astimezone(tzinfo)
        if round_trip.replace(tzinfo=None) != naive_local:
            continue
        utc_instant = aware.astimezone(timezone.utc)
        if utc_instant in seen:
            continue
        seen.add(utc_instant)
        candidates.append(aware)
    if not candidates:
        return None
    candidates.sort(key=lambda value: value.astimezone(timezone.utc))
    return candidates[0]


def _wall_time_for(job: JobSpec, scheduled_day: date, repeat_index: int) -> datetime:
    base = datetime.combine(scheduled_day, job.wall_time)
    return base + repeat_index * REPEAT_INTERVAL


def _advance_sequence(day: date, repeat_index: int, job: JobSpec) -> tuple[date, int]:
    if repeat_index + 1 < job.repeat_count:
        return day, repeat_index + 1
    return day + timedelta(days=1), 0


def _first_valid_occurrence(
    job_index: int,
    job: JobSpec,
    now_utc: datetime,
) -> ScheduledOccurrence:
    start_day = now_utc.astimezone(job.timezone_info).date()
    day = start_day
    repeat_index = 0
    for _ in range(MAX_SEQUENCE_LOOKAHEAD_DAYS * max(1, job.repeat_count)):
        naive_local = _wall_time_for(job, day, repeat_index)
        localized = _localize_wall_time(naive_local, job.timezone_info)
        if localized is None:
            LOGGER.warning(
                "Line %d: skipping nonexistent local time %s in timezone %s",
                job.line_no,
                naive_local.strftime("%Y-%m-%d %H:%M:%S"),
                job.timezone_text,
            )
        else:
            localized_utc = localized.astimezone(timezone.utc)
            if localized_utc >= now_utc:
                return ScheduledOccurrence(
                    job_index=job_index,
                    job=job,
                    scheduled_day=day,
                    repeat_index=repeat_index,
                    scheduled_local=localized,
                    scheduled_utc=localized_utc,
                )
        day, repeat_index = _advance_sequence(day, repeat_index, job)
    raise ConfigError(
        f"Line {job.line_no}: could not find a future occurrence within the lookahead window."
    )


def _next_occurrence(occurrence: ScheduledOccurrence) -> ScheduledOccurrence:
    job = occurrence.job
    day, repeat_index = _advance_sequence(
        occurrence.scheduled_day,
        occurrence.repeat_index,
        job,
    )
    for _ in range(MAX_SEQUENCE_LOOKAHEAD_DAYS * max(1, job.repeat_count)):
        naive_local = _wall_time_for(job, day, repeat_index)
        localized = _localize_wall_time(naive_local, job.timezone_info)
        if localized is None:
            LOGGER.warning(
                "Line %d: skipping nonexistent local time %s in timezone %s",
                job.line_no,
                naive_local.strftime("%Y-%m-%d %H:%M:%S"),
                job.timezone_text,
            )
        else:
            localized_utc = localized.astimezone(timezone.utc)
            return ScheduledOccurrence(
                job_index=occurrence.job_index,
                job=job,
                scheduled_day=day,
                repeat_index=repeat_index,
                scheduled_local=localized,
                scheduled_utc=localized_utc,
            )
        day, repeat_index = _advance_sequence(day, repeat_index, job)
    raise ConfigError(
        f"Line {job.line_no}: could not find the next occurrence within the lookahead window."
    )


def _expand_placeholders(command: Sequence[str], *, start_time_utc: datetime) -> tuple[str, ...]:
    stamp = _format_utc(start_time_utc)
    return tuple(part.replace("%t", stamp) for part in command)


def _terminate_process(proc: subprocess.Popen[str]) -> None:
    if proc.poll() is not None:
        return
    if os.name == "posix":
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
        deadline = time.monotonic() + 10.0
        while proc.poll() is None and time.monotonic() < deadline:
            time.sleep(0.2)
        if proc.poll() is None:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                return
    else:  # pragma: no cover - Windows fallback
        proc.terminate()


def _run_command(command: Sequence[str], stop_event: Event) -> int:
    start_utc = datetime.now(timezone.utc)
    expanded = _expand_placeholders(command, start_time_utc=start_utc)
    LOGGER.info("Starting command at %s: %s", _format_utc(start_utc), shlex.join(expanded))
    popen_kwargs: dict[str, object] = {}
    if os.name == "posix":
        popen_kwargs["start_new_session"] = True
    proc = subprocess.Popen(expanded, **popen_kwargs)  # type: ignore[arg-type]
    try:
        while True:
            rc = proc.poll()
            if rc is not None:
                elapsed = datetime.now(timezone.utc) - start_utc
                LOGGER.info(
                    "Command finished at %s after %.1fs with exit code %d",
                    _format_utc(datetime.now(timezone.utc)),
                    elapsed.total_seconds(),
                    rc,
                )
                return rc
            if stop_event.is_set():
                LOGGER.info("Stop requested; terminating running command.")
                _terminate_process(proc)
                rc = proc.wait()
                elapsed = datetime.now(timezone.utc) - start_utc
                LOGGER.info(
                    "Command stopped at %s after %.1fs with exit code %d",
                    _format_utc(datetime.now(timezone.utc)),
                    elapsed.total_seconds(),
                    rc,
                )
                return rc
            stop_event.wait(timeout=0.5)
    except KeyboardInterrupt:
        stop_event.set()
        LOGGER.info("KeyboardInterrupt received; terminating running command.")
        _terminate_process(proc)
        rc = proc.wait()
        elapsed = datetime.now(timezone.utc) - start_utc
        LOGGER.info(
            "Command interrupted at %s after %.1fs with exit code %d",
            _format_utc(datetime.now(timezone.utc)),
            elapsed.total_seconds(),
            rc,
        )
        raise


def _install_signal_handlers(stop_event: Event):
    previous = {}

    def _handler(signum, frame):  # noqa: ARG001
        stop_event.set()

    for sig in (signal.SIGINT, signal.SIGTERM):
        previous[sig] = signal.getsignal(sig)
        signal.signal(sig, _handler)
    return previous


def _restore_signal_handlers(previous) -> None:
    for sig, handler in previous.items():
        signal.signal(sig, handler)


def _sleep_until(target_utc: datetime, stop_event: Event) -> None:
    while not stop_event.is_set():
        now_utc = datetime.now(timezone.utc)
        remaining = (target_utc - now_utc).total_seconds()
        if remaining <= 0:
            return
        minutes, seconds = divmod(int(remaining), 60)
        LOGGER.info(
            "Waiting %dm %ds until %s",
            minutes,
            seconds,
            _format_utc(target_utc),
        )
        stop_event.wait(timeout=min(remaining, 60.0))


def run_scheduler(config_path: Path) -> int:
    jobs = _load_jobs(config_path)
    now_utc = datetime.now(timezone.utc)
    stop_event = Event()
    previous_signals = _install_signal_handlers(stop_event)

    try:
        LOGGER.info("Loaded %d jobs from %s", len(jobs), config_path)
        queue: list[tuple[datetime, int, ScheduledOccurrence]] = []
        sequence = 0
        for job_index, job in enumerate(jobs):
            occurrence = _first_valid_occurrence(job_index, job, now_utc)
            LOGGER.info(
                "Line %d scheduled next at %s (%s UTC) -> %s",
                job.line_no,
                _format_aware(occurrence.scheduled_local),
                _format_utc(occurrence.scheduled_utc),
                shlex.join(job.command),
            )
            heapq.heappush(queue, (occurrence.scheduled_utc, sequence, occurrence))
            sequence += 1

        while queue and not stop_event.is_set():
            scheduled_utc, _, occurrence = heapq.heappop(queue)
            now_utc = datetime.now(timezone.utc)
            if scheduled_utc > now_utc:
                _sleep_until(scheduled_utc, stop_event)
                if stop_event.is_set():
                    break
            LOGGER.info(
                "Running line %d at %s (scheduled %s)",
                occurrence.job.line_no,
                _format_utc(datetime.now(timezone.utc)),
                _format_utc(occurrence.scheduled_utc),
            )
            exit_code = _run_command(occurrence.job.command, stop_event)
            LOGGER.info("Line %d completed with exit code %d", occurrence.job.line_no, exit_code)
            next_occurrence = _next_occurrence(occurrence)
            heapq.heappush(queue, (next_occurrence.scheduled_utc, sequence, next_occurrence))
            sequence += 1

        if stop_event.is_set():
            LOGGER.info("Scheduler stopped by request.")
        return 130 if stop_event.is_set() else 0
    finally:
        _restore_signal_handlers(previous_signals)


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )
    try:
        return run_scheduler(args.config_file)
    except ConfigError as exc:
        LOGGER.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        LOGGER.info("Interrupted.")
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
