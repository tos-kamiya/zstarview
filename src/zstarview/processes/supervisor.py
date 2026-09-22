"""Non-blocking supervisor for one running job and one latest pending job."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from .protocol import JobRequest, JobResult, ProcessFailure


@dataclass(frozen=True, slots=True)
class SupervisorEvent:
    request: JobRequest
    result: JobResult | None = None
    failure: ProcessFailure | None = None
    directory: Path | None = None


@dataclass(slots=True)
class _RunningJob:
    request: JobRequest
    directory: Path
    result_path: Path
    process: subprocess.Popen[str]
    command: tuple[str, ...]
    timeout_s: float
    deadline: float
    terminate_grace_s: float
    kill_at: float | None = None
    timed_out: bool = False


class ProcessJobSupervisor:
    """Supervise isolated command jobs without a supervisor thread.

    The GUI can call :meth:`poll` from its event loop. A newer submission replaces
    the pending request, while the current process is allowed to finish or timeout.
    """

    def __init__(
        self,
        work_root: Path,
        *,
        session_id: str | None = None,
        max_artifact_bytes: int = 512 * 1024 * 1024,
        max_retained_jobs: int = 8,
        max_retained_bytes: int = 64 * 1024 * 1024,
        max_restarts: int = 0,
        restart_window_s: float = 300.0,
        restart_backoff_s: tuple[float, ...] = (0.5, 2.0),
    ) -> None:
        self.work_root = Path(work_root).resolve()
        self.session_id = session_id or uuid.uuid4().hex
        if max_artifact_bytes < 0:
            raise ValueError("max_artifact_bytes must be non-negative")
        if max_retained_jobs < 0 or max_retained_bytes < 0 or max_restarts < 0:
            raise ValueError("retention and restart limits must be non-negative")
        if restart_window_s <= 0 or any(delay < 0 for delay in restart_backoff_s):
            raise ValueError("restart window and backoff must be non-negative")
        self.max_artifact_bytes = max_artifact_bytes
        self.max_retained_jobs = max_retained_jobs
        self.max_retained_bytes = max_retained_bytes
        self.max_restarts = max_restarts
        self.restart_window_s = restart_window_s
        self.restart_backoff_s = restart_backoff_s
        self._running: _RunningJob | None = None
        self._pending: tuple[JobRequest, tuple[str, ...], float, float] | None = None
        self._pending_due = 0.0
        self._restart_count = 0
        self._restart_started_at = 0.0
        self._retained_directories: list[Path] = []
        self._closed = False

    @property
    def running_request(self) -> JobRequest | None:
        return None if self._running is None else self._running.request

    @property
    def has_pending_work(self) -> bool:
        return self._running is not None or self._pending is not None

    def submit(
        self,
        request: JobRequest,
        command: Sequence[str],
        *,
        timeout_s: float,
        terminate_grace_s: float = 1.0,
    ) -> None:
        if self._closed:
            raise RuntimeError("process supervisor is closed")
        if timeout_s <= 0 or terminate_grace_s < 0:
            raise ValueError("timeouts must be non-negative and timeout_s must be positive")
        if request.session_id != self.session_id:
            raise ValueError("request belongs to another supervisor session")
        if not command:
            raise ValueError("worker command must not be empty")
        self._pending = (request, tuple(str(part) for part in command), timeout_s, terminate_grace_s)
        self._pending_due = 0.0
        self._start_pending()

    def poll(self, *, now: float | None = None) -> tuple[SupervisorEvent, ...]:
        current = self._running
        if current is None:
            self._start_pending()
            return ()
        timestamp = time.monotonic() if now is None else now
        exit_code = current.process.poll()
        if exit_code is None and current.kill_at is None and timestamp >= current.deadline:
            current.process.terminate()
            current.timed_out = True
            current.kill_at = timestamp + current.terminate_grace_s
            return ()
        if exit_code is None and current.kill_at is not None and timestamp >= current.kill_at:
            current.process.kill()
            current.kill_at = timestamp + 0.01
            return ()
        if exit_code is None:
            return ()
        self._running = None
        event = self._finish(current, exit_code)
        if event.failure is not None and self._schedule_restart(current, event.failure, timestamp):
            self._discard_directory(current.directory)
            self._start_pending(now=timestamp)
            return ()
        self._start_pending()
        return (event,)

    def close(self) -> tuple[SupervisorEvent, ...]:
        self._closed = True
        self._pending = None
        if self._running is None:
            return ()
        self._running.process.terminate()
        try:
            self._running.process.wait(timeout=1.0)
        except subprocess.TimeoutExpired:
            self._running.process.kill()
            self._running.process.wait()
        current = self._running
        self._running = None
        event = self._finish(current, current.process.returncode)
        self._discard_directory(current.directory)
        return (event,)

    def release(self, event: SupervisorEvent) -> None:
        """Release a completed job directory after its artifacts are consumed."""
        directory = event.directory
        if directory is None:
            return
        if directory in self._retained_directories:
            self._retained_directories.remove(directory)
        self._discard_directory(directory)

    def _start_pending(self, *, now: float | None = None) -> None:
        if self._running is not None or self._pending is None or self._closed:
            return
        timestamp = time.monotonic() if now is None else now
        if timestamp < self._pending_due:
            return
        request, command, timeout_s, terminate_grace_s = self._pending
        self._pending = None
        if request.request_id != getattr(self, "_active_request_id", None):
            self._restart_count = 0
            self._restart_started_at = timestamp
            self._active_request_id = request.request_id
        self._prune_retained()
        directory = self.work_root / self.session_id / f"job-{request.request_id}-{self._restart_count}"
        directory.mkdir(parents=True, exist_ok=False)
        request_path = directory / "request.json"
        result_path = directory / "result.json"
        request.write(request_path)
        process = subprocess.Popen(
            (*command, "--request", os.fspath(request_path), "--result", os.fspath(result_path)),
            cwd=directory,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
        self._running = _RunningJob(
            request=request,
            directory=directory,
            result_path=result_path,
            process=process,
            command=command,
            timeout_s=timeout_s,
            deadline=time.monotonic() + timeout_s,
            terminate_grace_s=terminate_grace_s,
            kill_at=None,
        )

    def _finish(self, current: _RunningJob, exit_code: int | None) -> SupervisorEvent:
        stderr = current.process.stderr.read() if current.process.stderr is not None else ""
        if current.timed_out:
            failure = ProcessFailure(
                "timeout",
                "worker exceeded its deadline",
                exit_code=exit_code,
            )
            event = SupervisorEvent(current.request, failure=failure, directory=current.directory)
            self._retain_directory(current.directory)
            return event
        if exit_code != 0:
            failure = ProcessFailure(
                "worker_exit",
                f"worker exited with status {exit_code}: {stderr[-1000:]}",
                exit_code=exit_code,
            )
            event = SupervisorEvent(current.request, failure=failure, directory=current.directory)
            self._retain_directory(current.directory)
            return event
        try:
            result = JobResult.read(current.result_path, expected_request=current.request)
        except ProcessFailure as failure:
            event = SupervisorEvent(current.request, failure=failure, directory=current.directory)
            self._retain_directory(current.directory)
            return event
        for artifact in result.artifacts:
            if not artifact.relative_path or artifact.size_bytes > self.max_artifact_bytes:
                event = SupervisorEvent(current.request, failure=ProcessFailure("invalid_result", "artifact exceeds supervisor limits"), directory=current.directory)
                self._retain_directory(current.directory)
                return event
            artifact_path = (current.directory / artifact.relative_path).resolve()
            try:
                artifact_path.relative_to(current.directory.resolve())
            except ValueError:
                event = SupervisorEvent(current.request, failure=ProcessFailure("invalid_result", "artifact escapes job directory"), directory=current.directory)
                self._retain_directory(current.directory)
                return event
            if not artifact_path.is_file() or artifact_path.stat().st_size != artifact.size_bytes:
                event = SupervisorEvent(current.request, failure=ProcessFailure("invalid_result", "artifact does not match manifest"), directory=current.directory)
                self._retain_directory(current.directory)
                return event
        event = SupervisorEvent(current.request, result=result, directory=current.directory)
        self._retain_directory(current.directory)
        return event

    def _retain_directory(self, directory: Path) -> None:
        self._retained_directories.append(directory)
        self._prune_retained(protected=directory)

    def _schedule_restart(
        self, current: _RunningJob, failure: ProcessFailure, now: float
    ) -> bool:
        if failure.kind not in {"worker_exit", "timeout"} or self.max_restarts == 0:
            return False
        if now - self._restart_started_at > self.restart_window_s:
            self._restart_count = 0
            self._restart_started_at = now
        if self._restart_count >= self.max_restarts:
            return False
        delay_index = min(self._restart_count, len(self.restart_backoff_s) - 1)
        delay = self.restart_backoff_s[delay_index] if self.restart_backoff_s else 0.0
        self._restart_count += 1
        self._pending = (current.request, current.command, current.timeout_s, current.terminate_grace_s)
        self._pending_due = now + delay
        return True

    def _discard_directory(self, directory: Path) -> None:
        if directory in self._retained_directories:
            self._retained_directories.remove(directory)
        session_root = (self.work_root / self.session_id).resolve()
        resolved = directory.resolve()
        try:
            resolved.relative_to(session_root)
        except ValueError:
            return
        shutil.rmtree(resolved, ignore_errors=True)

    def _prune_retained(self, *, protected: Path | None = None) -> None:
        while len(self._retained_directories) > self.max_retained_jobs:
            candidate = next(
                (path for path in self._retained_directories if path != protected),
                None,
            )
            if candidate is None:
                break
            self._retained_directories.remove(candidate)
            self._discard_directory(candidate)
        total = 0
        for directory in self._retained_directories:
            if directory.exists():
                total += sum(path.stat().st_size for path in directory.rglob("*") if path.is_file())
        while total > self.max_retained_bytes and self._retained_directories:
            directory_to_remove = next(
                (path for path in self._retained_directories if path != protected),
                None,
            )
            if directory_to_remove is None:
                break
            self._retained_directories.remove(directory_to_remove)
            size = sum(
                path.stat().st_size
                for path in directory_to_remove.rglob("*")
                if path.is_file()
            ) if directory_to_remove.exists() else 0
            self._discard_directory(directory_to_remove)
            total -= size
