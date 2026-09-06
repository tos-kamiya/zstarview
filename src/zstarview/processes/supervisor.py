"""Non-blocking supervisor for one running job and one latest pending job."""

from __future__ import annotations

import os
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


@dataclass(slots=True)
class _RunningJob:
    request: JobRequest
    directory: Path
    result_path: Path
    process: subprocess.Popen[str]
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
    ) -> None:
        self.work_root = Path(work_root).resolve()
        self.session_id = session_id or uuid.uuid4().hex
        if max_artifact_bytes < 0:
            raise ValueError("max_artifact_bytes must be non-negative")
        self.max_artifact_bytes = max_artifact_bytes
        self._running: _RunningJob | None = None
        self._pending: tuple[JobRequest, tuple[str, ...], float, float] | None = None
        self._closed = False

    @property
    def running_request(self) -> JobRequest | None:
        return None if self._running is None else self._running.request

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
        return (self._finish(current, current.process.returncode),)

    def _start_pending(self) -> None:
        if self._running is not None or self._pending is None or self._closed:
            return
        request, command, timeout_s, terminate_grace_s = self._pending
        self._pending = None
        directory = self.work_root / self.session_id / f"job-{request.request_id}"
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
            return SupervisorEvent(current.request, failure=failure)
        if exit_code != 0:
            failure = ProcessFailure(
                "worker_exit",
                f"worker exited with status {exit_code}: {stderr[-1000:]}",
                exit_code=exit_code,
            )
            return SupervisorEvent(current.request, failure=failure)
        try:
            result = JobResult.read(current.result_path, expected_request=current.request)
        except ProcessFailure as failure:
            return SupervisorEvent(current.request, failure=failure)
        for artifact in result.artifacts:
            if not artifact.relative_path or artifact.size_bytes > self.max_artifact_bytes:
                return SupervisorEvent(
                    current.request,
                    failure=ProcessFailure("invalid_result", "artifact exceeds supervisor limits"),
                )
            artifact_path = (current.directory / artifact.relative_path).resolve()
            try:
                artifact_path.relative_to(current.directory.resolve())
            except ValueError:
                return SupervisorEvent(
                    current.request,
                    failure=ProcessFailure("invalid_result", "artifact escapes job directory"),
                )
            if not artifact_path.is_file() or artifact_path.stat().st_size != artifact.size_bytes:
                return SupervisorEvent(
                    current.request,
                    failure=ProcessFailure("invalid_result", "artifact does not match manifest"),
                )
        return SupervisorEvent(current.request, result=result)
