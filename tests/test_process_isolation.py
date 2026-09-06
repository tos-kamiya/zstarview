from __future__ import annotations

import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from zstarview.processes import JobRequest, ProcessJobSupervisor
from zstarview.tropical_cyclones.cache import TropicalCycloneCacheEntry, save_tropical_cyclone_cache
from zstarview.tropical_cyclones.models import TropicalCycloneSnapshotCollection


WORKER = """
import json, pathlib, sys, time
args = iter(sys.argv[1:])
options = dict(zip(args, args))
request_path = pathlib.Path(options['--request'])
result_path = pathlib.Path(options['--result'])
request = json.loads(request_path.read_text())
mode = request['payload'].get('mode')
if mode == 'crash':
    raise SystemExit(9)
if mode == 'sleep':
    time.sleep(10)
if mode == 'corrupt':
    result_path.write_text('{not json')
    raise SystemExit(0)
result = dict(request)
result.update({'status': 'ok', 'artifacts': []})
result_path.write_text(json.dumps(result))
"""


def _request(session_id: str, request_id: int, mode: str = "ok") -> JobRequest:
    return JobRequest(
        session_id=session_id,
        worker_epoch=1,
        request_id=request_id,
        job_kind="test",
        layer_generation=2,
        view_generation=request_id,
        input_revision="input-1",
        payload={"mode": mode},
    )


def _command() -> tuple[str, ...]:
    return (sys.executable, "-c", WORKER)


def _wait(supervisor: ProcessJobSupervisor):
    deadline = time.monotonic() + 3.0
    while time.monotonic() < deadline:
        events = supervisor.poll()
        if events:
            return events[0]
        time.sleep(0.01)
    raise AssertionError("worker did not finish")


def test_supervisor_accepts_result_and_replaces_pending_request(tmp_path: Path) -> None:
    supervisor = ProcessJobSupervisor(tmp_path, session_id="session")
    supervisor.submit(_request("session", 1, "sleep"), _command(), timeout_s=2.0)
    supervisor.submit(_request("session", 2), _command(), timeout_s=2.0)

    first = _wait(supervisor)
    assert first.failure is not None
    assert first.failure.kind == "timeout"
    second = _wait(supervisor)
    assert second.result is not None
    assert second.result.request_id == 2
    assert supervisor.close() == ()


def test_supervisor_reports_worker_crash(tmp_path: Path) -> None:
    supervisor = ProcessJobSupervisor(tmp_path, session_id="session")
    supervisor.submit(_request("session", 1, "crash"), _command(), timeout_s=2.0)
    event = _wait(supervisor)
    assert event.failure is not None
    assert event.failure.kind == "worker_exit"


def test_supervisor_reports_corrupt_result(tmp_path: Path) -> None:
    supervisor = ProcessJobSupervisor(tmp_path, session_id="session")
    supervisor.submit(_request("session", 1, "corrupt"), _command(), timeout_s=2.0)
    event = _wait(supervisor)
    assert event.failure is not None
    assert event.failure.kind == "invalid_result"


def test_tropical_cyclone_worker_uses_json_artifact(tmp_path: Path) -> None:
    cache_root = tmp_path / "cache"
    collection = TropicalCycloneSnapshotCollection(
        snapshots=(),
        source_url="https://example.invalid",
        refreshed_at_utc=datetime.now(timezone.utc),
    )
    save_tropical_cyclone_cache(
        TropicalCycloneCacheEntry(
            snapshot_collection=collection,
            cached_at_utc=datetime.now(timezone.utc),
        ),
        cache_root=cache_root,
    )
    session_id = "cyclone-session"
    request = JobRequest(
        session_id=session_id,
        worker_epoch=1,
        request_id=1,
        job_kind="tropical-cyclone-fetch",
        layer_generation=1,
        view_generation=0,
        input_revision="4",
        payload={"cache_root": str(cache_root)},
    )
    supervisor = ProcessJobSupervisor(tmp_path / "jobs", session_id=session_id)
    supervisor.submit(
        request,
        (sys.executable, "-m", "zstarview.tropical_cyclones.worker"),
        timeout_s=3.0,
    )
    event = _wait(supervisor)
    assert event.failure is None
    assert event.result is not None
    assert event.result.status == "ok"
    assert event.result.artifacts[0].schema == "zstarview.tropical-cyclone-payload.v1"
