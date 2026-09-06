"""Versioned, Qt-independent protocol for isolated jobs."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

IPC_PROTOCOL_VERSION = 1


class ProcessFailure(RuntimeError):
    """A worker did not produce an acceptable result."""

    def __init__(self, kind: str, message: str, *, exit_code: int | None = None) -> None:
        super().__init__(message)
        self.kind = kind
        self.exit_code = exit_code


def _json_object(value: dict[str, Any]) -> str:
    try:
        return json.dumps(value, ensure_ascii=True, allow_nan=False, sort_keys=True)
    except (TypeError, ValueError) as exc:
        raise ValueError("IPC payload must contain JSON-compatible finite values") from exc


@dataclass(frozen=True, slots=True)
class JobRequest:
    """Input envelope written to a job directory before worker startup."""

    session_id: str
    worker_epoch: int
    request_id: int
    job_kind: str
    layer_generation: int
    view_generation: int
    input_revision: str
    payload: dict[str, Any] = field(default_factory=dict)
    deadline_utc: str | None = None
    time_budget_s: float | None = None

    def to_dict(self) -> dict[str, Any]:
        result = {
            "protocol_version": IPC_PROTOCOL_VERSION,
            "session_id": self.session_id,
            "worker_epoch": self.worker_epoch,
            "request_id": self.request_id,
            "job_kind": self.job_kind,
            "layer_generation": self.layer_generation,
            "view_generation": self.view_generation,
            "input_revision": self.input_revision,
            "deadline_utc": self.deadline_utc,
            "time_budget_s": self.time_budget_s,
            "payload": self.payload,
        }
        _json_object(result)
        return result

    def write(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_text(_json_object(self.to_dict()), encoding="utf-8")
        temporary.replace(path)


@dataclass(frozen=True, slots=True)
class ProcessArtifact:
    """A result file owned by the worker job directory."""

    relative_path: str
    schema: str
    size_bytes: int
    shape: tuple[int, ...] | None = None
    dtype: str | None = None

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "ProcessArtifact":
        shape_value = value.get("shape")
        shape = None if shape_value is None else tuple(int(item) for item in shape_value)
        artifact = cls(
            relative_path=str(value["relative_path"]),
            schema=str(value["schema"]),
            size_bytes=int(value["size_bytes"]),
            shape=shape,
            dtype=None if value.get("dtype") is None else str(value["dtype"]),
        )
        if artifact.size_bytes < 0 or any(item < 0 for item in artifact.shape or ()):
            raise ProcessFailure("invalid_result", "artifact size or shape is invalid")
        return artifact


@dataclass(frozen=True, slots=True)
class JobResult:
    """Validated result envelope read after the worker exits or reports done."""

    session_id: str
    worker_epoch: int
    request_id: int
    job_kind: str
    layer_generation: int
    view_generation: int
    status: str
    artifacts: tuple[ProcessArtifact, ...] = ()
    data_time_utc: str | None = None
    coverage: dict[str, Any] | None = None
    processing_time_s: float | None = None
    error_kind: str | None = None
    error_message: str | None = None

    @classmethod
    def read(cls, path: Path, *, expected_request: JobRequest) -> "JobResult":
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise ProcessFailure("invalid_result", "result manifest is not valid JSON") from exc
        if not isinstance(value, dict) or value.get("protocol_version") != IPC_PROTOCOL_VERSION:
            raise ProcessFailure("invalid_result", "unsupported or missing IPC protocol version")
        for name in ("session_id", "worker_epoch", "request_id", "job_kind"):
            if value.get(name) != getattr(expected_request, name):
                raise ProcessFailure("stale_result", f"result {name} does not match request")
        for name in ("layer_generation", "view_generation"):
            if value.get(name) != getattr(expected_request, name):
                raise ProcessFailure("stale_result", f"result {name} does not match request")
        status = value.get("status")
        if status not in {"ok", "failed", "cancelled"}:
            raise ProcessFailure("invalid_result", "result status is invalid")
        processing_time = value.get("processing_time_s")
        if processing_time is not None and (
            not isinstance(processing_time, (int, float)) or not math.isfinite(processing_time)
        ):
            raise ProcessFailure("invalid_result", "processing time is invalid")
        artifacts_value = value.get("artifacts", [])
        if not isinstance(artifacts_value, list):
            raise ProcessFailure("invalid_result", "artifacts must be a list")
        if any(not isinstance(item, dict) for item in artifacts_value):
            raise ProcessFailure("invalid_result", "artifact entries must be objects")
        return cls(
            session_id=str(value["session_id"]),
            worker_epoch=int(value["worker_epoch"]),
            request_id=int(value["request_id"]),
            job_kind=str(value["job_kind"]),
            layer_generation=int(value.get("layer_generation", -1)),
            view_generation=int(value.get("view_generation", -1)),
            status=status,
            artifacts=tuple(ProcessArtifact.from_dict(item) for item in artifacts_value),
            data_time_utc=value.get("data_time_utc"),
            coverage=value.get("coverage"),
            processing_time_s=None if processing_time is None else float(processing_time),
            error_kind=value.get("error_kind"),
            error_message=value.get("error_message"),
        )
