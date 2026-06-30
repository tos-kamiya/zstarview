# -*- coding: utf-8 -*-
"""Structured diagnostics for cloud source acquisition."""

from __future__ import annotations

import datetime as dt
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Protocol

DiagnosticStatus = Literal["start", "ok", "failed", "info"]


@dataclass(frozen=True, slots=True)
class DiagnosticEvent:
    """A single cloud-source diagnostic event."""

    step: str
    status: DiagnosticStatus
    message: str
    fields: dict[str, Any] = field(default_factory=dict)
    timestamp_utc: dt.datetime = field(default_factory=lambda: dt.datetime.now(dt.timezone.utc))

    def to_json_dict(self) -> dict[str, Any]:
        return {
            "timestamp_utc": self.timestamp_utc.astimezone(dt.timezone.utc).isoformat().replace("+00:00", "Z"),
            "step": self.step,
            "status": self.status,
            "message": self.message,
            "fields": _json_safe(self.fields),
        }


class DiagnosticSink(Protocol):
    """Receiver for cloud-source diagnostic events."""

    def emit(self, event: DiagnosticEvent) -> None:
        """Record one diagnostic event."""


def emit_diagnostic(
    sink: DiagnosticSink | None,
    step: str,
    status: DiagnosticStatus,
    message: str,
    **fields: Any,
) -> None:
    """Emit a diagnostic event if a sink was provided."""
    if sink is None:
        return
    sink.emit(DiagnosticEvent(step=step, status=status, message=message, fields=fields))


class FileDiagnosticSink:
    """Write cloud-source diagnostics to JSONL and optional text files."""

    def __init__(self, jsonl_path: Path, text_path: Path | None = None) -> None:
        self.jsonl_path = Path(jsonl_path)
        self.text_path = None if text_path is None else Path(text_path)
        self.jsonl_path.parent.mkdir(parents=True, exist_ok=True)
        if self.text_path is not None:
            self.text_path.parent.mkdir(parents=True, exist_ok=True)

    def emit(self, event: DiagnosticEvent) -> None:
        payload = event.to_json_dict()
        with self.jsonl_path.open("a", encoding="utf-8") as fp:
            fp.write(json.dumps(payload, ensure_ascii=True, sort_keys=True))
            fp.write("\n")
        if self.text_path is not None:
            with self.text_path.open("a", encoding="utf-8") as fp:
                fp.write(_format_text_event(event))
                fp.write("\n")


def _format_text_event(event: DiagnosticEvent) -> str:
    timestamp = event.timestamp_utc.astimezone(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    details = " ".join(f"{key}={_format_text_value(value)}" for key, value in sorted(event.fields.items()))
    if details:
        return f"{timestamp} [{event.status.upper()}] {event.step}: {event.message} ({details})"
    return f"{timestamp} [{event.status.upper()}] {event.step}: {event.message}"


def _format_text_value(value: Any) -> str:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dt.datetime):
        return value.astimezone(dt.timezone.utc).isoformat().replace("+00:00", "Z")
    return str(value)


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dt.datetime):
        return value.astimezone(dt.timezone.utc).isoformat().replace("+00:00", "Z")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)
