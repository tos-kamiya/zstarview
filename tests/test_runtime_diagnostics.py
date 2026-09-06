from __future__ import annotations

import json

from zstarview.cli import diagnose_runtime
from zstarview.runtime_diagnostics import (
    collect_runtime_diagnostics,
    format_runtime_diagnostics,
)


def test_runtime_diagnostics_contains_process_and_dependency_fields(monkeypatch) -> None:
    monkeypatch.setenv("ZSTARVIEW_APP_REVISION", "test-revision")
    diagnostics = collect_runtime_diagnostics(session_id="session-1", worker_epoch=3)

    assert diagnostics["app_revision"] == "test-revision"
    assert diagnostics["session_id"] == "session-1"
    assert diagnostics["worker_epoch"] == 3
    assert diagnostics["pid"] > 0
    assert diagnostics["python_version"]
    assert diagnostics["gil_state"] in {"enabled", "disabled", "unknown"}
    assert "numpy" in diagnostics["dependencies"]


def test_runtime_diagnostics_format_is_ascii_json() -> None:
    rendered = format_runtime_diagnostics({"message": "ascii", "nested": {"value": "ok"}})

    assert json.loads(rendered) == {"message": "ascii", "nested": {"value": "ok"}}
    assert rendered.isascii()


def test_runtime_diagnostics_cli_emits_ascii_key_values(capsys, monkeypatch) -> None:
    monkeypatch.setattr(
        diagnose_runtime,
        "collect_runtime_diagnostics",
        lambda: {"cwd": "Matsue\u65e5\u672c", "dependencies": {"numpy": "1.0"}},
    )

    assert diagnose_runtime.main([]) == 0
    output = capsys.readouterr().out
    assert output.isascii()
    assert "cwd=Matsue\\u65e5\\u672c" in output
