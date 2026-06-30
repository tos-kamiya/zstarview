from __future__ import annotations

import json
import subprocess
from pathlib import Path

from zstarview.cli import diagnose_cloud_source
from zstarview.clouddisc.diagnostics import FileDiagnosticSink, emit_diagnostic


def test_file_diagnostic_sink_writes_jsonl_and_text(tmp_path: Path) -> None:
    jsonl_path = tmp_path / "diagnostics.jsonl"
    log_path = tmp_path / "diagnose.log"
    sink = FileDiagnosticSink(jsonl_path, log_path)

    emit_diagnostic(
        sink,
        "download_source",
        "ok",
        "Source ready",
        path=tmp_path / "source.nc",
    )

    payload = json.loads(jsonl_path.read_text(encoding="utf-8"))
    assert payload["step"] == "download_source"
    assert payload["status"] == "ok"
    assert payload["fields"]["path"] == str(tmp_path / "source.nc")
    assert "download_source: Source ready" in log_path.read_text(encoding="utf-8")


def test_diagnose_cloud_source_uses_isolated_default_cache(
    monkeypatch, tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "diagnosis"
    commands: list[list[str]] = []
    envs: list[dict[str, str]] = []

    def fake_run(cmd, *, cwd, check, env):  # noqa: ANN001, ANN202
        commands.append(list(cmd))
        envs.append(dict(env))
        worker_dir = output_dir / "worker"
        worker_dir.mkdir(parents=True, exist_ok=True)
        (worker_dir / "result.json").write_text(
            json.dumps(
                {
                    "request_id": 7,
                    "status": "succeeded",
                    "started_at_utc": "2026-06-30T00:00:00Z",
                    "finished_at_utc": "2026-06-30T00:00:01Z",
                    "work_dir": str(worker_dir),
                    "artifact_path": str(worker_dir / "source.pkl"),
                }
            ),
            encoding="utf-8",
        )
        (output_dir / "diagnostics.jsonl").write_text(
            json.dumps(
                {
                    "step": "build_altaz_grid",
                    "status": "ok",
                    "message": "Cloud source alt/az grid built",
                    "fields": {},
                }
            )
            + "\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(diagnose_cloud_source.subprocess, "run", fake_run)

    rc = diagnose_cloud_source.main(
        [
            "--output-dir",
            str(output_dir),
            "--request-id",
            "7",
            "--lat",
            "33.660109",
            "--lon",
            "-84.4102046",
            "--no-grid",
        ]
    )

    assert rc == 0
    cmd = commands[0]
    assert "zstarview.clouddisc.workers.cloud_source_worker" in cmd
    assert cmd[cmd.index("--work-dir") + 1] == str(output_dir / "worker")
    assert cmd[cmd.index("--cache-dir") + 1] == str(output_dir / "cache")
    assert cmd[cmd.index("--diagnostic-jsonl") + 1] == str(output_dir / "diagnostics.jsonl")
    assert "--skip-altaz-grid" in cmd
    assert envs[0]["ZSTARVIEW_DISABLE_FILE_LOGGING"] == "1"
    assert (output_dir / "diagnose.json").exists()
    assert "Diagnostic status: succeeded" in capsys.readouterr().out


def test_diagnose_cloud_source_rejects_nonempty_output_dir(tmp_path: Path) -> None:
    output_dir = tmp_path / "diagnosis"
    output_dir.mkdir()
    (output_dir / "old.txt").write_text("old", encoding="utf-8")

    try:
        diagnose_cloud_source.main(
            [
                "--output-dir",
                str(output_dir),
                "--lat",
                "33.0",
                "--lon",
                "-84.0",
            ]
        )
    except SystemExit as exc:
        assert "not empty" in str(exc)
    else:
        raise AssertionError("Expected SystemExit")
