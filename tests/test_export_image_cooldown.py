from __future__ import annotations

import json
import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

import zstarview.cli.export_image as export_image
import zstarview.cli.export_image_cooldown as cooldown
from zstarview.cli.args import parse_export_image_args


def test_recent_start_rejected_and_old_start_accepted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(cooldown.time, "time", lambda: 1000.0)
    with cooldown.export_cooldown(tmp_path):
        pass
    with pytest.raises(SystemExit, match="retry in 60 seconds"):
        with cooldown.export_cooldown(tmp_path):
            pytest.fail("must not start")
    monkeypatch.setattr(cooldown.time, "time", lambda: 1060.0)
    with cooldown.export_cooldown(tmp_path):
        pass


def test_bypass_leaves_record_untouched(tmp_path: Path) -> None:
    record = tmp_path / cooldown.EXPORT_COOLDOWN_FILE
    record.write_text("existing", encoding="utf-8")
    with cooldown.export_cooldown(tmp_path, bypass=True):
        pass
    assert record.read_text(encoding="utf-8") == "existing"
    args = parse_export_image_args(["--debug-ignore-export-cooldown", "-o", "test.png"])
    assert args.debug_ignore_export_cooldown


def test_changed_record_warns_even_when_export_fails(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.WARNING), pytest.raises(RuntimeError, match="failure"):
        with cooldown.export_cooldown(tmp_path):
            path = tmp_path / cooldown.EXPORT_COOLDOWN_FILE
            record = json.loads(path.read_text(encoding="utf-8"))
            record["token"] = "another-instance"
            path.write_text(json.dumps(record), encoding="utf-8")
            raise RuntimeError("failure")
    assert "may have overlapped" in caplog.text


def test_invalid_or_future_record_does_not_block(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(cooldown.time, "time", lambda: 1000.0)
    path = tmp_path / cooldown.EXPORT_COOLDOWN_FILE
    path.write_bytes(b"\xff")
    with cooldown.export_cooldown(tmp_path):
        pass
    path.write_text(json.dumps({"started_at": 2000, "token": "x"}), encoding="utf-8")
    with cooldown.export_cooldown(tmp_path):
        pass


def test_main_runs_inside_cooldown_and_debug_bypasses(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(export_image, "CACHE_PATH", tmp_path)
    args = SimpleNamespace(print_cache_dir=False, debug_ignore_export_cooldown=False)
    monkeypatch.setattr(export_image, "parse_export_image_args", lambda: args)
    calls: list[object] = []
    monkeypatch.setattr(export_image, "_run_export", calls.append)
    export_image.main()
    with pytest.raises(SystemExit, match="started recently"):
        export_image.main()
    args.debug_ignore_export_cooldown = True
    export_image.main()
    assert len(calls) == 2
