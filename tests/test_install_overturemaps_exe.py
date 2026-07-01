from __future__ import annotations

import importlib
from pathlib import Path


def _load_module():
    return importlib.import_module("zstarview.cli.install_overturemaps_exe")


def test_stage_overturemaps_executable_copies_and_overwrites(tmp_path: Path, monkeypatch) -> None:
    mod = _load_module()
    source = tmp_path / "overturemaps-1.0.1-windows-x86_64.exe"
    source.write_text("first", encoding="utf-8")
    cache_root = tmp_path / "cache"
    staged_path = cache_root / "overturemaps.exe"
    staged_path.parent.mkdir(parents=True, exist_ok=True)
    staged_path.write_text("old", encoding="utf-8")
    monkeypatch.setattr(mod, "staged_overturemaps_executable_path", lambda: staged_path)

    got = mod.stage_overturemaps_executable(source)

    assert got == staged_path
    assert staged_path.read_text(encoding="utf-8") == "first"

    source.write_text("second", encoding="utf-8")
    mod.stage_overturemaps_executable(source)

    assert staged_path.read_text(encoding="utf-8") == "second"


def test_main_reports_missing_source(tmp_path: Path, monkeypatch, capsys) -> None:
    mod = _load_module()
    cache_root = tmp_path / "cache"
    staged_path = cache_root / "overturemaps.exe"
    monkeypatch.setattr(mod, "staged_overturemaps_executable_path", lambda: staged_path)

    rc = mod.main([str(tmp_path / "missing.exe")])

    assert rc == 1
    assert not staged_path.exists()
    captured = capsys.readouterr()
    assert "Source executable does not exist" in captured.err
