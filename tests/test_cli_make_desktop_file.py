from __future__ import annotations

from pathlib import Path

from zstarview.cli import make_desktop_file


def test_make_desktop_file_writes_gui_launcher_in_current_directory(
    monkeypatch,
    capsys,
    tmp_path: Path,
) -> None:
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("sys.argv", ["zstarview-make-desktop-file"])

    make_desktop_file.main()

    out_path = tmp_path / "zstarview-gui.desktop"
    content = out_path.read_text(encoding="utf-8")

    assert out_path.is_file()
    assert "Exec=zstarview-gui" in content
    assert "StartupWMClass=zstarview" in content
    assert capsys.readouterr().out.strip() == f"Wrote: {out_path}"


def test_make_desktop_file_write_installs_gui_launcher(
    monkeypatch,
    capsys,
    tmp_path: Path,
) -> None:
    applications_dir = tmp_path / ".local" / "share" / "applications"
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr("sys.argv", ["zstarview-make-desktop-file", "--write"])
    monkeypatch.setattr(make_desktop_file.os, "system", lambda _cmd: 0)

    make_desktop_file.main()

    out_path = applications_dir / "zstarview-gui.desktop"
    content = out_path.read_text(encoding="utf-8")

    assert out_path.is_file()
    assert "Exec=zstarview-gui" in content
    assert capsys.readouterr().out.strip() == f"Wrote: {out_path}"
