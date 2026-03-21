from __future__ import annotations

import subprocess
from io import BytesIO

import pytest
from PySide6.QtGui import QImage

import zstarview.export_image as mod


def test_require_img2sixel_binary_rejects_missing_command(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(mod.shutil, "which", lambda name: None)

    with pytest.raises(SystemExit):
        mod._require_img2sixel_binary()


def test_write_sixel_to_stdout_pipes_png_bytes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    image = QImage(4, 4, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0xFF112233)
    captured: dict[str, object] = {}
    stdout_bytes = BytesIO()

    class FakeStdout:
        def __init__(self) -> None:
            self.buffer = stdout_bytes

    monkeypatch.setattr(mod.sys, "stdout", FakeStdout())

    def fake_run(cmd: list[str], *, input: bytes, stdout: int, stderr: int, check: bool) -> subprocess.CompletedProcess[bytes]:
        captured["cmd"] = cmd
        captured["input"] = input
        captured["stdout"] = stdout
        captured["stderr"] = stderr
        captured["check"] = check
        return subprocess.CompletedProcess(cmd, 0, b"\x1bPqSIXEL\x1b\\", b"")

    monkeypatch.setattr(mod.subprocess, "run", fake_run)

    assert mod._write_sixel_to_stdout(image, img2sixel_bin="/usr/bin/img2sixel") is True

    assert stdout_bytes.getvalue() == b"\x1bPqSIXEL\x1b\\"
    assert captured["cmd"] == ["/usr/bin/img2sixel", "-"]
    assert isinstance(captured["input"], bytes)
    assert bytes(captured["input"]).startswith(b"\x89PNG\r\n\x1a\n")


def test_write_sixel_to_stdout_reports_failed_command(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    image = QImage(2, 2, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0xFF000000)

    def fake_run(cmd: list[str], *, input: bytes, stdout: int, stderr: int, check: bool) -> subprocess.CompletedProcess[bytes]:
        return subprocess.CompletedProcess(cmd, 1, b"", b"img2sixel failed")

    monkeypatch.setattr(mod.subprocess, "run", fake_run)

    assert mod._write_sixel_to_stdout(image, img2sixel_bin="img2sixel") is False


def test_write_png_to_stdout_writes_png_bytes(monkeypatch: pytest.MonkeyPatch) -> None:
    image = QImage(3, 3, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0xFFABCDEF)
    stdout_bytes = BytesIO()

    class FakeStdout:
        def __init__(self) -> None:
            self.buffer = stdout_bytes

    monkeypatch.setattr(mod.sys, "stdout", FakeStdout())

    assert mod._write_png_to_stdout(image) is True
    assert stdout_bytes.getvalue().startswith(b"\x89PNG\r\n\x1a\n")
