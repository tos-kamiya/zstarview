from __future__ import annotations

import subprocess
from types import SimpleNamespace
from io import BytesIO, StringIO

from astropy.time import Time
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


def test_write_export_overlay_summary_to_stderr_emits_gui_metadata(monkeypatch: pytest.MonkeyPatch) -> None:
    stderr_buffer = StringIO()
    monkeypatch.setattr(mod.sys, "stderr", stderr_buffer)

    viewer = mod.ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=12.0,
        location_height_label="Tower height",
        location_height_m=634.0,
        show_observer_height=True,
    )

    mod._write_export_overlay_summary_to_stderr(
        viewer_data=viewer,
        celestial_data=type("DummyCelestial", (), {"time": Time("2026-02-27T00:00:00", format="isot", scale="utc")})(),
        vmag_limit=6.0,
    )

    text = stderr_buffer.getvalue()
    assert "Tokyo\n" in text
    assert "Tower height 634 m\n" in text
    assert "Observer height 12 m\n" in text
    assert "Alt 45°  Az 180° (S)\n" in text
    assert text.endswith("Vmag limit 6.0\n")


def test_main_writes_overlay_summary_before_sixel(monkeypatch: pytest.MonkeyPatch) -> None:
    events: list[str] = []
    viewer = mod.ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    catalogs = SimpleNamespace(
        star_catalog_lod6_np=object(),
        star_catalog_full_np=object(),
        dso_catalog_np=None,
    )
    user_options = SimpleNamespace(
        vmag_limit=6.0,
        sky_disc_alpha=0.0,
        cloud_disc_alpha=0.0,
        terrain_horizon_opacity=0.0,
        urban_outline_opacity=0.0,
        aircraft_opacity=0.0,
        visual_preset="night",
    )
    runtime_options = SimpleNamespace(delta_t=0.0)

    monkeypatch.setattr(
        mod,
        "parse_export_image_args",
        lambda: SimpleNamespace(
            sixel=True,
            output=None,
            image_size=(4, 4),
            layer_timeout_seconds=30.0,
            allow_partial_data=False,
        ),
    )
    monkeypatch.setattr(mod, "setup_root_logger", lambda: None)
    monkeypatch.setattr(mod, "_require_img2sixel_binary", lambda: "/usr/bin/img2sixel")
    monkeypatch.setattr(mod, "_build_window_inputs_from_args", lambda _args: (catalogs, viewer, user_options, runtime_options))
    monkeypatch.setattr(mod, "setup_app", lambda _name: SimpleNamespace(setQuitOnLastWindowClosed=lambda _flag: None))
    monkeypatch.setattr(mod, "_load_fonts", lambda: (object(), object()))
    monkeypatch.setattr(mod, "_build_compositor", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(
        mod,
        "compute_sky_snapshot",
        lambda **_kwargs: {"celestial": SimpleNamespace(time=Time("2026-02-27T00:00:00", format="isot", scale="utc")), "sky_disc": None},
    )
    monkeypatch.setattr(mod, "_build_render_style", lambda **_kwargs: SimpleNamespace(vmag_limit=6.0))
    monkeypatch.setattr(mod, "_render_image", lambda **_kwargs: object())
    monkeypatch.setattr(
        mod,
        "_write_export_overlay_summary_to_stderr",
        lambda **_kwargs: events.append("summary"),
    )
    monkeypatch.setattr(
        mod,
        "_write_sixel_to_stdout",
        lambda _image, *, img2sixel_bin: events.append(f"sixel:{img2sixel_bin}") or True,
    )

    mod.main()

    assert events == ["summary", "sixel:/usr/bin/img2sixel"]
