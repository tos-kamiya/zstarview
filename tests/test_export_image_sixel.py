from __future__ import annotations

import subprocess
from types import SimpleNamespace
from io import BytesIO, StringIO

from astropy.time import Time
import pytest
from PySide6.QtGui import QImage

import zstarview.cli.export_image as mod


def test_require_img2sixel_binary_rejects_missing_command(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(mod.shutil, "which", lambda name: None)

    with pytest.raises(SystemExit):
        mod._require_img2sixel_binary()


def test_response_indicates_sixel_support_detects_code_four() -> None:
    assert mod._response_indicates_sixel_support(b"\x1b[?1;2;4c") is True
    assert mod._response_indicates_sixel_support(b"\x1b[?1;2c") is False


def test_require_sixel_terminal_support_accepts_device_attributes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("TERM", "xterm-256color")

    events: list[object] = []

    def fake_open(path: str, flags: int) -> int:
        events.append(("open", path, flags))
        return 3

    def fake_tcgetattr(fd: int) -> list[object]:
        assert fd == 3
        return [0, 0, 0, 0xFFFF, 0, 0, [0] * 8]

    def fake_tcsetattr(fd: int, when: int, attrs: list[object]) -> None:
        assert fd == 3
        events.append(
            (
                "tcsetattr",
                when,
                attrs[3],
                attrs[6][mod.termios.VMIN],
                attrs[6][mod.termios.VTIME],
            )
        )

    def fake_select(
        read_fds: list[int],
        _write_fds: list[int],
        _error_fds: list[int],
        _timeout: float,
    ):
        assert read_fds == [3]
        return ([3], [], [])

    read_calls = {"count": 0}

    def fake_read(fd: int, size: int) -> bytes:
        assert fd == 3
        assert size == 1024
        read_calls["count"] += 1
        return b"\x1b[?1;2;4c" if read_calls["count"] == 1 else b""

    def fake_write(fd: int, data: bytes) -> int:
        assert fd == 3
        events.append(("write", data))
        return len(data)

    def fake_close(fd: int) -> None:
        events.append(("close", fd))

    monkeypatch.setattr(mod.os, "open", fake_open)
    monkeypatch.setattr(mod.os, "write", fake_write)
    monkeypatch.setattr(mod.os, "read", fake_read)
    monkeypatch.setattr(mod.os, "close", fake_close)
    monkeypatch.setattr(mod.select, "select", fake_select)
    monkeypatch.setattr(mod.termios, "tcgetattr", fake_tcgetattr)
    monkeypatch.setattr(mod.termios, "tcsetattr", fake_tcsetattr)

    mod._require_sixel_terminal_support()

    assert events[0][0] == "open"
    assert ("write", b"\x1b[c") in events
    assert any(
        event[0] == "tcsetattr"
        and event[2] == (0xFFFF & ~(mod.termios.ECHO | mod.termios.ICANON))
        for event in events
        if event[0] == "tcsetattr"
    )
    assert any(
        event[0] == "tcsetattr" and event[2] == 0xFFFF
        for event in events
        if event[0] == "tcsetattr"
    )
    assert ("close", 3) in events


def test_require_sixel_terminal_support_rejects_non_sixel_response(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("TERM", "xterm-256color")

    def fake_open(path: str, flags: int) -> int:
        return 3

    def fake_tcgetattr(fd: int) -> list[object]:
        return [0, 0, 0, 0xFFFF, 0, 0, [0] * 8]

    def fake_tcsetattr(fd: int, when: int, attrs: list[object]) -> None:
        return None

    def fake_select(
        read_fds: list[int],
        _write_fds: list[int],
        _error_fds: list[int],
        _timeout: float,
    ):
        return ([3], [], [])

    def fake_read(fd: int, size: int) -> bytes:
        return b"\x1b[?1;2c"

    def fake_write(fd: int, data: bytes) -> int:
        return len(data)

    monkeypatch.setattr(mod.os, "open", fake_open)
    monkeypatch.setattr(mod.os, "write", fake_write)
    monkeypatch.setattr(mod.os, "read", fake_read)
    monkeypatch.setattr(mod.os, "close", lambda fd: None)
    monkeypatch.setattr(mod.select, "select", fake_select)
    monkeypatch.setattr(mod.termios, "tcgetattr", fake_tcgetattr)
    monkeypatch.setattr(mod.termios, "tcsetattr", fake_tcsetattr)

    with pytest.raises(SystemExit):
        mod._require_sixel_terminal_support()


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

    def fake_run(
        cmd: list[str], *, input: bytes, stdout: int, stderr: int, check: bool
    ) -> subprocess.CompletedProcess[bytes]:
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

    def fake_run(
        cmd: list[str], *, input: bytes, stdout: int, stderr: int, check: bool
    ) -> subprocess.CompletedProcess[bytes]:
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


def test_write_export_overlay_summary_to_stderr_emits_gui_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
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
    target = mod.SearchJumpTarget(
        label="Ceres",
        kind="jpl_small_body",
        sort_key=(0.0, "Ceres"),
        object_key="2000001",
        alt_deg=12.25,
        az_deg=34.5,
    )

    mod._write_export_overlay_summary_to_stderr(
        viewer_data=viewer,
        celestial_data=type(
            "DummyCelestial",
            (),
            {"time": Time("2026-02-27T00:00:00", format="isot", scale="utc")},
        )(),
        vmag_limit=6.0,
        search_overlay_target=target,
    )

    text = stderr_buffer.getvalue()
    assert text.startswith("Tokyo\n")
    assert "Tower height 634 m\n" in text
    assert "Observer height 12 m\n" in text
    assert "2026-02-27 00:00:00 UTC\n" in text
    assert "Alt 45°  Az 180° (S)\n" in text
    assert "Vmag limit 6.0\n" in text
    assert (
        "Search target label=Ceres | id=2000001 | kind=jpl_small_body | alt=12.2 deg | "
        "az=34.5 deg\n"
    ) in text
    assert text.endswith("az=34.5 deg\n")


def test_main_writes_overlay_summary_before_sixel(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events: list[str] = []
    viewer = mod.ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    viewer._search_overlay_target = mod.SearchJumpTarget(
        label="Ceres",
        kind="jpl_small_body",
        sort_key=(0.0, "Ceres"),
        object_key="2000001",
    )
    catalogs = SimpleNamespace(
        star_catalog_np=object(),
        star_catalog_lod6_indices=object(),
        star_catalog_meta=None,
        dso_catalog_np=None,
    )
    user_options = SimpleNamespace(
        vmag_limit=6.0,
        sky_disc_alpha=0.0,
        cloud_disc_alpha=0.0,
        satellite_opacity=0.0,
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
    monkeypatch.setattr(
        mod,
        "_build_window_inputs_from_args",
        lambda _args: (catalogs, viewer, user_options, runtime_options),
    )
    monkeypatch.setattr(
        mod,
        "setup_app",
        lambda _name: SimpleNamespace(setQuitOnLastWindowClosed=lambda _flag: None),
    )
    monkeypatch.setattr(mod, "_load_fonts", lambda: (object(), object()))
    monkeypatch.setattr(mod, "_build_compositor", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(mod, "_require_sixel_terminal_support", lambda: None)
    monkeypatch.setattr(
        mod,
        "compute_sky_snapshot",
        lambda **_kwargs: {
            "celestial": SimpleNamespace(
                time=Time("2026-02-27T00:00:00", format="isot", scale="utc")
            ),
            "sky_disc": None,
        },
    )
    monkeypatch.setattr(
        mod, "_build_render_style", lambda **_kwargs: SimpleNamespace(vmag_limit=6.0)
    )
    monkeypatch.setattr(mod, "_render_image", lambda **_kwargs: object())
    monkeypatch.setattr(
        mod,
        "_write_export_overlay_summary_to_stderr",
        lambda **kwargs: events.append(
            f"summary:{getattr(kwargs['search_overlay_target'], 'label', '')}"
        ),
    )
    monkeypatch.setattr(
        mod,
        "_write_sixel_to_stdout",
        lambda _image, *, img2sixel_bin: (
            events.append(f"sixel:{img2sixel_bin}") or True
        ),
    )

    mod.main()

    assert events == ["summary:Ceres", "sixel:/usr/bin/img2sixel"]


def test_main_rejects_unsupported_sixel_terminal_before_loading_inputs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events: list[str] = []
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
    monkeypatch.setattr(
        mod,
        "_require_sixel_terminal_support",
        lambda: (_ for _ in ()).throw(SystemExit(1)),
    )
    monkeypatch.setattr(
        mod, "_build_window_inputs_from_args", lambda _args: events.append("build")
    )

    with pytest.raises(SystemExit):
        mod.main()

    assert events == []
