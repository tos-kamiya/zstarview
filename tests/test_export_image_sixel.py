from __future__ import annotations

import json
import subprocess
from io import BytesIO, StringIO
from types import SimpleNamespace

import pytest
from astropy.time import Time
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
        cmd: list[str], *, input: bytes, capture_output: bool, check: bool
    ) -> subprocess.CompletedProcess[bytes]:
        captured["cmd"] = cmd
        captured["input"] = input
        captured["capture_output"] = capture_output
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
        cmd: list[str], *, input: bytes, capture_output: bool, check: bool
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


def test_encode_image_as_png_bytes_embeds_export_metadata() -> None:
    image = QImage(3, 3, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0xFFABCDEF)
    viewer = mod.ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Matsue",
        view_center=(35.0, 120.0),
        observer_height_m=12.0,
        height_add_m=1.7,
        ground_elevation_m=35.0,
        location_height_label="Building",
        location_height_m=20.0,
    )
    place_location = mod.ResolvedLocation(
        display_name="Matsue, Shimane, Japan",
        lat=35.47,
        lon=133.05,
        tz="Asia/Tokyo",
        persistence_key="matsue",
        observer_height_m=1.7,
        kind="place",
        persistence_value={"resolver": "nominatim", "query": "Matsue"},
        ground_elevation_m=35.0,
        location_height_label=None,
        location_height_m=0.0,
        height_add_m=1.7,
        cc="JP",
    )
    target = mod.SearchJumpTarget(
        label="Ceres",
        kind="jpl_small_body",
        sort_key=(0.0, "Ceres"),
        object_key="2000001",
        alt_deg=12.25,
        az_deg=34.5,
    )
    payload = mod._build_export_image_metadata_payload(
        app_version="1.31.18",
        viewer_data=viewer,
        celestial_data=SimpleNamespace(
            time=Time("2026-02-27T00:00:00", format="isot", scale="utc"),
            planets=[],
        ),
        style=SimpleNamespace(vmag_limit=6.0),
        place_query="Matsue",
        place_location=place_location,
        search_overlay_target=target,
        cloud_coverage_ratio=0.625,
        urban_outline_source="PLATEAU",
        urban_outline_count=10301,
    )

    png_bytes = mod._encode_image_as_png_bytes(image, metadata_payload=payload)
    loaded = QImage()
    assert loaded.loadFromData(png_bytes, "PNG")
    assert mod.EXPORT_IMAGE_METADATA_TEXT_KEY in loaded.textKeys()
    loaded_payload = json.loads(loaded.text(mod.EXPORT_IMAGE_METADATA_TEXT_KEY))
    assert loaded_payload == payload


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
        height_add_m=12.0,
        ground_elevation_m=35.0,
        location_height_label="Tower height",
        location_height_m=634.0,
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
        cloud_coverage_ratio=0.625,
        search_overlay_target=target,
    )

    text = stderr_buffer.getvalue()
    assert text.startswith(
        "Tokyo\nLat: 35.00000, Lon: 139.00000 | Height: ground 35 m, building 634 m, add 12 m\n"
    )
    assert "2026-02-27 00:00:00 UTC\n" in text
    assert "Alt 45°  Az 180° (S)\n" in text
    assert "Vmag limit 6.0\n" in text
    assert "Cloud coverage 62.5%\n" in text
    assert (
        "Search target label=Ceres | id=2000001 | kind=jpl_small_body | alt=12.2 deg | "
        "az=34.5 deg\n"
    ) in text
    assert text.index("Vmag limit 6.0\n") < text.index("Cloud coverage 62.5%\n")
    assert text.index("Cloud coverage 62.5%\n") < text.index(
        "Search target label=Ceres | id=2000001 | kind=jpl_small_body | alt=12.2 deg | "
        "az=34.5 deg\n"
    )
    assert text.endswith("az=34.5 deg\n")


def test_main_writes_overlay_summary_before_sixel(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events: list[str] = []
    search_overlay_target = mod.SearchJumpTarget(
        label="Ceres",
        kind="jpl_small_body",
        sort_key=(0.0, "Ceres"),
        object_key="2000001",
    )
    viewer = mod.ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
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
        sky_disc_style="smooth",
        cloud_disc_alpha=0.0,
        geo_satellite=False,
        water_overlay_opacity=0.0,
        satellite_opacity=0.0,
        terrain_horizon_opacity=0.0,
        urban_outline_opacity=0.0,
        aircraft_opacity=0.0,
        overlay_font_size=11,
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
            include_direction_grid=False,
            print_cache_dir=False,
            clear_long_lived_cache=False,
            sky_disc_style="smooth",
        ),
    )
    monkeypatch.setattr(mod, "setup_root_logger", lambda: None)
    monkeypatch.setattr(mod, "_require_img2sixel_binary", lambda: "/usr/bin/img2sixel")
    monkeypatch.setattr(
        mod,
        "_build_window_inputs_from_args",
        lambda _args: (
            catalogs,
            viewer,
            user_options,
            runtime_options,
            search_overlay_target,
            None,
        ),
    )
    monkeypatch.setattr(
        mod,
        "setup_app",
        lambda _name: SimpleNamespace(setQuitOnLastWindowClosed=lambda _flag: None),
    )
    monkeypatch.setattr(mod, "_load_fonts", lambda *_args: (object(), object()))
    monkeypatch.setattr(mod, "_build_compositor", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(mod, "_require_sixel_terminal_support", lambda: None)
    monkeypatch.setattr(
        mod,
        "compute_sky_snapshot",
        lambda **_kwargs: {
            "celestial": SimpleNamespace(
                time=Time("2026-02-27T00:00:00", format="isot", scale="utc"),
                planets=[SimpleNamespace(name="sun", alt=-10.0)],
            ),
            "sky_disc": None,
        },
    )
    monkeypatch.setattr(
        mod, "_build_render_style", lambda **_kwargs: SimpleNamespace(vmag_limit=6.0)
    )
    monkeypatch.setattr(mod, "_fetch_water_overlay_dots_layer", lambda **_kwargs: None)
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
        lambda _image, *, img2sixel_bin, **_kwargs: (
            events.append(f"sixel:{img2sixel_bin}") or True
        ),
    )

    mod.main()

    assert events == ["summary:Ceres", "sixel:/usr/bin/img2sixel"]


def test_main_aborts_when_cloud_layer_is_unavailable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events: list[str] = []
    viewer = mod.ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
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
        sky_disc_style="smooth",
        cloud_disc_alpha=0.25,
        geo_satellite=False,
        water_overlay_opacity=0.0,
        satellite_opacity=0.0,
        terrain_horizon_opacity=0.0,
        urban_outline_opacity=0.0,
        aircraft_opacity=0.0,
        overlay_font_size=11,
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
            include_direction_grid=False,
            print_cache_dir=False,
            clear_long_lived_cache=False,
            sky_disc_style="smooth",
        ),
    )
    monkeypatch.setattr(mod, "setup_root_logger", lambda: None)
    monkeypatch.setattr(mod, "_require_img2sixel_binary", lambda: "/usr/bin/img2sixel")
    monkeypatch.setattr(mod, "_require_sixel_terminal_support", lambda: None)
    monkeypatch.setattr(
        mod,
        "_build_window_inputs_from_args",
        lambda _args: (
            catalogs,
            viewer,
            user_options,
            runtime_options,
            None,
            None,
        ),
    )
    monkeypatch.setattr(
        mod,
        "setup_app",
        lambda _name: SimpleNamespace(setQuitOnLastWindowClosed=lambda _flag: None),
    )
    monkeypatch.setattr(mod, "_load_fonts", lambda *_args: (object(), object()))
    monkeypatch.setattr(mod, "_build_compositor", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(
        mod,
        "compute_sky_snapshot",
        lambda **_kwargs: {
            "celestial": SimpleNamespace(
                time=Time("2026-02-27T00:00:00", format="isot", scale="utc"),
                planets=[SimpleNamespace(name="sun", alt=-10.0)],
            ),
            "sky_disc": None,
        },
    )
    monkeypatch.setattr(
        mod,
        "_fetch_cloud_layer",
        lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("No supported satellite for this region")),
    )
    monkeypatch.setattr(mod, "_fetch_water_overlay_dots_layer", lambda **_kwargs: None)
    monkeypatch.setattr(
        mod,
        "_build_render_style",
        lambda **_kwargs: SimpleNamespace(vmag_limit=6.0),
    )
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

    with pytest.raises(SystemExit):
        mod.main()

    assert events == []


def test_main_reports_partial_data_note_when_terrain_layer_aborts(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    messages: list[str] = []
    viewer = mod.ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
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
        sky_disc_style="smooth",
        cloud_disc_alpha=0.0,
        geo_satellite=False,
        water_overlay_opacity=0.0,
        satellite_opacity=0.0,
        terrain_horizon_opacity=0.05,
        urban_outline_opacity=0.0,
        aircraft_opacity=0.0,
        overlay_font_size=11,
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
            include_direction_grid=False,
            print_cache_dir=False,
            clear_long_lived_cache=False,
            sky_disc_style="smooth",
        ),
    )
    monkeypatch.setattr(mod, "setup_root_logger", lambda: None)
    monkeypatch.setattr(mod, "_require_img2sixel_binary", lambda: "/usr/bin/img2sixel")
    monkeypatch.setattr(mod, "_require_sixel_terminal_support", lambda: None)
    monkeypatch.setattr(
        mod,
        "_build_window_inputs_from_args",
        lambda _args: (
            catalogs,
            viewer,
            user_options,
            runtime_options,
            None,
            None,
        ),
    )
    monkeypatch.setattr(
        mod,
        "setup_app",
        lambda _name: SimpleNamespace(setQuitOnLastWindowClosed=lambda _flag: None),
    )
    monkeypatch.setattr(mod, "_load_fonts", lambda *_args: (object(), object()))
    monkeypatch.setattr(mod, "_build_compositor", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(
        mod,
        "compute_sky_snapshot",
        lambda **_kwargs: {
            "celestial": SimpleNamespace(
                time=Time("2026-02-27T00:00:00", format="isot", scale="utc"),
                planets=[SimpleNamespace(name="sun", alt=-10.0)],
            ),
            "sky_disc": None,
        },
    )
    monkeypatch.setattr(
        mod,
        "_fetch_terrain_horizon_layer",
        lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("terrain timed out")),
    )
    monkeypatch.setattr(
        mod.logger,
        "error",
        lambda message, *args, **kwargs: messages.append(
            message % args if args else message
        ),
    )

    with pytest.raises(SystemExit):
        mod.main()

    assert any("--allow-partial-data" in message for message in messages)


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
            include_direction_grid=False,
            print_cache_dir=False,
            clear_long_lived_cache=False,
            sky_disc_style="smooth",
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
        mod,
        "_build_window_inputs_from_args",
        lambda _args: events.append("build"),
    )

    with pytest.raises(SystemExit):
        mod.main()

    assert events == []
