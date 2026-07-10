from __future__ import annotations

import logging
from types import SimpleNamespace

from PySide6.QtWidgets import QApplication, QWidget

import zstarview.gui.window as window_module
from zstarview.gui.window import SkyWindow, StartupLogOverlay
from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.startup_log import BufferedStartupLogHandler

_app = QApplication.instance() or QApplication([])


def test_buffered_startup_log_handler_flushes_pending_lines() -> None:
    handler = BufferedStartupLogHandler()
    logger = logging.getLogger("zstarview.tests.startup-log")
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False

    lines: list[tuple[str, int]] = []
    try:
        logger.info("first message")
        assert lines == []

        handler.set_consumer(lambda line, levelno: lines.append((line, levelno)))
        assert len(lines) == 1
        assert "first message" in lines[0][0]
        assert lines[0][1] == logging.INFO

        logger.info("second message")
        assert len(lines) == 2
        assert "second message" in lines[1][0]
        assert lines[1][1] == logging.INFO
    finally:
        logger.removeHandler(handler)
        logger.propagate = True


def test_startup_log_overlay_appends_lines() -> None:
    parent = QWidget()
    overlay = StartupLogOverlay(parent)

    overlay.append_line("startup line", logging.INFO)

    assert "startup line" in overlay.toPlainText()


def test_startup_log_overlay_colors_errors() -> None:
    parent = QWidget()
    overlay = StartupLogOverlay(parent)

    overlay.append_line("startup failed", logging.ERROR)

    html = overlay.toHtml()
    assert "startup failed" in html
    assert "#ff8080" in html


def test_startup_log_overlay_can_use_theme_text_color() -> None:
    parent = QWidget()
    theme = THEME_STYLES_BY_PRESET["day"]
    overlay = StartupLogOverlay(parent, text_rgb=theme.text.foreground_rgb)

    overlay.append_line("startup line", logging.INFO)
    html = overlay.toHtml().lower()

    red, green, blue = theme.text.foreground_rgb
    assert f"#{red:02x}{green:02x}{blue:02x}" in html


def test_startup_log_overlay_uses_splash_info_color_for_day_theme(monkeypatch) -> None:
    captured: dict[str, tuple[int, ...] | None] = {"text_rgb": None, "background_rgba": None}

    class _FakeOverlay:
        def __init__(self, _parent, *, text_rgb=None, background_rgba=None):
            captured["text_rgb"] = text_rgb
            captured["background_rgba"] = background_rgba
            self.text_rgb = text_rgb

    dummy = SimpleNamespace(
        visual_preset="day",
        _client_widget=QWidget(),
        _startup_log_overlay=None,
        _layout_startup_log_overlay=lambda: None,
    )
    monkeypatch.setattr(window_module, "StartupLogOverlay", _FakeOverlay)

    overlay = SkyWindow._ensure_startup_log_overlay(dummy)

    assert overlay is dummy._startup_log_overlay
    assert captured["text_rgb"] == THEME_STYLES_BY_PRESET["day"].splash.info_text_rgb
    assert captured["background_rgba"] == (
        *THEME_STYLES_BY_PRESET["day"].window_background.base_rgb,
        180,
    )


def test_atlas_startup_log_overlay_is_opaque(monkeypatch) -> None:
    captured: dict[str, tuple[int, ...] | None] = {"background_rgba": None}

    class _FakeOverlay:
        def __init__(self, _parent, *, text_rgb=None, background_rgba=None):
            captured["background_rgba"] = background_rgba

    dummy = SimpleNamespace(
        visual_preset="white",
        presentation_id="instrument",
        _client_widget=QWidget(),
        _startup_log_overlay=None,
        _layout_startup_log_overlay=lambda: None,
    )
    monkeypatch.setattr(window_module, "StartupLogOverlay", _FakeOverlay)

    SkyWindow._ensure_startup_log_overlay(dummy)

    assert captured["background_rgba"] == (
        *THEME_STYLES_BY_PRESET["white"].window_background.base_rgb,
        255,
    )


def test_startup_log_overlay_append_line_does_not_process_events(monkeypatch) -> None:
    parent = QWidget()
    overlay = StartupLogOverlay(parent)
    calls: list[int] = []

    def _process_events() -> None:
        calls.append(1)

    monkeypatch.setattr(QApplication, "processEvents", staticmethod(_process_events))

    overlay.append_line("startup line", logging.INFO)

    assert calls == []
