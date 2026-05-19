from __future__ import annotations

import logging

from PySide6.QtWidgets import QApplication, QWidget

from zstarview.gui.window import StartupLogOverlay
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


def test_startup_log_overlay_append_line_does_not_process_events(monkeypatch) -> None:
    parent = QWidget()
    overlay = StartupLogOverlay(parent)
    calls: list[int] = []

    def _process_events() -> None:
        calls.append(1)

    monkeypatch.setattr(QApplication, "processEvents", staticmethod(_process_events))

    overlay.append_line("startup line", logging.INFO)

    assert calls == []
