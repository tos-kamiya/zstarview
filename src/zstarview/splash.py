import logging
import threading
import sys
from pathlib import Path
from typing import Callable, Tuple

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QGuiApplication, QIcon, QLinearGradient, QPainter, QPixmap
from PySide6.QtWidgets import QApplication, QSplashScreen

from .__about__ import __version__
from .paths import APP_AUTHOR, APP_DISPLAY_NAME, APP_ICON_FILE, APP_ID, CACHE_PATH

logger = logging.getLogger(__name__)

_cache_path = Path(CACHE_PATH)
_cache_path.mkdir(parents=True, exist_ok=True)

SPLASH_INFO_COLOR = QColor(18, 29, 48)
SPLASH_WARN_COLOR = QColor(130, 82, 20)
SPLASH_ERROR_COLOR = QColor(146, 34, 34)


def _get_splash_palette(visual_preset: str) -> tuple[list[QColor], QColor, QColor]:
    """Return gradient colors, frame color, and default message color for splash."""
    if visual_preset == "night":
        return (
            [QColor(12, 14, 20), QColor(8, 10, 14), QColor(4, 6, 9)],
            QColor(70, 76, 92),
            QColor(228, 236, 250),
        )
    if visual_preset == "black":
        return (
            [QColor(6, 6, 6), QColor(3, 3, 3), QColor(0, 0, 0)],
            QColor(56, 56, 64),
            QColor(244, 248, 255),
        )
    if visual_preset == "white":
        return (
            [QColor(252, 252, 252), QColor(234, 234, 234), QColor(206, 206, 206)],
            QColor(158, 178, 206),
            QColor(19, 31, 50),
        )
    return (
        [QColor(240, 248, 255), QColor(226, 240, 252), QColor(206, 228, 246)],
        QColor(158, 182, 206),
        QColor(18, 29, 48),
    )


def _build_splash_pixmap(
    visual_preset: str, width: int = 400, height: int = 200
) -> QPixmap:
    """Create a splash background matching the selected visual preset."""
    grad_colors, frame_color, _ = _get_splash_palette(visual_preset)
    pixmap = QPixmap(width, height)
    pixmap.fill(Qt.GlobalColor.transparent)
    painter = QPainter(pixmap)
    grad = QLinearGradient(0, 0, width, height)
    grad.setColorAt(0.0, grad_colors[0])
    grad.setColorAt(0.55, grad_colors[1])
    grad.setColorAt(1.0, grad_colors[2])
    painter.fillRect(0, 0, width, height, grad)
    painter.setPen(frame_color)
    painter.drawRect(0, 0, width - 1, height - 1)
    painter.end()
    return pixmap


class SplashLogHandler(logging.Handler):
    """A temporary log handler to display logs on the splash screen."""

    def __init__(self, show_fn: Callable[[str, QColor], None], info_color: QColor):
        super().__init__()
        self.show_fn = show_fn
        self._main_thread_id = threading.get_ident()
        self._info_color = info_color
        self.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))

    def emit(self, record: logging.LogRecord) -> None:
        try:
            if threading.get_ident() != self._main_thread_id:
                return
            msg = self.format(record)
            color = (
                self._info_color
                if record.levelno < logging.WARNING
                else SPLASH_WARN_COLOR
                if record.levelno < logging.ERROR
                else SPLASH_ERROR_COLOR
            )
            self.show_fn(msg, color)
        except Exception:
            self.handleError(record)


def setup_app(app_name: str = APP_DISPLAY_NAME) -> QApplication:
    """Initialize and configure the QApplication instance."""
    app = QApplication(sys.argv)
    QGuiApplication.setDesktopFileName(APP_ID)
    app.setApplicationName(app_name)
    app.setApplicationDisplayName(app_name)
    app.setOrganizationName(APP_AUTHOR)
    app.setApplicationVersion(__version__)
    app.setWindowIcon(QIcon(APP_ICON_FILE))
    return app


def setup_splash_and_attach_logger(
    app: QApplication,
    app_name: str,
    root_logger: logging.Logger,
    visual_preset: str,
) -> Tuple[QSplashScreen, SplashLogHandler, Callable[[str], None]]:
    """Create a splash screen and attach a temporary log handler."""
    splash = QSplashScreen(QPixmap(400, 200), Qt.WindowType.WindowStaysOnTopHint)
    _, _, info_color = _get_splash_palette(visual_preset)
    splash.setPixmap(_build_splash_pixmap(visual_preset, 400, 200))
    splash.show()

    context_line = ""

    def set_splash_context(line: str) -> None:
        nonlocal context_line
        context_line = line.strip()

    def show_splash_message(message: str, color: QColor) -> None:
        header = f"{app_name} ver. {__version__}"
        if context_line:
            header += f"\n{context_line}"
        splash.showMessage(
            f"{header}\n{message}",
            Qt.AlignmentFlag.AlignCenter,
            color,
        )
        app.processEvents()

    splash_handler = SplashLogHandler(show_splash_message, info_color)
    root_logger.addHandler(splash_handler)
    return splash, splash_handler, set_splash_context
