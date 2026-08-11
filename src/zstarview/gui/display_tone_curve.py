"""Display-only tone compensation and its visual calibration dialog."""

from __future__ import annotations

import argparse

import numpy as np
from PySide6.QtCore import QRect, Qt, Signal
from PySide6.QtGui import QColor, QImage, QMouseEvent, QPainter, QPaintEvent
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

DisplayToneCurve = tuple[int, int] | None
DISPLAY_TONE_CURVE_LUT_VERSION = 1
CALIBRATION_BLACK_VALUES = tuple(range(0, 33, 2))
CALIBRATION_WHITE_VALUES = tuple(range(223, 256, 2))


def parse_display_tone_curve(value: str) -> DisplayToneCurve:
    """Parse ``off`` or a validated ``BLACK,WHITE`` pair."""
    text = (value or "").strip().lower()
    if text == "off":
        return None
    parts = [part.strip() for part in text.split(",")]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            "Display tone curve must be 'off' or BLACK,WHITE."
        )
    try:
        black, white = (int(part) for part in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Display tone curve BLACK and WHITE must be integers."
        ) from exc
    if not 0 <= black < white <= 255:
        raise argparse.ArgumentTypeError(
            "Display tone curve must satisfy 0 <= BLACK < WHITE <= 255."
        )
    return black, white


def format_display_tone_curve(value: DisplayToneCurve) -> str:
    if value is None:
        return "off"
    return f"{value[0]},{value[1]}"


def build_display_tone_curve_lut(black: int, white: int) -> np.ndarray:
    """Build a monotone LUT which moves edge tones into the visible range."""
    if not 0 <= black < white <= 255:
        raise ValueError("display tone curve limits are invalid")
    output = np.arange(256, dtype=np.float64)

    def smooth_segment(x0: int, y0: int, x1: int, y1: int) -> None:
        indices = np.arange(x0, x1 + 1)
        t = (indices.astype(np.float64) - x0) / float(x1 - x0)
        smooth = t * t * (3.0 - 2.0 * t)
        output[indices] = float(y0) + float(y1 - y0) * smooth

    dark_join = max(64, black)
    light_join = min(191, white)
    if dark_join < light_join:
        smooth_segment(1, black, dark_join, dark_join)
        smooth_segment(light_join, light_join, 254, white)
    else:
        smooth_segment(1, black, 254, white)
    lut = np.rint(output).astype(np.uint8)
    lut[0] = 0
    lut[255] = 255
    return np.maximum.accumulate(lut)


def apply_display_tone_curve(image: QImage, curve: tuple[int, int]) -> QImage:
    """Return a detached image with the same LUT applied to RGB; preserve alpha."""
    lut = build_display_tone_curve_lut(*curve)
    rgba_image = image.convertToFormat(QImage.Format.Format_RGBA8888)
    height, width = rgba_image.height(), rgba_image.width()
    view = np.frombuffer(rgba_image.bits(), dtype=np.uint8).reshape(
        height, rgba_image.bytesPerLine()
    )
    pixels = view[:, : width * 4].reshape(height, width, 4)
    pixels[:, :, :3] = lut[pixels[:, :, :3]]
    return rgba_image


class DisplayToneCalibrationPattern(QWidget):
    """Uncorrected numbered grayscale patches on exact black/white fields."""

    selection_changed = Signal(int, int)

    def __init__(self, curve: tuple[int, int] = (12, 247), parent=None) -> None:
        super().__init__(parent)
        self.black, self.white = curve
        self._last_selected: str | None = None
        self.setMinimumSize(640, 400)
        self.setMouseTracking(True)

    def paintEvent(self, _event: QPaintEvent) -> None:
        painter = QPainter(self)
        half = self.height() // 2
        self._paint_section(
            painter,
            QRect(0, 0, self.width(), half),
            CALIBRATION_BLACK_VALUES,
            selected=self.black,
            instruction=(
                "BLACK: On the RGB 0,0,0 background, click the darkest (leftmost) "
                "numbered square whose boundary you can clearly distinguish."
            ),
            foreground=QColor(255, 255, 255),
            selection_color=QColor(0, 255, 255),
        )
        self._paint_section(
            painter,
            QRect(0, half, self.width(), self.height() - half),
            CALIBRATION_WHITE_VALUES,
            selected=self.white,
            instruction=(
                "WHITE: On the RGB 255,255,255 background, click the brightest "
                "(rightmost) numbered square whose boundary you can clearly distinguish."
            ),
            foreground=QColor(0, 0, 0),
            selection_color=QColor(180, 0, 180),
        )

    @staticmethod
    def _patch_rects(rect: QRect, values: tuple[int, ...]) -> list[QRect]:
        margin = 12
        gap = 8
        patch_height = max(24, min(64, rect.height() - 116))
        top = rect.top() + 90
        usable_width = max(1, rect.width() - margin * 2)
        patches: list[QRect] = []
        for index, value in enumerate(values):
            del value
            left = margin + round(index * usable_width / len(values))
            right = margin + round((index + 1) * usable_width / len(values))
            patches.append(
                QRect(
                    rect.left() + left,
                    top,
                    max(1, right - left - gap),
                    patch_height,
                )
            )
        return patches

    def _paint_section(
        self,
        painter: QPainter,
        rect: QRect,
        values: tuple[int, ...],
        *,
        selected: int,
        instruction: str,
        foreground: QColor,
        selection_color: QColor,
    ) -> None:
        background_value = 0 if values is CALIBRATION_BLACK_VALUES else 255
        painter.fillRect(
            rect,
            QColor(background_value, background_value, background_value),
        )
        painter.setPen(foreground)
        status_name = "BLACK" if values is CALIBRATION_BLACK_VALUES else "WHITE"
        changed = "  (just selected)" if self._last_selected == status_name.lower() else ""
        painter.drawText(
            QRect(rect.left() + 12, rect.top() + 6, rect.width() - 24, 46),
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop | Qt.TextFlag.TextWordWrap,
            f"{instruction}  Selected {status_name}={selected}{changed}",
        )
        patches = self._patch_rects(rect, values)
        for patch, value in zip(patches, values, strict=True):
            painter.fillRect(
                patch,
                QColor(value, value, value),
            )
            painter.setPen(foreground)
            painter.drawText(
                QRect(patch.left() - 1, patch.bottom() + 2, patch.width() + 2, 20),
                Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop,
                str(value),
            )
            if value == selected:
                painter.setPen(selection_color)
                painter.drawText(
                    QRect(patch.left(), patch.bottom() + 18, patch.width(), 18),
                    Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop,
                    "^",
                )

    def mousePressEvent(self, event: QMouseEvent) -> None:
        black_section = event.position().y() < self.height() / 2
        values = CALIBRATION_BLACK_VALUES if black_section else CALIBRATION_WHITE_VALUES
        top = 0 if black_section else self.height() // 2
        height = self.height() // 2 if black_section else self.height() - top
        rects = self._patch_rects(
            QRect(0, top, self.width(), height),
            values,
        )
        selected = next(
            (
                value
                for value, rect in zip(values, rects, strict=True)
                if rect.contains(event.position().toPoint())
            ),
            None,
        )
        if selected is None:
            return
        if black_section:
            self.black = selected
            self._last_selected = "black"
        else:
            self.white = selected
            self._last_selected = "white"
        self.selection_changed.emit(self.black, self.white)
        self.update()


class DisplayToneCalibrationDialog(QDialog):
    def __init__(
        self,
        curve: tuple[int, int] = (12, 247),
        *,
        show_copy: bool = False,
        fullscreen: bool = False,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Display tone curve calibration")
        if fullscreen:
            self.setWindowState(self.windowState() | Qt.WindowState.WindowFullScreen)
        else:
            self.resize(800, 600)
        layout = QVBoxLayout(self)
        self.pattern = DisplayToneCalibrationPattern(curve, self)
        self.pattern.selection_changed.connect(self._update_result)
        layout.addWidget(self.pattern, 1)
        row = QHBoxLayout()
        self.result_label = QLabel(self.option_fragment(), self)
        row.addWidget(self.result_label)
        row.addStretch(1)
        if show_copy:
            copy_button = QPushButton("Copy option", self)
            copy_button.clicked.connect(self.copy_option)
            row.addWidget(copy_button)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        row.addWidget(buttons)
        layout.addLayout(row)

    def curve(self) -> tuple[int, int]:
        return self.pattern.black, self.pattern.white

    def option_fragment(self) -> str:
        black, white = self.curve()
        return f"--display-tone-curve {black},{white}"

    def copy_option(self) -> None:
        from PySide6.QtWidgets import QApplication

        QApplication.clipboard().setText(self.option_fragment())

    def _update_result(self, _black: int, _white: int) -> None:
        self.result_label.setText(self.option_fragment())
