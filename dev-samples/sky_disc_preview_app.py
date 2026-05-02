#!/usr/bin/env python3
"""Standalone preview app for the sky-disc renderer.

This sample is intentionally kept outside the main zstarview UI so sky-disc
color behavior can be checked in isolation.

Usage:
    uv run dev-samples/sky_disc_preview_app.py

Keyboard:
    Esc  Quit
    R    Reset controls
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from PySide6.QtCore import Qt, Signal  # noqa: E402
from PySide6.QtGui import QColor, QFont, QKeyEvent, QPainter  # noqa: E402
from PySide6.QtWidgets import (  # noqa: E402
    QApplication,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QPushButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)

from zstarview.render.geometry import get_screen_geometry  # noqa: E402
from zstarview.render.sky_disc import draw_sky_color_disc  # noqa: E402


@dataclass(frozen=True)
class PreviewDefaults:
    sun_alt_deg: float = 7.5
    sun_az_deg: float = 90.0
    alpha: float = 1.0
    exposure: float = 1.14
    saturation: float = 1.35
    disc_opacity: float = 1.0
    view_center_alt_deg: float = 90.0
    view_center_az_deg: float = 0.0
    content_fov_deg: float = 90.0


class FloatSliderRow(QWidget):
    value_changed = Signal(float)

    def __init__(
        self,
        label: str,
        *,
        minimum: float,
        maximum: float,
        step: float,
        value: float,
        suffix: str = "",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._scale = 1.0 / step if step > 0.0 else 10.0
        self._suffix = suffix

        title = QLabel(label, self)
        title.setMinimumWidth(100)

        self._slider = QSlider(Qt.Orientation.Horizontal, self)
        self._slider.setRange(int(round(minimum * self._scale)), int(round(maximum * self._scale)))
        self._slider.valueChanged.connect(self._on_slider_value_changed)

        self._value_label = QLabel(self)
        self._value_label.setMinimumWidth(72)
        self._value_label.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.addWidget(title)
        layout.addWidget(self._slider, 1)
        layout.addWidget(self._value_label)

        self.set_value(value)

    def _format_value(self, raw_value: float) -> str:
        if math.isclose(raw_value, round(raw_value), abs_tol=1e-6):
            return f"{raw_value:.0f}{self._suffix}"
        return f"{raw_value:.2f}{self._suffix}"

    def _on_slider_value_changed(self, slider_value: int) -> None:
        value = float(slider_value) / self._scale
        self._value_label.setText(self._format_value(value))
        self.value_changed.emit(value)

    def value(self) -> float:
        return float(self._slider.value()) / self._scale

    def set_value(self, value: float) -> None:
        slider_value = int(round(float(value) * self._scale))
        self._slider.setValue(slider_value)


class SkyDiscPreviewWidget(QWidget):
    def __init__(self, defaults: PreviewDefaults, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumSize(420, 420)
        self.setAutoFillBackground(False)
        self._sun_alt_deg = float(defaults.sun_alt_deg)
        self._sun_az_deg = float(defaults.sun_az_deg)
        self._alpha = float(defaults.alpha)
        self._exposure = float(defaults.exposure)
        self._saturation = float(defaults.saturation)
        self._disc_opacity = float(defaults.disc_opacity)
        self._view_center_alt_deg = float(defaults.view_center_alt_deg)
        self._view_center_az_deg = float(defaults.view_center_az_deg)
        self._content_fov_deg = float(defaults.content_fov_deg)

    def set_sun_alt_deg(self, value: float) -> None:
        self._sun_alt_deg = float(value)
        self.update()

    def set_sun_az_deg(self, value: float) -> None:
        self._sun_az_deg = float(value)
        self.update()

    def set_alpha(self, value: float) -> None:
        self._alpha = float(value)
        self.update()

    def set_exposure(self, value: float) -> None:
        self._exposure = float(value)
        self.update()

    def set_saturation(self, value: float) -> None:
        self._saturation = float(value)
        self.update()

    def set_disc_opacity(self, value: float) -> None:
        self._disc_opacity = float(value)
        self.update()

    def set_view_center_alt_deg(self, value: float) -> None:
        self._view_center_alt_deg = float(value)
        self.update()

    def set_view_center_az_deg(self, value: float) -> None:
        self._view_center_az_deg = float(value)
        self.update()

    def set_content_fov_deg(self, value: float) -> None:
        self._content_fov_deg = float(value)
        self.update()

    def reset_to_defaults(self, defaults: PreviewDefaults) -> None:
        self._sun_alt_deg = float(defaults.sun_alt_deg)
        self._sun_az_deg = float(defaults.sun_az_deg)
        self._alpha = float(defaults.alpha)
        self._exposure = float(defaults.exposure)
        self._saturation = float(defaults.saturation)
        self._disc_opacity = float(defaults.disc_opacity)
        self._view_center_alt_deg = float(defaults.view_center_alt_deg)
        self._view_center_az_deg = float(defaults.view_center_az_deg)
        self._content_fov_deg = float(defaults.content_fov_deg)
        self.update()

    def paintEvent(self, _event) -> None:
        width = max(2, self.width())
        height = max(2, self.height())
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor(7, 9, 13))

        geometry = get_screen_geometry(width, height, self._view_center_alt_deg)
        image = draw_sky_color_disc(
            geometry,
            view_center=(self._view_center_alt_deg, self._view_center_az_deg),
            sun_altaz=(self._sun_alt_deg, self._sun_az_deg),
            alpha=self._alpha,
            exposure=self._exposure,
            saturation=self._saturation,
            disc_opacity=self._disc_opacity,
            eclipse_factor=1.0,
            content_fov_deg=self._content_fov_deg,
            image_size=(width, height),
        )
        painter.drawImage(0, 0, image)

        overlay = (
            f"sun_alt={self._sun_alt_deg:.1f} deg  "
            f"sun_az={self._sun_az_deg:.1f} deg  "
            f"alpha={self._alpha:.2f}  "
            f"exposure={self._exposure:.2f}  "
            f"saturation={self._saturation:.2f}"
        )
        painter.setPen(QColor(255, 255, 255, 210))
        font = QFont()
        font.setPointSize(10)
        painter.setFont(font)
        painter.drawText(14, 22, overlay)
        painter.end()


class SkyDiscPreviewWindow(QMainWindow):
    def __init__(self, defaults: PreviewDefaults) -> None:
        super().__init__()
        self._defaults = defaults
        self.setWindowTitle("Sky Disc Preview")
        self.resize(1100, 860)

        root = QWidget(self)
        self.setCentralWidget(root)

        self._preview = SkyDiscPreviewWidget(defaults, self)

        controls = QFrame(self)
        controls.setFrameShape(QFrame.Shape.StyledPanel)
        controls_layout = QVBoxLayout(controls)
        controls_layout.setContentsMargins(16, 16, 16, 16)
        controls_layout.setSpacing(12)

        title = QLabel("Sky Disc Controls", controls)
        title_font = QFont()
        title_font.setPointSize(14)
        title_font.setBold(True)
        title.setFont(title_font)
        controls_layout.addWidget(title)

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        form.setFormAlignment(Qt.AlignmentFlag.AlignTop)
        form.setHorizontalSpacing(10)
        form.setVerticalSpacing(10)

        self._sun_alt = FloatSliderRow("Sun altitude", minimum=-12.0, maximum=60.0, step=0.1, value=defaults.sun_alt_deg, suffix=" deg")
        self._sun_az = FloatSliderRow("Sun azimuth", minimum=0.0, maximum=360.0, step=0.1, value=defaults.sun_az_deg, suffix=" deg")
        self._view_alt = FloatSliderRow("View center alt", minimum=0.0, maximum=90.0, step=0.1, value=defaults.view_center_alt_deg, suffix=" deg")
        self._view_az = FloatSliderRow("View center az", minimum=0.0, maximum=360.0, step=0.1, value=defaults.view_center_az_deg, suffix=" deg")
        self._content_fov = FloatSliderRow("Content FOV", minimum=40.0, maximum=120.0, step=0.1, value=defaults.content_fov_deg, suffix=" deg")
        self._alpha = FloatSliderRow("Alpha", minimum=0.0, maximum=1.0, step=0.01, value=defaults.alpha)
        self._exposure = FloatSliderRow("Exposure", minimum=0.6, maximum=1.6, step=0.01, value=defaults.exposure)
        self._saturation = FloatSliderRow("Saturation", minimum=0.6, maximum=1.8, step=0.01, value=defaults.saturation)
        self._disc_opacity = FloatSliderRow("Disc opacity", minimum=0.0, maximum=1.0, step=0.01, value=defaults.disc_opacity)

        form.addRow(self._sun_alt)
        form.addRow(self._sun_az)
        form.addRow(self._view_alt)
        form.addRow(self._view_az)
        form.addRow(self._content_fov)
        form.addRow(self._alpha)
        form.addRow(self._exposure)
        form.addRow(self._saturation)
        form.addRow(self._disc_opacity)
        controls_layout.addLayout(form)

        reset_button = QPushButton("Reset controls", controls)
        reset_button.clicked.connect(self._reset_controls)
        controls_layout.addWidget(reset_button)
        controls_layout.addStretch(1)

        layout = QHBoxLayout(root)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)
        layout.addWidget(self._preview, 1)
        layout.addWidget(controls, 0)

        self._sun_alt.value_changed.connect(self._preview.set_sun_alt_deg)
        self._sun_az.value_changed.connect(self._preview.set_sun_az_deg)
        self._view_alt.value_changed.connect(self._preview.set_view_center_alt_deg)
        self._view_az.value_changed.connect(self._preview.set_view_center_az_deg)
        self._content_fov.value_changed.connect(self._preview.set_content_fov_deg)
        self._alpha.value_changed.connect(self._preview.set_alpha)
        self._exposure.value_changed.connect(self._preview.set_exposure)
        self._saturation.value_changed.connect(self._preview.set_saturation)
        self._disc_opacity.value_changed.connect(self._preview.set_disc_opacity)

    def _reset_controls(self) -> None:
        self._sun_alt.set_value(self._defaults.sun_alt_deg)
        self._sun_az.set_value(self._defaults.sun_az_deg)
        self._view_alt.set_value(self._defaults.view_center_alt_deg)
        self._view_az.set_value(self._defaults.view_center_az_deg)
        self._content_fov.set_value(self._defaults.content_fov_deg)
        self._alpha.set_value(self._defaults.alpha)
        self._exposure.set_value(self._defaults.exposure)
        self._saturation.set_value(self._defaults.saturation)
        self._disc_opacity.set_value(self._defaults.disc_opacity)
        self._preview.reset_to_defaults(self._defaults)

    def keyPressEvent(self, event: QKeyEvent) -> None:
        if event.key() in (Qt.Key.Key_Escape, Qt.Key.Key_Q):
            self.close()
            event.accept()
            return
        if event.key() == Qt.Key.Key_R:
            self._reset_controls()
            event.accept()
            return
        super().keyPressEvent(event)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Preview the sky-disc renderer in isolation.")
    parser.add_argument("--sun-alt", type=float, default=7.5, help="Sun altitude in degrees.")
    parser.add_argument("--sun-az", type=float, default=90.0, help="Sun azimuth in degrees.")
    parser.add_argument("--alpha", type=float, default=1.0, help="Sky-disc alpha used by the color model.")
    parser.add_argument("--exposure", type=float, default=1.14, help="Sky-disc exposure.")
    parser.add_argument("--saturation", type=float, default=1.35, help="Sky-disc saturation.")
    parser.add_argument("--disc-opacity", type=float, default=1.0, help="Final disc opacity.")
    parser.add_argument("--view-alt", type=float, default=90.0, help="View-center altitude in degrees.")
    parser.add_argument("--view-az", type=float, default=0.0, help="View-center azimuth in degrees.")
    parser.add_argument("--content-fov", type=float, default=90.0, help="Content field of view in degrees.")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    app = QApplication(sys.argv)
    window = SkyDiscPreviewWindow(
        PreviewDefaults(
            sun_alt_deg=float(args.sun_alt),
            sun_az_deg=float(args.sun_az),
            alpha=float(args.alpha),
            exposure=float(args.exposure),
            saturation=float(args.saturation),
            disc_opacity=float(args.disc_opacity),
            view_center_alt_deg=float(args.view_alt),
            view_center_az_deg=float(args.view_az),
            content_fov_deg=float(args.content_fov),
        )
    )
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
