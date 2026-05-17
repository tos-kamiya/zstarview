from __future__ import annotations

from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QWidget,
)

from ..paths import OBSERVER_MAX_ALT_DEG, OBSERVER_MIN_ALT_DEG
from .view_direction import clamp_view_center_alt_az


class ViewDirectionDialog(QDialog):
    def __init__(
        self,
        view_center: tuple[float, float],
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Set View Center")
        self.setModal(True)
        self.resize(320, 140)

        layout = QFormLayout(self)

        alt_deg, az_deg = clamp_view_center_alt_az(view_center[0], view_center[1])

        self._alt_spin = QDoubleSpinBox(self)
        self._alt_spin.setRange(OBSERVER_MIN_ALT_DEG, OBSERVER_MAX_ALT_DEG)
        self._alt_spin.setDecimals(2)
        self._alt_spin.setSingleStep(1.0)
        self._alt_spin.setValue(alt_deg)
        self._alt_spin.setSuffix(" deg")
        layout.addRow("Altitude", self._alt_spin)

        self._az_spin = QDoubleSpinBox(self)
        self._az_spin.setRange(0.0, 360.0)
        self._az_spin.setDecimals(2)
        self._az_spin.setSingleStep(1.0)
        self._az_spin.setValue(az_deg)
        self._az_spin.setSuffix(" deg")
        layout.addRow("Azimuth", self._az_spin)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

    def selected_view_center(self) -> tuple[float, float]:
        return (
            float(self._alt_spin.value()),
            float(self._az_spin.value()),
        )
