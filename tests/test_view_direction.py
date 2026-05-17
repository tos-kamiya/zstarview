from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from zstarview.gui.view_direction import clamp_view_center_alt_az, resolve_view_direction_step
from zstarview.gui.view_direction_dialog import ViewDirectionDialog

app = QApplication.instance() or QApplication([])


def test_resolve_view_direction_step_uses_fine_step_with_shift() -> None:
    assert resolve_view_direction_step(Qt.KeyboardModifier.NoModifier, 5.0) == 5.0
    assert (
        resolve_view_direction_step(Qt.KeyboardModifier.ShiftModifier, 5.0)
        == 1.0
    )


def test_clamp_view_center_alt_az_clamps_altitude_and_wraps_azimuth() -> None:
    alt_deg, az_deg = clamp_view_center_alt_az(-60.0, 725.0)

    assert alt_deg == -45.0
    assert az_deg == 5.0


def test_view_direction_dialog_initializes_from_existing_center() -> None:
    dialog = ViewDirectionDialog((-60.0, 725.0))

    assert dialog.selected_view_center() == (-45.0, 5.0)
