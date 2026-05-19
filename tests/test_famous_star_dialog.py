from __future__ import annotations

import polars as pl
from PySide6.QtWidgets import QApplication

from zstarview.gui.famous_star_dialog import NamedStarJumpDialog
from zstarview.gui.famous_star_shortcuts import build_named_star_shortcuts

_app = QApplication.instance() or QApplication([])


def test_named_star_jump_dialog_includes_satellite_tab_and_correct_title() -> None:
    star_catalog = pl.DataFrame(
        {
            "Name": ["Sirius", "Vega", "Canopus"],
            "RAh": [6.75, 18.6, 6.4],
            "Dec": [-16.7, 38.7, -52.7],
            "Vmag": [-1.44, 0.03, -0.62],
        }
    )
    dialog = NamedStarJumpDialog(build_named_star_shortcuts(star_catalog, max_vmag=2.0, include_satellites=False))

    assert dialog.windowTitle() == "Jump to Named Star"
    assert dialog._tabs.count() == 4
    assert dialog._tabs.tabText(0).startswith("North (")
    assert dialog._tabs.tabText(1).startswith("Equatorial (")
    assert dialog._tabs.tabText(2).startswith("South (")
    assert dialog._tabs.tabText(3).startswith("Artificial Satellites (")

    satellites_widget = dialog._tabs.widget(3)
    satellites_widget.setCurrentRow(0)
    assert dialog.selected_star() is not None
    assert dialog.selected_star().name == "ISS"
