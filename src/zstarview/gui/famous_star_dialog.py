# -*- coding: utf-8 -*-
"""Dialog for selecting a named star to jump to."""
from __future__ import annotations

from typing import Dict, List, Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QListWidget,
    QListWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .famous_star_shortcuts import (
    DEC_BAND_EQUATOR,
    DEC_BAND_NORTH,
    DEC_BAND_SOUTH,
    NamedStarShortcut,
)


class NamedStarJumpDialog(QDialog):
    def __init__(self, stars_by_band: Dict[str, List[NamedStarShortcut]], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Jump to Famous Star")
        self.setModal(True)
        self.resize(460, 420)

        layout = QVBoxLayout(self)
        self._tabs = QTabWidget(self)
        self._lists: Dict[str, QListWidget] = {}

        tab_defs = [
            (DEC_BAND_NORTH, "North"),
            (DEC_BAND_EQUATOR, "Equatorial"),
            (DEC_BAND_SOUTH, "South"),
        ]
        for key, label in tab_defs:
            lw = QListWidget(self)
            lw.itemDoubleClicked.connect(self._on_item_double_clicked)
            self._lists[key] = lw
            self._tabs.addTab(lw, f"{label} ({len(stars_by_band.get(key, []))})")
            for star in stars_by_band.get(key, []):
                suffix = star.subtitle or (f"Vmag {star.vmag:.2f}" if star.kind == "star" else "")
                item = QListWidgetItem(f"{star.name}  ({suffix})" if suffix else star.name, lw)
                item.setData(Qt.ItemDataRole.UserRole, star)

        layout.addWidget(self._tabs)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel, parent=self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def selected_star(self) -> Optional[NamedStarShortcut]:
        widget = self._tabs.currentWidget()
        if not isinstance(widget, QListWidget):
            return None
        item = widget.currentItem()
        if item is None:
            return None
        star = item.data(Qt.ItemDataRole.UserRole)
        return star if isinstance(star, NamedStarShortcut) else None

    def _on_item_double_clicked(self, _item: QListWidgetItem) -> None:
        self.accept()
