# -*- coding: utf-8 -*-
"""Search-first dialog for named star jump."""
from __future__ import annotations

from typing import List, Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .famous_star_shortcuts import NamedStarShortcut


class NamedStarSearchDialog(QDialog):
    def __init__(self, stars: List[NamedStarShortcut], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Search Named Stars")
        self.setModal(True)
        self.resize(460, 460)

        layout = QVBoxLayout(self)
        self._search = QLineEdit(self)
        self._search.setPlaceholderText("Type star name...")
        self._search.textChanged.connect(self._apply_filter)
        layout.addWidget(self._search)

        self._list = QListWidget(self)
        self._list.itemDoubleClicked.connect(self._on_item_double_clicked)
        for star in stars:
            item = QListWidgetItem(f"{star.name}  (Vmag {star.vmag:.2f})", self._list)
            item.setData(Qt.ItemDataRole.UserRole, star)
        layout.addWidget(self._list)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel, parent=self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self._select_first_visible()

    def selected_star(self) -> Optional[NamedStarShortcut]:
        item = self._list.currentItem()
        if item is None:
            return None
        star = item.data(Qt.ItemDataRole.UserRole)
        return star if isinstance(star, NamedStarShortcut) else None

    def _apply_filter(self, text: str) -> None:
        query = text.strip().casefold()
        for i in range(self._list.count()):
            item = self._list.item(i)
            star = item.data(Qt.ItemDataRole.UserRole)
            if not isinstance(star, NamedStarShortcut):
                item.setHidden(True)
                continue
            if not query:
                item.setHidden(False)
                continue
            item.setHidden(query not in star.name.casefold())
        self._select_first_visible()

    def _select_first_visible(self) -> None:
        for i in range(self._list.count()):
            item = self._list.item(i)
            if not item.isHidden():
                self._list.setCurrentItem(item)
                return
        self._list.setCurrentItem(None)

    def _on_item_double_clicked(self, _item: QListWidgetItem) -> None:
        self.accept()
