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

from .famous_star_shortcuts import SearchJumpTarget


class NamedStarSearchDialog(QDialog):
    def __init__(self, targets: List[SearchJumpTarget], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Search Stars and Asterisms")
        self.setModal(True)
        self.resize(460, 460)

        layout = QVBoxLayout(self)
        self._search = QLineEdit(self)
        self._search.setPlaceholderText("Type star or asterism name...")
        self._search.textChanged.connect(self._apply_filter)
        layout.addWidget(self._search)

        self._list = QListWidget(self)
        self._list.itemDoubleClicked.connect(self._on_item_double_clicked)
        for target in targets:
            suffix = f"  ({target.subtitle})" if target.subtitle else ""
            item = QListWidgetItem(f"{target.label}{suffix}", self._list)
            item.setData(Qt.ItemDataRole.UserRole, target)
        layout.addWidget(self._list)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel, parent=self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self._select_first_visible()

    def selected_target(self) -> Optional[SearchJumpTarget]:
        item = self._list.currentItem()
        if item is None:
            return None
        target = item.data(Qt.ItemDataRole.UserRole)
        return target if isinstance(target, SearchJumpTarget) else None

    def _apply_filter(self, text: str) -> None:
        query = text.strip().casefold()
        for i in range(self._list.count()):
            item = self._list.item(i)
            target = item.data(Qt.ItemDataRole.UserRole)
            if not isinstance(target, SearchJumpTarget):
                item.setHidden(True)
                continue
            if not query:
                item.setHidden(False)
                continue
            haystack = f"{target.label} {target.subtitle}".casefold()
            item.setHidden(query not in haystack)
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
