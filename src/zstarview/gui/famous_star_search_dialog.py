# -*- coding: utf-8 -*-
"""Search-first dialog for named star and small-body jump targets."""
from __future__ import annotations

import threading
from dataclasses import replace
from typing import Callable, List, Optional, Sequence

from PySide6.QtCore import QEvent, Qt, Signal
from PySide6.QtGui import QKeyEvent
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .famous_star_shortcuts import SearchJumpTarget


class NamedStarSearchDialog(QDialog):
    jpl_search_finished = Signal(int, object, str)

    def __init__(
        self,
        targets: List[SearchJumpTarget],
        parent: QWidget | None = None,
        jpl_search_callback: Callable[[str], Sequence[SearchJumpTarget]] | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Search Stars and Asterisms")
        self.setModal(True)
        self.resize(560, 560)

        self._jpl_search_callback = jpl_search_callback
        self._jpl_targets: list[SearchJumpTarget] = []
        self._jpl_search_request_id = 0
        self.jpl_search_finished.connect(self._on_jpl_search_finished)

        outer = QVBoxLayout(self)
        self._tabs = QTabWidget(self)
        self._tabs.currentChanged.connect(lambda *_: self._sync_ok_button())
        outer.addWidget(self._tabs)

        self._local_tab = QWidget(self)
        local_layout = QVBoxLayout(self._local_tab)
        self._search = QLineEdit(self._local_tab)
        self._search.setPlaceholderText("Type star or asterism name...")
        self._search.textChanged.connect(self._apply_local_filter)
        local_layout.addWidget(self._search)

        self._list = QListWidget(self._local_tab)
        self._list.itemDoubleClicked.connect(self._on_item_double_clicked)
        for target in targets:
            suffix = f"  ({target.subtitle})" if target.subtitle else ""
            item = QListWidgetItem(f"{target.label}{suffix}", self._list)
            item.setData(Qt.ItemDataRole.UserRole, target)
        local_layout.addWidget(self._list)
        self._tabs.addTab(self._local_tab, "Stars and Asterisms")

        self._jpl_tab = QWidget(self)
        jpl_layout = QVBoxLayout(self._jpl_tab)
        jpl_layout.addWidget(
            QLabel("Search JPL small bodies and keep the selected target on the map.", self._jpl_tab)
        )
        self._jpl_search = QLineEdit(self._jpl_tab)
        self._jpl_search.setPlaceholderText("Type Ceres, Eris, Haumea, Makemake...")
        self._jpl_search.returnPressed.connect(self._start_jpl_search)
        self._jpl_search.installEventFilter(self)
        jpl_layout.addWidget(self._jpl_search)

        self._jpl_search_button = QPushButton("Search", self._jpl_tab)
        self._jpl_search_button.clicked.connect(self._start_jpl_search)
        jpl_layout.addWidget(self._jpl_search_button)

        self._jpl_status = QLabel("", self._jpl_tab)
        self._jpl_status.setVisible(False)
        jpl_layout.addWidget(self._jpl_status)

        self._jpl_list = QListWidget(self._jpl_tab)
        self._jpl_list.itemDoubleClicked.connect(self._on_item_double_clicked)
        jpl_layout.addWidget(self._jpl_list)

        keep_row = QHBoxLayout()
        self._jpl_keep_marker = QCheckBox("Keep marker", self._jpl_tab)
        self._jpl_keep_label = QCheckBox("Keep label", self._jpl_tab)
        keep_row.addWidget(self._jpl_keep_marker)
        keep_row.addWidget(self._jpl_keep_label)
        keep_row.addStretch(1)
        jpl_layout.addLayout(keep_row)

        self._tabs.addTab(self._jpl_tab, "JPL Small Bodies")

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        self._ok_button = buttons.button(QDialogButtonBox.StandardButton.Ok)
        if self._ok_button is not None:
            self._ok_button.setEnabled(False)
        outer.addWidget(buttons)

        self._select_first_visible_local()

    def selected_target(self) -> Optional[SearchJumpTarget]:
        current_tab = self._tabs.currentIndex()
        if current_tab == 1:
            item = self._jpl_list.currentItem()
            if item is None:
                return None
            target = item.data(Qt.ItemDataRole.UserRole)
            if not isinstance(target, SearchJumpTarget):
                return None
            return replace(
                target,
                persistent_keep_marker=self._jpl_keep_marker.isChecked(),
                persistent_keep_label=self._jpl_keep_label.isChecked(),
            )

        item = self._list.currentItem()
        if item is None:
            return None
        target = item.data(Qt.ItemDataRole.UserRole)
        return target if isinstance(target, SearchJumpTarget) else None

    def _apply_local_filter(self, text: str) -> None:
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
        self._select_first_visible_local()

    def _select_first_visible_local(self) -> None:
        for i in range(self._list.count()):
            item = self._list.item(i)
            if not item.isHidden():
                self._list.setCurrentItem(item)
                self._sync_ok_button()
                return
        self._list.setCurrentItem(None)
        self._sync_ok_button()

    def _select_first_visible_jpl(self) -> None:
        for i in range(self._jpl_list.count()):
            item = self._jpl_list.item(i)
            if not item.isHidden():
                self._jpl_list.setCurrentItem(item)
                self._sync_ok_button()
                return
        self._jpl_list.setCurrentItem(None)
        self._sync_ok_button()

    def _sync_ok_button(self) -> None:
        if self._ok_button is not None:
            self._ok_button.setEnabled(self.selected_target() is not None)

    def _on_item_double_clicked(self, _item: QListWidgetItem) -> None:
        self.accept()

    def _start_jpl_search(self) -> None:
        query = self._jpl_search.text().strip()
        if not query:
            return
        if self._jpl_search_callback is None:
            self._jpl_status.setText("JPL search is unavailable.")
            self._jpl_status.setVisible(True)
            return
        self._jpl_search_request_id += 1
        request_id = self._jpl_search_request_id
        self._jpl_search_button.setEnabled(False)
        self._jpl_status.setText(f"Searching JPL for '{query}'...")
        self._jpl_status.setVisible(True)
        worker = threading.Thread(target=self._run_jpl_search, args=(request_id, query), daemon=True)
        worker.start()

    def _run_jpl_search(self, request_id: int, query: str) -> None:
        try:
            targets = list(self._jpl_search_callback(query)) if self._jpl_search_callback is not None else []
        except Exception as exc:  # pragma: no cover - exercised through signal delivery
            self.jpl_search_finished.emit(request_id, [], f"JPL search failed: {exc}")
            return
        status_text = f"Found {len(targets)} JPL result(s)" if targets else "No JPL results found"
        self.jpl_search_finished.emit(request_id, targets, status_text)

    def _on_jpl_search_finished(self, request_id: int, payload: object, status_text: str) -> None:
        if request_id != self._jpl_search_request_id:
            return
        self._jpl_search_button.setEnabled(True)
        self._jpl_targets = [target for target in payload if isinstance(target, SearchJumpTarget)]
        self._jpl_list.clear()
        for target in self._jpl_targets:
            suffix = f"  ({target.subtitle})" if target.subtitle else ""
            item = QListWidgetItem(f"{target.label}{suffix}", self._jpl_list)
            item.setData(Qt.ItemDataRole.UserRole, target)
        self._jpl_status.setText(status_text)
        self._jpl_status.setVisible(bool(status_text))
        self._select_first_visible_jpl()

    def eventFilter(self, watched: object, event: object) -> bool:
        if watched is self._jpl_search and isinstance(event, QKeyEvent) and event.type() == QEvent.Type.KeyPress:
            if event.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
                self._start_jpl_search()
                return True
        return super().eventFilter(watched, event)
