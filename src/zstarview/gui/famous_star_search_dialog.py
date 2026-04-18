# -*- coding: utf-8 -*-
"""Search-first dialog for named star, asterism, and small-body targets."""
from __future__ import annotations

import threading
from dataclasses import replace
from typing import Callable, List, Optional, Sequence

from PySide6.QtCore import Qt, Signal
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

        self._local_targets = list(targets)
        self._jpl_search_callback = jpl_search_callback
        self._jpl_search_request_id = 0
        self._jpl_search_in_progress = False
        self._local_result_count = 0
        self.jpl_search_finished.connect(self._on_jpl_search_finished)

        outer = QVBoxLayout(self)

        intro = QLabel(
            "Search stars and asterisms first. If there is no local match, the dialog can fall back to JPL small bodies.",
            self,
        )
        intro.setWordWrap(True)
        outer.addWidget(intro)

        search_row = QHBoxLayout()
        self._search = QLineEdit(self)
        self._search.setPlaceholderText("Type star, asterism, or small-body name...")
        self._search.textChanged.connect(self._apply_local_filter)
        self._search.returnPressed.connect(self.accept)
        search_row.addWidget(self._search)
        outer.addLayout(search_row)

        self._jpl_search_button = QPushButton("Search JPL database", self)
        self._jpl_search_button.clicked.connect(self._start_jpl_search)
        self._jpl_search_button.setEnabled(False)
        outer.addWidget(self._jpl_search_button)

        self._list = QListWidget(self)
        self._list.itemDoubleClicked.connect(self._on_item_double_clicked)
        outer.addWidget(self._list)

        self._status = QLabel("", self)
        self._status.setVisible(False)
        outer.addWidget(self._status)

        keep_row = QHBoxLayout()
        self._jpl_keep_marker = QCheckBox("Keep marker", self)
        self._jpl_keep_label = QCheckBox("Keep label", self)
        keep_row.addWidget(self._jpl_keep_marker)
        keep_row.addWidget(self._jpl_keep_label)
        keep_row.addStretch(1)
        outer.addLayout(keep_row)

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

        self._apply_local_filter("")

    def selected_target(self) -> Optional[SearchJumpTarget]:
        item = self._list.currentItem()
        if item is None:
            return None
        target = item.data(Qt.ItemDataRole.UserRole)
        if not isinstance(target, SearchJumpTarget):
            return None
        if target.kind != "jpl_small_body":
            return target
        return replace(
            target,
            persistent_keep_marker=self._jpl_keep_marker.isChecked(),
            persistent_keep_label=self._jpl_keep_label.isChecked(),
        )

    def accept(self) -> None:  # noqa: D401 - Qt override
        target = self.selected_target()
        if target is not None:
            super().accept()
            return

        query = self._search.text().strip()
        if not query:
            return
        if self._local_result_count > 0:
            self._select_first_visible()
            if self.selected_target() is not None:
                super().accept()
            return
        self._start_jpl_search()

    def _apply_local_filter(self, text: str) -> None:
        query = text.strip().casefold()
        self._list.clear()

        matching_targets: list[SearchJumpTarget] = []
        for target in self._local_targets:
            haystack = f"{target.label} {target.subtitle}".casefold()
            if not query or query in haystack:
                matching_targets.append(target)

        self._local_result_count = len(matching_targets)
        for target in matching_targets:
            suffix = f"  ({target.subtitle})" if target.subtitle else ""
            item = QListWidgetItem(f"{target.label}{suffix}", self._list)
            item.setData(Qt.ItemDataRole.UserRole, target)

        if self._local_result_count > 0:
            self._set_status(
                f"Found {self._local_result_count} local result(s). Select one and press OK."
            )
            self._select_first_visible()
        elif query:
            self._set_status("No local match. Use the JPL database button below.")
        else:
            self._set_status("")
        self._sync_ok_button()
        self._sync_jpl_button()

    def _select_first_visible(self) -> None:
        for i in range(self._list.count()):
            item = self._list.item(i)
            if not item.isHidden():
                self._list.setCurrentItem(item)
                self._sync_ok_button()
                return
        self._list.setCurrentItem(None)
        self._sync_ok_button()

    def _set_status(self, text: str) -> None:
        self._status.setText(text)
        self._status.setVisible(bool(text))

    def _sync_jpl_button(self) -> None:
        if not hasattr(self, "_jpl_search_button"):
            return
        query = self._search.text().strip()
        enabled = bool(query) and self._local_result_count == 0 and not self._jpl_search_in_progress
        self._jpl_search_button.setEnabled(enabled)
        if not query:
            self._jpl_search_button.setText("Search JPL database")
        else:
            self._jpl_search_button.setText(f"Search JPL database for '{query}'")

    def _sync_ok_button(self) -> None:
        if self._ok_button is not None:
            self._ok_button.setEnabled(self.selected_target() is not None)

    def _on_item_double_clicked(self, _item: QListWidgetItem) -> None:
        self.accept()

    def _start_jpl_search(self) -> None:
        query = self._search.text().strip()
        if not query or self._jpl_search_in_progress:
            return
        if self._jpl_search_callback is None:
            self._set_status("JPL search is unavailable.")
            self._sync_ok_button()
            self._sync_jpl_button()
            return
        self._jpl_search_request_id += 1
        request_id = self._jpl_search_request_id
        self._jpl_search_in_progress = True
        if self._ok_button is not None:
            self._ok_button.setEnabled(False)
        self._sync_jpl_button()
        self._set_status(f"Searching JPL for '{query}'...")
        worker = threading.Thread(
            target=self._run_jpl_search,
            args=(request_id, query),
            daemon=True,
        )
        worker.start()

    def _run_jpl_search(self, request_id: int, query: str) -> None:
        try:
            targets = (
                list(self._jpl_search_callback(query))
                if self._jpl_search_callback is not None
                else []
            )
        except Exception as exc:  # pragma: no cover - signal delivery path
            self.jpl_search_finished.emit(request_id, [], f"JPL search failed: {exc}")
            return
        status_text = (
            f"Found {len(targets)} JPL result(s)"
            if targets
            else "No JPL results found"
        )
        self.jpl_search_finished.emit(request_id, targets, status_text)

    def _on_jpl_search_finished(
        self, request_id: int, payload: object, status_text: str
    ) -> None:
        if request_id != self._jpl_search_request_id:
            return
        self._jpl_search_in_progress = False
        targets = [target for target in payload if isinstance(target, SearchJumpTarget)]
        self._list.clear()
        for target in targets:
            suffix = f"  ({target.subtitle})" if target.subtitle else ""
            item = QListWidgetItem(f"{target.label}{suffix}", self._list)
            item.setData(Qt.ItemDataRole.UserRole, target)
        self._set_status(status_text)
        if targets:
            self._list.setCurrentRow(0)
        else:
            self._list.setCurrentItem(None)
        self._sync_ok_button()
        self._sync_jpl_button()
