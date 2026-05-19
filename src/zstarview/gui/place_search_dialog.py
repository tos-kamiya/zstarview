# -*- coding: utf-8 -*-
"""Dialog for online place search targets."""
from __future__ import annotations

from dataclasses import replace
from typing import Callable, Optional, Sequence

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
    QVBoxLayout,
    QWidget,
)

from ..search.models import SearchJumpTarget
from .worker_pool import submit_gui_work


class PlaceSearchDialog(QDialog):
    place_search_finished = Signal(int, object, str)

    def __init__(
        self,
        place_search_callback: Callable[[str], Sequence[SearchJumpTarget]],
        parent: QWidget | None = None,
        *,
        cli_view_center_alt_specified: bool = False,
        cli_view_center_az_specified: bool = False,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Search Places")
        self.setModal(True)
        self.resize(460, 460)

        self._place_search_callback = place_search_callback
        self._place_targets: list[SearchJumpTarget] = []
        self._place_search_request_id = 0
        self._cli_view_center_alt_specified = bool(cli_view_center_alt_specified)
        self._cli_view_center_az_specified = bool(cli_view_center_az_specified)
        self.place_search_finished.connect(self._on_place_search_finished)

        layout = QVBoxLayout(self)
        self._search = QLineEdit(self)
        self._search.setPlaceholderText("Type place, station, or facility name...")
        self._search.returnPressed.connect(self._start_place_search)
        self._search.installEventFilter(self)
        layout.addWidget(self._search)

        self._place_search_button = QPushButton("Search", self)
        self._place_search_button.clicked.connect(self._start_place_search)
        layout.addWidget(self._place_search_button)

        self._status = QLabel("", self)
        self._status.setVisible(False)
        layout.addWidget(self._status)

        cli_keep_row = QHBoxLayout()
        self._cli_view_center_keep = QCheckBox("Keep CLI-specified Alt/Az", self)
        cli_view_center_enabled = (
            self._cli_view_center_alt_specified or self._cli_view_center_az_specified
        )
        self._cli_view_center_keep.setEnabled(cli_view_center_enabled)
        self._cli_view_center_keep.setChecked(cli_view_center_enabled)
        cli_keep_row.addWidget(self._cli_view_center_keep)
        cli_keep_row.addStretch(1)
        layout.addLayout(cli_keep_row)

        self._list = QListWidget(self)
        self._list.itemDoubleClicked.connect(self._on_item_double_clicked)
        layout.addWidget(self._list)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel, parent=self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        self._ok_button = buttons.button(QDialogButtonBox.StandardButton.Ok)
        if self._ok_button is not None:
            self._ok_button.setEnabled(False)
        layout.addWidget(buttons)

    def selected_target(self) -> Optional[SearchJumpTarget]:
        item = self._list.currentItem()
        if item is None:
            return None
        target = item.data(Qt.ItemDataRole.UserRole)
        if not isinstance(target, SearchJumpTarget):
            return None
        return replace(
            target,
            preserve_cli_view_center=self._cli_view_center_keep.isChecked(),
        )

    def _rebuild_list(self) -> None:
        self._list.clear()
        for target in self._place_targets:
            suffix = f"  ({target.subtitle})" if target.subtitle else ""
            item = QListWidgetItem(f"{target.label}{suffix}", self._list)
            item.setData(Qt.ItemDataRole.UserRole, target)
        self._select_first_visible()
        self._sync_ok_button()

    def _select_first_visible(self) -> None:
        for i in range(self._list.count()):
            item = self._list.item(i)
            if not item.isHidden():
                self._list.setCurrentItem(item)
                return
        self._list.setCurrentItem(None)

    def _sync_ok_button(self) -> None:
        if self._ok_button is not None:
            self._ok_button.setEnabled(self.selected_target() is not None)

    def _on_item_double_clicked(self, _item: QListWidgetItem) -> None:
        self.accept()

    def _start_place_search(self) -> None:
        query = self._search.text().strip()
        if not query:
            return
        self._place_search_request_id += 1
        request_id = self._place_search_request_id
        self._place_search_button.setEnabled(False)
        self._status.setText(f"Searching places for '{query}'...")
        self._status.setVisible(True)
        submit_gui_work(self._run_place_search, request_id=request_id, query=query)

    def _run_place_search(self, request_id: int, query: str) -> None:
        try:
            targets = list(self._place_search_callback(query))
        except Exception as exc:  # pragma: no cover - exercised through signal delivery
            self.place_search_finished.emit(request_id, [], f"Place search failed: {exc}")
            return
        status_text = f"Found {len(targets)} place result(s)" if targets else "No places found"
        self.place_search_finished.emit(request_id, targets, status_text)

    def _on_place_search_finished(self, request_id: int, payload: object, status_text: str) -> None:
        if request_id != self._place_search_request_id:
            return
        self._place_search_button.setEnabled(True)
        self._place_targets = [target for target in payload if isinstance(target, SearchJumpTarget)]
        self._status.setText(status_text)
        self._status.setVisible(bool(status_text))
        self._rebuild_list()

    def eventFilter(self, watched: object, event: object) -> bool:
        if watched is self._search and isinstance(event, QKeyEvent) and event.type() == QEvent.Type.KeyPress:
            if event.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
                self._start_place_search()
                return True
        return super().eventFilter(watched, event)
