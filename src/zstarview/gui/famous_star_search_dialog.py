"""Search-first dialog for named star, asterism, satellite, and JPL targets."""
from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import replace

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
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from ..search.constants import SOLAR_SYSTEM_BODY_QUERIES
from ..search.models import SearchJumpTarget
from ..search.query import parse_search_query, search_target_matches_query
from .application_services import ApplicationServices

_JPL_BYPASS_QUERIES = SOLAR_SYSTEM_BODY_QUERIES


class NamedStarSearchDialog(QDialog):
    jpl_search_finished = Signal(int, object, str)

    def __init__(
        self,
        targets: list[SearchJumpTarget],
        parent: QWidget | None = None,
        *,
        services: ApplicationServices | None = None,
        cli_view_center_alt_specified: bool = False,
        cli_view_center_az_specified: bool = False,
        satellite_search_callback: Callable[[str], Sequence[SearchJumpTarget]] | None = None,
        jpl_search_callback: Callable[[str], Sequence[SearchJumpTarget]] | None = None,
    ) -> None:
        super().__init__(parent)
        self._services = services or ApplicationServices()
        self.setWindowTitle("Search Objects")
        self.setModal(True)
        self.resize(560, 560)

        self._local_targets = list(targets)
        self._satellite_search_callback = satellite_search_callback
        self._jpl_search_callback = jpl_search_callback
        self._jpl_search_request_id = 0
        self._jpl_search_in_progress = False
        self._local_result_count = 0
        self._clear_persistent_marker_on_accept = False
        self._cli_view_center_alt_specified = bool(cli_view_center_alt_specified)
        self._cli_view_center_az_specified = bool(cli_view_center_az_specified)
        self.jpl_search_finished.connect(self._on_jpl_search_finished)

        outer = QVBoxLayout(self)

        intro = QLabel(
            "Search stars and asterisms first. If there is no local match, the dialog can resolve known artificial satellites and then fall back to JPL bodies. JPL search shows up to 500 results.",
            self,
        )
        intro.setWordWrap(True)
        outer.addWidget(intro)

        search_row = QHBoxLayout()
        self._search = QLineEdit(self)
        self._search.setPlaceholderText(
            "Type star, asterism, satellite, or JPL body name..."
        )
        self._search.textChanged.connect(self._apply_local_filter)
        self._search.returnPressed.connect(self.accept)
        search_row.addWidget(self._search)
        outer.addLayout(search_row)

        self._jpl_search_button = QPushButton("Search satellites / JPL (up to 500)", self)
        self._jpl_search_button.clicked.connect(self._start_jpl_search)
        self._jpl_search_button.setEnabled(False)
        outer.addWidget(self._jpl_search_button)

        self._list = QListWidget(self)
        self._list.itemDoubleClicked.connect(self._on_item_double_clicked)
        self._list.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        self._list.setMinimumHeight(220)
        outer.addWidget(self._list)

        self._status = QLabel("", self)
        self._status.setVisible(False)
        outer.addWidget(self._status)

        keep_row = QHBoxLayout()
        self._jpl_keep_marker = QCheckBox("Keep marker", self)
        keep_row.addWidget(self._jpl_keep_marker)
        keep_row.addStretch(1)
        outer.addLayout(keep_row)

        cli_keep_row = QHBoxLayout()
        self._cli_view_center_keep = QCheckBox("Keep CLI-specified Alt/Az", self)
        cli_view_center_enabled = (
            self._cli_view_center_alt_specified or self._cli_view_center_az_specified
        )
        self._cli_view_center_keep.setEnabled(cli_view_center_enabled)
        self._cli_view_center_keep.setChecked(cli_view_center_enabled)
        cli_keep_row.addWidget(self._cli_view_center_keep)
        cli_keep_row.addStretch(1)
        outer.addLayout(cli_keep_row)

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

    def selected_target(self) -> SearchJumpTarget | None:
        if self._clear_persistent_marker_on_accept:
            return None
        item = self._list.currentItem()
        if item is None:
            return None
        target = item.data(Qt.ItemDataRole.UserRole)
        if not isinstance(target, SearchJumpTarget):
            return None
        return replace(
            target,
            persistent_keep_marker=self._jpl_keep_marker.isChecked(),
            preserve_cli_view_center=self._cli_view_center_keep.isChecked(),
        )

    def accept(self) -> None:
        query = self._search.text().strip()
        if not query:
            self._clear_persistent_marker_on_accept = True
            super().accept()
            return
        target = self.selected_target()
        if target is not None:
            super().accept()
            return
        if self._local_result_count > 0:
            self._select_first_visible()
            if self.selected_target() is not None:
                super().accept()
            return
        self._start_jpl_search()

    def should_clear_persistent_marker(self) -> bool:
        return self._clear_persistent_marker_on_accept

    def _apply_local_filter(self, text: str) -> None:
        query_text = text.strip()
        query_spec = parse_search_query(text)
        self._list.clear()

        if not query_text:
            matching_targets = list(self._local_targets)
        else:
            matching_targets = [
                target
                for target in self._local_targets
                if search_target_matches_query(target, query_spec)
            ]

        self._local_result_count = 0 if not query_text else len(matching_targets)
        for target in matching_targets:
            suffix = f"  ({target.subtitle})" if target.subtitle else ""
            item = QListWidgetItem(f"{target.label}{suffix}", self._list)
            item.setData(Qt.ItemDataRole.UserRole, target)
            item.setHidden(not query_text)

        if self._local_result_count > 0:
            self._set_status(
                f"Found {self._local_result_count} local result(s). Select one and press OK."
            )
            self._select_first_visible()
        elif query_text and query_spec.normalized in _JPL_BYPASS_QUERIES:
            self._set_status(
                "Solar-system bodies are already handled by the solar-system view."
            )
        elif query_text and query_spec.normalized:
            self._set_status("No local match. Use the satellites / JPL button below.")
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
        query_spec = parse_search_query(query)
        enabled = (
            bool(query)
            and query_spec.normalized not in _JPL_BYPASS_QUERIES
            and not self._jpl_search_in_progress
        )
        self._jpl_search_button.setEnabled(enabled)
        if not query:
            self._jpl_search_button.setText("Search satellites / JPL (up to 500)")
        elif query_spec.normalized in _JPL_BYPASS_QUERIES:
            self._jpl_search_button.setText("Solar-system bodies are already shown")
        else:
            self._jpl_search_button.setText(
                f"Search satellites / JPL (up to 500) for '{query}'"
            )

    def _sync_ok_button(self) -> None:
        if self._ok_button is not None:
            self._ok_button.setEnabled(self.selected_target() is not None)

    def _on_item_double_clicked(self, _item: QListWidgetItem) -> None:
        self.accept()

    def _start_jpl_search(self) -> None:
        query = self._search.text().strip()
        query_spec = parse_search_query(query)
        if not query or self._jpl_search_in_progress:
            return
        if query_spec.normalized in _JPL_BYPASS_QUERIES:
            self._set_status(
                "Solar-system bodies are already handled by the solar-system view."
            )
            self._sync_jpl_button()
            return
        if self._jpl_search_callback is None:
            self._set_status("Satellite/JPL search is unavailable.")
            self._sync_ok_button()
            self._sync_jpl_button()
            return
        self._jpl_search_request_id += 1
        request_id = self._jpl_search_request_id
        self._jpl_search_in_progress = True
        if self._ok_button is not None:
            self._ok_button.setEnabled(False)
        self._sync_jpl_button()
        self._set_status(f"Searching satellites / JPL for '{query}'...")
        self._services.submit(self._run_jpl_search, request_id=request_id, query=query)

    def _run_jpl_search(self, request_id: int, query: str) -> None:
        try:
            if self._satellite_search_callback is not None:
                satellite_targets = list(self._satellite_search_callback(query))
                if satellite_targets:
                    self.jpl_search_finished.emit(
                        request_id,
                        satellite_targets,
                        f"Found {len(satellite_targets)} satellite result(s)",
                    )
                    return
            targets = (
                list(self._jpl_search_callback(query))
                if self._jpl_search_callback is not None
                else []
            )
        except Exception as exc:  # pragma: no cover - signal delivery path
            self.jpl_search_finished.emit(
                request_id,
                [],
                f"Satellite/JPL search failed: {exc}",
            )
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
