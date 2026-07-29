"""Dialog for selecting a named star to jump to."""
from __future__ import annotations

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
    SATELLITE_JUMP_SHORTCUTS,
    NamedStarShortcut,
)


class NamedStarJumpDialog(QDialog):
    def __init__(self, stars_by_band: dict[str, list[NamedStarShortcut]], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Jump to Named Star")
        self.setModal(True)
        self.resize(460, 420)

        layout = QVBoxLayout(self)
        self._tabs = QTabWidget(self)
        self._lists: dict[str, QListWidget] = {}

        tab_defs = [
            (DEC_BAND_NORTH, "North"),
            (DEC_BAND_EQUATOR, "Equatorial"),
            (DEC_BAND_SOUTH, "South"),
            ("satellites", "Artificial Satellites"),
        ]
        for key, label in tab_defs:
            lw = QListWidget(self)
            lw.itemDoubleClicked.connect(self._on_item_double_clicked)
            self._lists[key] = lw
            stars = list(stars_by_band.get(key, [])) if key != "satellites" else list(SATELLITE_JUMP_SHORTCUTS)
            self._tabs.addTab(lw, f"{label} ({len(stars)})")
            for star in stars:
                suffix = star.subtitle or (f"Vmag {star.vmag:.2f}" if star.kind == "star" else "")
                item = QListWidgetItem(f"{star.name}  ({suffix})" if suffix else star.name, lw)
                item.setData(Qt.ItemDataRole.UserRole, star)

        layout.addWidget(self._tabs)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel, parent=self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def selected_star(self) -> NamedStarShortcut | None:
        widget = self._tabs.currentWidget()
        if isinstance(widget, QListWidget):
            item = widget.currentItem()
            if item is not None:
                star = item.data(Qt.ItemDataRole.UserRole)
                if isinstance(star, NamedStarShortcut):
                    return star
        for list_widget in self._lists.values():
            item = list_widget.currentItem()
            if item is None:
                continue
            star = item.data(Qt.ItemDataRole.UserRole)
            if isinstance(star, NamedStarShortcut):
                return star
        return None

    def _on_item_double_clicked(self, _item: QListWidgetItem) -> None:
        self.accept()
