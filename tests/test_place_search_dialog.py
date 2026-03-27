from __future__ import annotations

from types import SimpleNamespace

from PySide6.QtCore import QEvent, Qt
from PySide6.QtGui import QKeyEvent

from zstarview.gui.place_search_dialog import PlaceSearchDialog


def test_place_search_dialog_enter_in_search_field_triggers_search() -> None:
    calls: list[str] = []
    dummy = SimpleNamespace(
        _search=object(),
        _start_place_search=lambda: calls.append("search"),
    )

    handled = PlaceSearchDialog.eventFilter(
        dummy,
        dummy._search,
        QKeyEvent(QEvent.Type.KeyPress, Qt.Key.Key_Return, Qt.KeyboardModifier.NoModifier),
    )

    assert handled is True
    assert calls == ["search"]
