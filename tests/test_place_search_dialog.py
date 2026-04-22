from __future__ import annotations

from types import SimpleNamespace

from PySide6.QtCore import QEvent, Qt
from PySide6.QtGui import QKeyEvent
from PySide6.QtWidgets import QApplication

from zstarview.gui.place_search_dialog import PlaceSearchDialog
from zstarview.gui.famous_star_search_dialog import NamedStarSearchDialog
from zstarview.gui.famous_star_shortcuts import SearchJumpTarget

_app = QApplication.instance() or QApplication([])


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


def test_named_star_search_dialog_uses_local_match_before_jpl() -> None:
    dialog = NamedStarSearchDialog(
        [
            SearchJumpTarget(
                label="Sirius",
                kind="star",
                sort_key=(0.0, "sirius"),
                subtitle="Vmag -1.44",
            )
        ],
        jpl_search_callback=lambda _query: [],
    )
    called: list[str] = []
    dialog._start_jpl_search = lambda: called.append("jpl")  # type: ignore[method-assign]

    dialog._search.setText("siri")
    dialog.accept()

    assert called == []
    assert dialog.selected_target() is not None
    assert dialog.selected_target().label == "Sirius"


def test_named_star_search_dialog_falls_back_to_jpl_when_local_empty() -> None:
    dialog = NamedStarSearchDialog(
        [
            SearchJumpTarget(
                label="Sirius",
                kind="star",
                sort_key=(0.0, "sirius"),
                subtitle="Vmag -1.44",
            )
        ],
        jpl_search_callback=lambda _query: [],
    )
    called: list[str] = []
    dialog._start_jpl_search = lambda: called.append("jpl")  # type: ignore[method-assign]

    dialog._search.setText("ceres")
    dialog.accept()

    assert called == ["jpl"]
    assert "satellites / JPL" in dialog._jpl_search_button.text()


def test_named_star_search_dialog_keeps_jpl_button_enabled_for_local_matches() -> None:
    dialog = NamedStarSearchDialog(
        [
            SearchJumpTarget(
                label="Sirius",
                kind="star",
                sort_key=(0.0, "sirius"),
                subtitle="Vmag -1.44",
            )
        ],
        jpl_search_callback=lambda _query: [],
    )

    dialog._search.setText("siri")

    assert dialog._jpl_search_button.isEnabled() is True


def test_named_star_search_dialog_applies_keep_marker_to_jpl_result() -> None:
    dialog = NamedStarSearchDialog(
        [
            SearchJumpTarget(
                label="Mars",
                kind="jpl_body",
                sort_key=(0.0, "mars"),
                subtitle="major body",
            )
        ],
        jpl_search_callback=lambda _query: [],
    )
    dialog._jpl_keep_marker.setChecked(True)
    item = dialog._list.item(0)
    assert item is not None
    dialog._list.setCurrentItem(item)

    target = dialog.selected_target()

    assert target is not None
    assert target.persistent_keep_marker is True


def test_named_star_search_dialog_blocks_solar_system_body_jpl_fallback() -> None:
    dialog = NamedStarSearchDialog(
        [],
        jpl_search_callback=lambda _query: [],
    )
    called: list[str] = []
    dialog._start_jpl_search = lambda: called.append("jpl")  # type: ignore[method-assign]

    dialog._search.setText("Mars")

    assert called == []
    assert dialog._jpl_search_button.isEnabled() is False
    assert "already shown" in dialog._jpl_search_button.text()


def test_named_star_search_dialog_empty_query_clears_persistent_marker() -> None:
    dialog = NamedStarSearchDialog(
        [
            SearchJumpTarget(
                label="Sirius",
                kind="star",
                sort_key=(0.0, "sirius"),
                subtitle="Vmag -1.44",
            )
        ],
        jpl_search_callback=lambda _query: [],
    )
    dialog._jpl_keep_marker.setChecked(True)
    dialog._search.setText("")

    dialog.accept()

    assert dialog.should_clear_persistent_marker() is True
    assert dialog.selected_target() is None
