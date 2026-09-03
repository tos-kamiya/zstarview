from __future__ import annotations

from types import SimpleNamespace

from PySide6.QtCore import QEvent, Qt
from PySide6.QtGui import QKeyEvent
from PySide6.QtWidgets import QApplication

from zstarview.gui.famous_star_search_dialog import NamedStarSearchDialog
from zstarview.gui.famous_star_shortcuts import SearchJumpTarget
from zstarview.gui.place_search_dialog import PlaceSearchDialog

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


def test_place_search_dialog_disables_cli_view_center_checkbox_without_cli_axes() -> None:
    dialog = PlaceSearchDialog(
        lambda _query, _countrycode, _language: [],
        cli_view_center_alt_specified=False,
        cli_view_center_az_specified=False,
    )

    assert dialog._cli_view_center_keep.isEnabled() is False
    assert dialog._cli_view_center_keep.isChecked() is False


def test_place_search_dialog_enables_cli_view_center_checkbox_with_cli_axes() -> None:
    dialog = PlaceSearchDialog(
        lambda _query, _countrycode, _language: [],
        cli_view_center_alt_specified=True,
        cli_view_center_az_specified=False,
    )

    assert dialog._cli_view_center_keep.isEnabled() is True
    assert dialog._cli_view_center_keep.isChecked() is True


def test_place_search_dialog_selected_target_carries_cli_view_center_choice() -> None:
    dialog = PlaceSearchDialog(
        lambda _query, _countrycode, _language: [],
        cli_view_center_alt_specified=True,
        cli_view_center_az_specified=False,
    )
    dialog._cli_view_center_keep.setChecked(False)
    dialog._place_targets = [
        SearchJumpTarget(
            label="Tokyo Station",
            kind="place",
            sort_key=(0.0, "tokyo station"),
            subtitle="Place / railway / station",
            latitude_deg=35.681236,
            longitude_deg=139.767125,
        )
    ]
    dialog._rebuild_list()

    target = dialog.selected_target()

    assert target is not None
    assert target.preserve_cli_view_center is False


def test_place_search_dialog_passes_query_countrycode_and_language_to_callback(monkeypatch) -> None:
    services = SimpleNamespace(
        submit=lambda target, *args, **kwargs: target(*args, **kwargs)
    )
    seen: list[tuple[str, str | None, str]] = []

    dialog = PlaceSearchDialog(
        lambda query, countrycode, language: seen.append((query, countrycode, language)) or [],
        services=services,
        initial_query="Matsue Station",
        initial_countrycode="jp",
        initial_language="ja",
    )

    dialog._start_place_search()

    assert seen == [("Matsue Station", "jp", "ja")]


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


def test_named_star_search_dialog_disables_cli_view_center_checkbox_without_cli_axes() -> None:
    dialog = NamedStarSearchDialog(
        [
            SearchJumpTarget(
                label="Sirius",
                kind="star",
                sort_key=(0.0, "sirius"),
                subtitle="Vmag -1.44",
            )
        ],
        cli_view_center_alt_specified=False,
        cli_view_center_az_specified=False,
        jpl_search_callback=lambda _query: [],
    )

    assert dialog._cli_view_center_keep.isEnabled() is False
    assert dialog._cli_view_center_keep.isChecked() is False


def test_named_star_search_dialog_enables_cli_view_center_checkbox_with_cli_axes() -> None:
    dialog = NamedStarSearchDialog(
        [
            SearchJumpTarget(
                label="Sirius",
                kind="star",
                sort_key=(0.0, "sirius"),
                subtitle="Vmag -1.44",
            )
        ],
        cli_view_center_alt_specified=True,
        cli_view_center_az_specified=False,
        jpl_search_callback=lambda _query: [],
    )

    assert dialog._cli_view_center_keep.isEnabled() is True
    assert dialog._cli_view_center_keep.isChecked() is True


def test_named_star_search_dialog_selected_target_carries_cli_view_center_choice() -> None:
    dialog = NamedStarSearchDialog(
        [
            SearchJumpTarget(
                label="Sirius",
                kind="star",
                sort_key=(0.0, "sirius"),
                subtitle="Vmag -1.44",
            )
        ],
        cli_view_center_alt_specified=True,
        cli_view_center_az_specified=False,
        jpl_search_callback=lambda _query: [],
    )
    dialog._cli_view_center_keep.setChecked(False)
    dialog._search.setText("siri")

    target = dialog.selected_target()

    assert target is not None
    assert target.preserve_cli_view_center is False


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
    assert target.preserve_cli_view_center is False


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
