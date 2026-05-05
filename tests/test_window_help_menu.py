from __future__ import annotations

from PySide6.QtWidgets import QApplication, QMenu

import zstarview.gui.window as window_module

_app = QApplication.instance() or QApplication([])


class _Recorder:
    def __init__(self) -> None:
        self.created: list[tuple[str, bool]] = []

    def _add_menu_action(self, menu: QMenu, text: str, **kwargs):
        action = window_module.QAction(text, menu)
        action.setEnabled(kwargs.get("enabled", True))
        menu.addAction(action)
        self.created.append((text, action.isEnabled()))
        return action


def test_add_help_menu_actions_adds_version_and_credits() -> None:
    recorder = _Recorder()
    menu = QMenu("Help")

    window_module.SkyWindowCoreMixin._add_help_menu_actions(recorder, menu)

    texts = [action.text() for action in menu.actions()]
    assert texts == [
        f"Version {window_module.__version__}",
        "Code, Data Licenses, and Credits ->",
    ]
    assert menu.actions()[0].isEnabled() is False
    assert menu.actions()[1].isEnabled() is True


def test_attach_help_menu_to_file_menu_places_help_before_exit() -> None:
    recorder = _Recorder()
    file_menu = QMenu("File")
    help_menu = QMenu("Help")

    recorder.file_menu = file_menu
    recorder.help_menu = help_menu

    window_module.SkyWindowCoreMixin._attach_help_menu_to_file_menu(recorder)

    actions = file_menu.actions()
    assert [action.text() for action in actions] == ["", "Help", ""]
    assert actions[1].menu() is help_menu


def test_attach_help_menu_to_frameless_menu_places_help_before_exit() -> None:
    recorder = _Recorder()
    menu = QMenu("Menu")
    help_menu = QMenu("Help")
    square_action = window_module.QAction("Square Client Area", menu)
    default_size_action = window_module.QAction("Default Window Size", menu)
    fullscreen_action = window_module.QAction("Fullscreen", menu)
    exit_action = window_module.QAction("Exit", menu)

    recorder.menu = menu
    recorder.help_menu = help_menu

    window_module.SkyWindowCoreMixin._attach_help_menu_to_frameless_menu(
        recorder,
        square_action,
        default_size_action,
        fullscreen_action,
        exit_action,
    )

    actions = menu.actions()
    assert [action.text() for action in actions] == [
        "",
        "Square Client Area",
        "Default Window Size",
        "Fullscreen",
        "",
        "Help",
        "",
        "Exit",
    ]
    assert actions[5].menu() is help_menu


def test_open_code_data_licenses_and_credits_launches_browser(monkeypatch) -> None:
    captured: dict[str, str] = {}

    def fake_open_url(url) -> bool:
        captured["url"] = url.toString()
        return True

    monkeypatch.setattr(window_module.QDesktopServices, "openUrl", fake_open_url)

    window_module.open_code_data_licenses_and_credits()

    assert (
        captured["url"]
        == window_module.GITHUB_CODE_DATA_LICENSES_AND_CREDITS_URL
    )
