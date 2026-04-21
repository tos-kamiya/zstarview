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
        "Code, Data Licenses, and Credits...",
    ]
    assert menu.actions()[0].isEnabled() is False
    assert menu.actions()[1].isEnabled() is True


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
