from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

import zstarview.gui.window as window_module
from zstarview.__about__ import __version__
from zstarview.gui.license_dialog import LicenseDialog, load_license_markdown

_app = QApplication.instance() or QApplication([])


def test_open_code_data_licenses_and_credits_opens_local_dialog(monkeypatch) -> None:
    events: list[str] = []

    class FakeLicenseDialog:
        def exec(self) -> None:
            events.append("exec")

    monkeypatch.setattr(window_module, "LicenseDialog", FakeLicenseDialog)

    window_module.open_code_data_licenses_and_credits()

    assert events == ["exec"]


def test_license_dialog_is_selectable_and_copies_all_text() -> None:
    dialog = LicenseDialog()
    flags = dialog.browser.textInteractionFlags()

    assert flags & Qt.TextInteractionFlag.TextSelectableByMouse
    assert flags & Qt.TextInteractionFlag.TextSelectableByKeyboard
    assert flags & Qt.TextInteractionFlag.LinksAccessibleByMouse
    assert flags & Qt.TextInteractionFlag.LinksAccessibleByKeyboard
    assert dialog.browser.openExternalLinks() is True
    assert f"zstarview {__version__}" in dialog.browser.toPlainText()
    assert "Open-Meteo Weather Forecast API" in dialog.browser.toPlainText()

    QApplication.clipboard().clear()
    dialog.copy_all_button.click()

    assert QApplication.clipboard().text() == dialog.browser.toPlainText()


def test_bundled_license_markdown_contains_release_version() -> None:
    markdown = load_license_markdown()

    assert f"zstarview {__version__}" in markdown
    assert "{{ZSTARVIEW_VERSION}}" not in markdown


def test_open_open_meteo_terms_launches_browser(monkeypatch) -> None:
    captured: dict[str, str] = {}

    def fake_open_url(url) -> bool:
        captured["url"] = url.toString()
        return True

    monkeypatch.setattr(window_module.QDesktopServices, "openUrl", fake_open_url)

    window_module.open_open_meteo_terms()

    assert captured["url"] == window_module.OPEN_METEO_TERMS_URL
