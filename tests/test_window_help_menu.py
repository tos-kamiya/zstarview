from __future__ import annotations

from PySide6.QtWidgets import QApplication

import zstarview.gui.window as window_module

_app = QApplication.instance() or QApplication([])


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


def test_open_open_meteo_terms_launches_browser(monkeypatch) -> None:
    captured: dict[str, str] = {}

    def fake_open_url(url) -> bool:
        captured["url"] = url.toString()
        return True

    monkeypatch.setattr(window_module.QDesktopServices, "openUrl", fake_open_url)

    window_module.open_open_meteo_terms()

    assert captured["url"] == window_module.OPEN_METEO_TERMS_URL
