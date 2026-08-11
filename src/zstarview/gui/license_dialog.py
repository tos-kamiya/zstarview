from __future__ import annotations

from importlib.resources import files

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QDialog,
    QDialogButtonBox,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from ..__about__ import __version__

LICENSE_RESOURCE_NAME = "licenses_and_data_sources.md"


def load_license_markdown() -> str:
    """Load the license snapshot bundled with this zstarview release."""
    markdown = (
        files("zstarview.data")
        .joinpath(LICENSE_RESOURCE_NAME)
        .read_text(encoding="utf-8")
    )
    return markdown.replace("{{ZSTARVIEW_VERSION}}", __version__)


class LicenseDialog(QDialog):
    """Display selectable, version-matched license and data-source details."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Licenses and Data Sources")
        self.resize(900, 640)

        layout = QVBoxLayout(self)
        self.browser = QTextBrowser(self)
        self.browser.setReadOnly(True)
        self.browser.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse
            | Qt.TextInteractionFlag.TextSelectableByKeyboard
            | Qt.TextInteractionFlag.LinksAccessibleByMouse
            | Qt.TextInteractionFlag.LinksAccessibleByKeyboard
        )
        self.browser.setOpenExternalLinks(True)
        self.browser.setMarkdown(load_license_markdown())
        layout.addWidget(self.browser, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close, self)
        self.copy_all_button = buttons.addButton(
            "Copy All", QDialogButtonBox.ButtonRole.ActionRole
        )
        self.copy_all_button.clicked.connect(self.copy_all)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def copy_all(self) -> None:
        QApplication.clipboard().setText(self.browser.toPlainText())
