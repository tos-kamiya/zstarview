from __future__ import annotations

from PySide6.QtWidgets import QDoubleSpinBox
from PySide6.QtWidgets import QApplication

from zstarview.gui.startup_dialog import _coerce_float_value, _coerce_int_value
from zstarview.gui.startup_dialog import StartupDialog

_app = QApplication.instance() or QApplication([])


def test_startup_dialog_numeric_coercion_uses_fallbacks() -> None:
    assert _coerce_float_value("observer_height_m", None, 0.0) == 1.7
    assert _coerce_float_value("other", None, 2.5) == 2.5
    assert _coerce_int_value("other", None, 7) == 7


def test_startup_dialog_tabs_follow_requested_order() -> None:
    dialog = StartupDialog()

    assert dialog.width() == 640
    assert dialog.height() == 380
    assert dialog._tabs.count() == 5
    assert dialog._tabs.tabText(0) == "Location & Time"
    assert dialog._tabs.tabText(1) == "Stars"
    assert dialog._tabs.tabText(2) == "Overlays"
    assert dialog._tabs.tabText(3) == "General"
    assert dialog._tabs.tabText(4) == "Search Objects at Startup"
    assert set(dialog._overlay_sections) == {
        "Sky",
        "Clouds",
        "Aircraft and Satellites",
        "Ground and Guides",
        "Urban Outline",
    }
    assert dialog._overlay_sections["Sky"].is_expanded() is True
    dialog._overlay_sections["Sky"]._button.setChecked(False)
    assert dialog._overlay_sections["Sky"].is_expanded() is False
    assert "overlay_font_size" in dialog._widgets
    assert dialog._reset_button.text() == "Reset to Default Values"
    terrain_widget = dialog._widgets["terrain_horizon_opacity"]
    assert isinstance(terrain_widget, QDoubleSpinBox)
    assert terrain_widget.decimals() == 3
    assert terrain_widget.value() == 0.003
