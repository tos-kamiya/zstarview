from __future__ import annotations

from zstarview.gui.startup_dialog import _coerce_float_value, _coerce_int_value


def test_startup_dialog_numeric_coercion_uses_fallbacks() -> None:
    assert _coerce_float_value("observer_height_m", None, 0.0) == 1.7
    assert _coerce_float_value("other", None, 2.5) == 2.5
    assert _coerce_int_value("other", None, 7) == 7
