from types import SimpleNamespace

import pytest
from PySide6.QtCore import QPoint, QRect

import zstarview.gui.window as window_module
from zstarview.gui.window import SkyWindowCoreMixin


@pytest.mark.parametrize("size", [(800, 600), (2400, 1600)])
def test_outside_cursor_fades_even_when_global_position_is_stale(
    monkeypatch, size
) -> None:
    # Wayland can keep reporting the last in-window global position after leave.
    monkeypatch.setattr(window_module.QCursor, "pos", lambda: QPoint(100, 100))
    clock = [0.0]
    monkeypatch.setattr(window_module.time, "monotonic", lambda: clock[0])
    window = SimpleNamespace(
        _landscape_cursor_inside=False,
        _landscape_annotation_opacity=1.0,
        _landscape_fade_was_active=True,
        _landscape_fade_last_tick=0.0,
        _landscape_cursor_at_edge_or_outside=lambda: (
            SkyWindowCoreMixin._landscape_cursor_at_edge_or_outside(window)
        ),
        _landscape_mode=lambda: True,
        frameGeometry=lambda: QRect(0, 0, *size),
        request_client_update=lambda: None,
    )
    window._landscape_fade_triggered = lambda: (
        SkyWindowCoreMixin._landscape_fade_triggered(window)
    )
    for tick in range(1, 52):
        clock[0] = tick / 10.0
        SkyWindowCoreMixin._update_landscape_annotation_fade(window)
        if tick == 25:
            assert window._landscape_annotation_opacity == pytest.approx(0.5)
    assert window._landscape_annotation_opacity == 0.0

    window._landscape_cursor_inside = True
    window._landscape_fade_triggered = lambda: False
    clock[0] += 0.1
    SkyWindowCoreMixin._update_landscape_annotation_fade(window)
    assert window._landscape_annotation_opacity == 1.0
