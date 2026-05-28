from __future__ import annotations

from types import SimpleNamespace

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QWidget

import zstarview.gui.window as window_module

_app = QApplication.instance() or QApplication([])


class _TimerStub:
    def __init__(self, active: bool) -> None:
        self._active = active

    def isActive(self) -> bool:  # noqa: N802 - Qt naming
        return self._active

    def stop(self) -> None:
        self._active = False


def test_shutdown_message_overlay_uses_centered_text() -> None:
    parent = QWidget()
    overlay = window_module.ShutdownMessageOverlay(parent)

    assert overlay.parent() is parent
    assert overlay.text() == "Shutting down (closing sub-processes)... please wait."
    assert overlay.alignment() & Qt.AlignmentFlag.AlignCenter
    assert overlay.wordWrap() is True
    style_sheet = overlay.styleSheet()
    assert "border: 3px solid" in style_sheet
    assert "border-left: none" not in style_sheet
    assert "border-right: none" not in style_sheet


def test_shutdown_message_geometry_is_a_centered_band() -> None:
    dummy = SimpleNamespace(
        client_width=lambda: 800,
        height=lambda: 600,
        _shutdown_overlay=SimpleNamespace(sizeHint=lambda: SimpleNamespace(height=lambda: 72)),
    )

    geometry = window_module.SkyWindowCoreMixin._shutdown_message_geometry(dummy)

    assert geometry.x() == 20
    assert geometry.width() == 760
    assert geometry.height() == 72
    assert geometry.y() == 264


def test_request_application_quit_shows_shutdown_message_before_quit(
    monkeypatch,
) -> None:
    call_order: list[str] = []

    dummy = SimpleNamespace(
        _show_shutdown_message=lambda: call_order.append("show"),
    )

    monkeypatch.setattr(
        window_module.QApplication,
        "quit",
        lambda: call_order.append("quit"),
    )

    window_module.SkyWindowCoreMixin._request_application_quit(dummy)

    assert call_order == ["show", "quit"]


def test_begin_shutdown_shows_message_around_controller_shutdown() -> None:
    call_order: list[str] = []

    dummy = SimpleNamespace(
        _is_shutting_down=False,
        _show_shutdown_message=lambda: call_order.append("show"),
        _hide_shutdown_message=lambda: call_order.append("hide"),
        _sky_worker=SimpleNamespace(
            shutdown=lambda: call_order.append("worker"),
        ),
        _cloud_controller=None,
        _geosatellite_controller=None,
        _satellite_controller=None,
        _aircraft_controller=None,
        _jpl_small_body_controller=None,
        _terrain_horizon_controller=None,
        _water_overlay_controller=None,
        _urban_outline_controller=None,
        _sky_data_update_timer=_TimerStub(False),
        _asterism_check_timer=_TimerStub(False),
        _cloud_update_timer=_TimerStub(False),
        _satellite_update_timer=_TimerStub(False),
        _overlay_projection_timer=_TimerStub(False),
        _aircraft_update_timer=_TimerStub(False),
        _persistent_search_update_timer=_TimerStub(False),
        _interaction_idle_timer=_TimerStub(False),
        _viewport_interaction_idle_timer=_TimerStub(False),
    )

    window_module.SkyWindowCoreMixin._begin_shutdown(dummy)

    assert call_order == ["show", "worker", "hide"]
