
from types import SimpleNamespace
from unittest.mock import Mock

from PySide6.QtCore import QPoint, Qt
from PySide6.QtWidgets import QWidget

import zstarview.gui.window as window_module
import zstarview.gui.window_render as window_render_module
import zstarview.gui.window_updates as window_updates_module
from tests._window_render_support import _DummyTimer, _WindowStub
from zstarview.gui.draggable_window import DraggableWindow
from zstarview.gui.window import SkyWindow
from zstarview.gui.window_state import SkyWindowState
from zstarview.types import ViewerData


def test_handle_client_key_press_rotates_view_immediately() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy._sync_view_altitude_actions = Mock()
    dummy.request_client_update = Mock()
    dummy._viewport_interaction_idle_timer = _DummyTimer(active=False)
    dummy._rotate_view = lambda *args, **kwargs: SkyWindow._rotate_view(
        dummy, *args, **kwargs
    )
    dummy._begin_viewport_interaction_mode = lambda *args, **kwargs: (
        SkyWindow._begin_viewport_interaction_mode(dummy, *args, **kwargs)
    )
    dummy._update_viewport_interaction_stars = lambda: (
        SkyWindow._update_viewport_interaction_stars(dummy)
    )
    dummy._viewport_rotation_keys = lambda: SkyWindow._viewport_rotation_keys(dummy)

    event = SimpleNamespace(
        key=lambda: Qt.Key.Key_Left,
        modifiers=lambda: Qt.KeyboardModifier.NoModifier,
        isAutoRepeat=lambda: False,
        accept=Mock(),
    )

    SkyWindow._handle_client_key_press(dummy, event)

    assert dummy.viewer_data.view_center == (20.0, 25.0)
    assert dummy.state.render_view_center == (20.0, 25.0)
    assert dummy.state.viewport_interaction_mode is True
    assert dummy._viewport_interaction_idle_timer.started_with == []
    assert dummy._viewport_rotation_keys_down == {Qt.Key.Key_Left}
    dummy._sync_view_altitude_actions.assert_called_once()
    dummy.request_client_update.assert_called_once()
    event.accept.assert_called_once()


def test_set_view_center_leaves_viewport_fast_mode_after_dialog_change() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=True,
        viewport_interaction_stars=object(),
    )
    dummy._sync_view_altitude_actions = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.request_client_update = Mock()
    dummy.start_background_cloud_update = Mock()
    dummy.start_background_terrain_horizon_update = Mock()
    dummy._end_viewport_interaction_mode = lambda *args, **kwargs: (
        SkyWindow._end_viewport_interaction_mode(dummy, *args, **kwargs)
    )
    dummy._finalize_view_direction_change = lambda: (
        SkyWindow._finalize_view_direction_change(dummy)
    )
    dummy.request_cloud_projection_update = Mock()

    SkyWindow._set_view_center(
        dummy,
        25.0,
        45.0,
        interactive_viewport=False,
        start_viewport_idle_timer=False,
    )

    assert dummy.viewer_data.view_center == (25.0, 45.0)
    assert dummy.state.render_view_center == (25.0, 45.0)
    assert dummy.state.viewport_interaction_mode is False
    assert dummy.state.viewport_interaction_stars is None
    dummy.request_sky_data_update.assert_called_once_with(
        reason="view-change-idle",
        allow_during_viewport_interaction=True,
    )
    dummy.request_cloud_projection_update.assert_called_once_with(
        reason="view-change-idle"
    )
    dummy.start_background_terrain_horizon_update.assert_called_once_with(
        reason="view-change-idle"
    )
    dummy.request_client_update.assert_called_once()
    dummy._sync_view_altitude_actions.assert_called_once()


def test_open_view_direction_dialog_shows_fast_frame_before_release(
    monkeypatch,
) -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=False,
    )
    set_calls: list[tuple[tuple[float, float], dict[str, object]]] = []
    end_calls: list[str] = []

    class _FakeDialog:
        def __init__(self, view_center: tuple[float, float], parent) -> None:
            self._view_center = view_center
            self._parent = parent

        def exec(self) -> int:
            return 1

        def selected_view_center(self) -> tuple[float, float]:
            return (25.0, 45.0)

    monkeypatch.setattr(window_module, "ViewDirectionDialog", _FakeDialog)
    def fake_single_shot(ms, func) -> None:
        end_calls.append(f"timer:{ms}")
        func()

    monkeypatch.setattr(window_module.QTimer, "singleShot", fake_single_shot)
    dummy._set_view_center = lambda *args, **kwargs: set_calls.append(
        (args, dict(kwargs))
    )
    dummy._end_viewport_interaction_mode = lambda *_args, **kwargs: end_calls.append(
        str(kwargs.get("reason"))
    )
    dummy._finalize_view_direction_change = lambda: (
        SkyWindow._finalize_view_direction_change(dummy)
    )
    dummy._finalize_view_direction_dialog_change = lambda: (
        SkyWindow._finalize_view_direction_dialog_change(dummy)
    )

    SkyWindow._open_view_direction_dialog(dummy)

    assert set_calls == [
        (
            (25.0, 45.0),
            {
                "interactive_viewport": True,
                "start_viewport_idle_timer": False,
            },
        )
    ]
    assert end_calls == ["timer:0", "view-change-release"]


def test_handle_client_key_release_ends_viewport_interaction_mode() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=True,
    )
    dummy._viewport_rotation_keys_down = {Qt.Key.Key_Left}
    dummy._viewport_rotation_keys = lambda: SkyWindow._viewport_rotation_keys(dummy)
    dummy._end_viewport_interaction_mode = lambda *args, **kwargs: (
        SkyWindow._end_viewport_interaction_mode(dummy, *args, **kwargs)
    )
    dummy.request_sky_data_update = Mock()
    dummy.start_background_cloud_update = Mock()
    dummy.start_background_terrain_horizon_update = Mock()
    dummy.reproject_tropical_cyclone_overlay = Mock()
    dummy.request_client_update = Mock()

    event = SimpleNamespace(
        key=lambda: Qt.Key.Key_Left,
        isAutoRepeat=lambda: False,
        accept=Mock(),
    )

    SkyWindow._handle_client_key_release(dummy, event)

    assert dummy._viewport_rotation_keys_down == set()
    assert dummy.state.viewport_interaction_release_pending is True
    assert dummy.state.viewport_interaction_mode is True
    dummy.request_sky_data_update.assert_called_once_with(
        reason="viewport-interaction-release",
        allow_during_viewport_interaction=True,
    )
    dummy.start_background_cloud_update.assert_not_called()
    dummy.start_background_terrain_horizon_update.assert_not_called()
    dummy.request_client_update.assert_not_called()
    event.accept.assert_called_once()


def test_handle_client_mouse_move_is_ignored_during_startup_block() -> None:
    dummy = _WindowStub()
    dummy._startup_input_blocked = lambda: True
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        mouse_pos=None,
    )
    dummy.request_client_update = Mock()

    event = SimpleNamespace(pos=lambda: QPoint(10, 20), accept=Mock())

    SkyWindow._handle_client_mouse_move(dummy, event)

    assert dummy.state.mouse_pos is None
    dummy.request_client_update.assert_not_called()
    event.accept.assert_called_once()


def test_handle_client_mouse_move_coalesces_repaints() -> None:
    dummy = _WindowStub()
    dummy._startup_input_blocked = lambda: False
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        mouse_pos=None,
    )
    dummy._hover_repaint_timer = _DummyTimer(active=False)
    dummy._note_landscape_cursor_motion = Mock()
    dummy.request_client_update = Mock()

    first_event = SimpleNamespace(pos=lambda: QPoint(10, 20), accept=Mock())
    second_event = SimpleNamespace(pos=lambda: QPoint(30, 40), accept=Mock())

    SkyWindow._handle_client_mouse_move(dummy, first_event)
    SkyWindow._handle_client_mouse_move(dummy, second_event)

    assert dummy.state.mouse_pos == QPoint(30, 40)
    assert dummy._note_landscape_cursor_motion.call_count == 2
    assert dummy._hover_repaint_timer.started_with == [0]
    dummy.request_client_update.assert_not_called()
    first_event.accept.assert_called_once()
    second_event.accept.assert_called_once()


def test_handle_client_leave_cancels_hover_repaint_and_updates() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        mouse_pos=QPoint(10, 20),
    )
    dummy._hover_repaint_timer = _DummyTimer(active=True)
    dummy.request_client_update = Mock()
    event = SimpleNamespace(accept=Mock())

    SkyWindow._handle_client_leave(dummy, event)

    assert dummy.state.mouse_pos is None
    assert dummy._hover_repaint_timer.isActive() is False
    dummy.request_client_update.assert_called_once()
    event.accept.assert_called_once()


def test_render_hud_state_ignores_mouse_position_during_startup_block() -> None:
    dummy = _WindowStub()
    dummy._startup_input_blocked = lambda: True
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        mouse_pos=QPoint(10, 20),
        overlay_info_bottom_left=False,
    )
    dummy.observation_info_pinned = False
    dummy.client_height = lambda: 300
    dummy._status_line_message = lambda: "status"

    hud = window_render_module.SkyWindowRenderMixin._render_hud_state(dummy)

    assert hud.mouse_pos is None
    assert hud.overlay_info_bottom_left is False


def test_end_viewport_interaction_mode_marks_idle_reason() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=True,
    )
    dummy.request_sky_data_update = Mock()
    dummy.request_cloud_projection_update = Mock()
    dummy.start_background_terrain_horizon_update = Mock()
    dummy.reproject_tropical_cyclone_overlay = Mock()
    dummy.request_client_update = Mock()

    SkyWindow._end_viewport_interaction_mode(dummy)

    assert dummy.state.viewport_interaction_mode is False
    assert dummy.state.viewport_interaction_stars is None
    assert dummy.state.viewport_interaction_release_pending is False
    assert dummy.state.viewport_interaction_completion_reason is None
    dummy.request_sky_data_update.assert_called_once_with(
        reason="viewport-interaction-idle",
        allow_during_viewport_interaction=True,
    )
    dummy.request_cloud_projection_update.assert_called_once_with(
        reason="view-change-idle"
    )
    dummy.start_background_terrain_horizon_update.assert_called_once_with(
        reason="view-change-idle"
    )
    dummy.reproject_tropical_cyclone_overlay.assert_called_once_with()
    dummy.request_client_update.assert_called_once()


def test_end_viewport_interaction_mode_release_reprojects_tropical_cyclone() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=True,
    )
    dummy.request_sky_data_update = Mock()
    dummy.request_cloud_projection_update = Mock()
    dummy.start_background_terrain_horizon_update = Mock()
    dummy.reproject_tropical_cyclone_overlay = Mock()
    dummy.request_client_update = Mock()

    SkyWindow._end_viewport_interaction_mode(
        dummy,
        reason="viewport-interaction-release",
    )

    dummy.request_sky_data_update.assert_called_once_with(
        reason="viewport-interaction-release",
        allow_during_viewport_interaction=True,
    )
    dummy.reproject_tropical_cyclone_overlay.assert_called_once_with(
        allow_during_viewport_interaction=True,
    )
    dummy.request_cloud_projection_update.assert_not_called()
    dummy.start_background_terrain_horizon_update.assert_not_called()
    dummy.request_client_update.assert_not_called()


def test_end_viewport_interaction_mode_release_clears_interaction_when_sky_update_is_busy() -> (
    None
):
    dummy = _WindowStub()
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=True,
    )
    dummy.request_sky_data_update = Mock(return_value=False)
    dummy.request_cloud_projection_update = Mock()
    dummy.start_background_terrain_horizon_update = Mock()
    dummy.reproject_tropical_cyclone_overlay = Mock()
    dummy.request_client_update = Mock()

    SkyWindow._end_viewport_interaction_mode(
        dummy,
        reason="viewport-interaction-release",
    )

    assert dummy.state.viewport_interaction_mode is False
    assert dummy.state.viewport_interaction_stars is None
    assert dummy.state.viewport_interaction_release_pending is False
    dummy.request_sky_data_update.assert_called_once_with(
        reason="viewport-interaction-release",
        allow_during_viewport_interaction=True,
    )
    dummy.request_cloud_projection_update.assert_called_once_with(
        reason="view-change-release",
    )
    dummy.start_background_terrain_horizon_update.assert_called_once_with(
        reason="view-change-release",
    )
    dummy.reproject_tropical_cyclone_overlay.assert_called_once_with()
    dummy.request_client_update.assert_called_once()


def test_tropical_cyclone_layer_is_disabled_when_opacity_is_zero() -> None:
    dummy = _WindowStub(
        show_tropical_cyclone_overlay=True,
        tropical_cyclone_opacity=0.0,
        _tropical_cyclone_controller=object(),
    )

    assert (
        window_updates_module.SkyWindowUpdatesMixin._tropical_cyclone_layer_enabled(
            dummy
        )
        is False
    )


def test_background_press_ignores_drag_exclusions() -> None:
    class _ProbeWindow(DraggableWindow):
        pass

    probe = _ProbeWindow()
    root = QWidget()
    root.resize(200, 200)
    root.show()
    excluded = QWidget(root)
    excluded.setGeometry(160, 160, 20, 20)
    excluded.show()
    probe.add_drag_exclusion(excluded)

    event = SimpleNamespace(
        button=lambda: Qt.MouseButton.LeftButton,
        position=lambda: SimpleNamespace(toPoint=lambda: QPoint(170, 170)),
        globalPosition=lambda: SimpleNamespace(toPoint=lambda: QPoint(170, 170)),
        accept=Mock(),
    )

    assert probe._begin_drag(None, event, root=root) is False
    assert probe._drag_press_pending is False


def _rotation_stub() -> _WindowStub:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy._sync_view_altitude_actions = Mock()
    dummy.request_client_update = Mock()
    dummy._viewport_interaction_idle_timer = _DummyTimer(active=False)
    dummy._rotate_view = lambda *args, **kwargs: SkyWindow._rotate_view(
        dummy, *args, **kwargs
    )
    dummy._begin_viewport_interaction_mode = lambda *args, **kwargs: (
        SkyWindow._begin_viewport_interaction_mode(dummy, *args, **kwargs)
    )
    dummy._update_viewport_interaction_stars = lambda: None
    dummy._viewport_rotation_keys = lambda: SkyWindow._viewport_rotation_keys(dummy)
    return dummy


def _arrow_event(key: Qt.Key, *, auto_repeat: bool) -> SimpleNamespace:
    return SimpleNamespace(
        key=lambda: key,
        modifiers=lambda: Qt.KeyboardModifier.NoModifier,
        isAutoRepeat=lambda: auto_repeat,
        accept=Mock(),
    )


def test_arrow_auto_repeat_rotates_once_per_painted_frame() -> None:
    dummy = _rotation_stub()

    SkyWindow._handle_client_key_press(dummy, _arrow_event(Qt.Key.Key_Right, auto_repeat=False))
    assert dummy.viewer_data.view_center == (20.0, 35.0)

    # Repeats before the next paint are accepted but do not rotate; each
    # one re-requests the paint.
    dummy.request_client_update.reset_mock()
    for _ in range(3):
        event = _arrow_event(Qt.Key.Key_Right, auto_repeat=True)
        SkyWindow._handle_client_key_press(dummy, event)
        event.accept.assert_called_once()
    assert dummy.viewer_data.view_center == (20.0, 35.0)
    assert dummy.request_client_update.call_count == 3

    # After a paint, one repeat rotates one more step.
    dummy._viewport_rotation_paint_pending = False
    SkyWindow._handle_client_key_press(dummy, _arrow_event(Qt.Key.Key_Right, auto_repeat=True))
    SkyWindow._handle_client_key_press(dummy, _arrow_event(Qt.Key.Key_Right, auto_repeat=True))
    assert dummy.viewer_data.view_center == (20.0, 40.0)


def test_arrow_initial_press_rotates_even_when_paint_is_pending() -> None:
    dummy = _rotation_stub()
    dummy._viewport_rotation_paint_pending = True

    SkyWindow._handle_client_key_press(dummy, _arrow_event(Qt.Key.Key_Up, auto_repeat=False))
    SkyWindow._handle_client_key_press(dummy, _arrow_event(Qt.Key.Key_Left, auto_repeat=False))

    assert dummy.viewer_data.view_center == (25.0, 25.0)
    assert dummy._viewport_rotation_keys_down == {Qt.Key.Key_Up, Qt.Key.Key_Left}


def test_paint_event_allows_next_arrow_repeat() -> None:
    dummy = _rotation_stub()
    dummy._viewport_rotation_paint_pending = True
    dummy._startup_splash_visible = lambda: True

    from zstarview.gui.window_render import SkyWindowRenderMixin

    SkyWindowRenderMixin.paintEvent(dummy, None)  # type: ignore[arg-type]

    assert dummy._viewport_rotation_paint_pending is False


def test_arrow_auto_repeat_rotates_when_paint_does_not_arrive() -> None:
    dummy = _rotation_stub()
    SkyWindow._handle_client_key_press(dummy, _arrow_event(Qt.Key.Key_Right, auto_repeat=False))
    # Simulate a stale pending flag: no paint for longer than the wait limit.
    dummy._viewport_rotation_step_time -= 1.0
    SkyWindow._handle_client_key_press(dummy, _arrow_event(Qt.Key.Key_Right, auto_repeat=True))
    assert dummy.viewer_data.view_center == (20.0, 40.0)
