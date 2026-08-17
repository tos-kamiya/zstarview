from __future__ import annotations

from zstarview.gui.window_actions import SkyWindowActionsMixin


class _Action:
    def __init__(self, checked: bool) -> None:
        self._checked = checked

    def isChecked(self) -> bool:
        return self._checked

    def setChecked(self, checked: bool) -> None:
        self._checked = bool(checked)


class _Window:
    _moon_toggle_target = SkyWindowActionsMixin._moon_toggle_target
    _moon_toggle_active = SkyWindowActionsMixin._moon_toggle_active
    toggle_moon_option = SkyWindowActionsMixin.toggle_moon_option

    def __init__(self, style: str, scale: int) -> None:
        self._configured_moon_style = style
        self._configured_moon_scale = scale
        self.moon_style = style
        self.moon_scale = scale
        self._action_moon_option = _Action(self._moon_toggle_active())
        self.update_requests = 0

    def request_client_update(self) -> None:
        self.update_requests += 1


def test_default_moon_option_toggles_temporary_sphere_five() -> None:
    window = _Window("marker", 1)

    window.toggle_moon_option()
    assert (window.moon_style, window.moon_scale) == ("sphere", 5)
    assert window._action_moon_option.isChecked() is True

    window.toggle_moon_option()
    assert (window.moon_style, window.moon_scale) == ("marker", 1)
    assert window._action_moon_option.isChecked() is False
    assert window.update_requests == 2


def test_nondefault_moon_option_toggles_temporary_off() -> None:
    window = _Window("image", 4)

    window.toggle_moon_option()
    assert (window.moon_style, window.moon_scale) == ("marker", 1)
    assert window._action_moon_option.isChecked() is False

    window.toggle_moon_option()
    assert (window.moon_style, window.moon_scale) == ("image", 4)
    assert window._action_moon_option.isChecked() is True


def test_nondefault_marker_scale_is_restored_after_toggle() -> None:
    window = _Window("marker", 3)

    window.toggle_moon_option()
    assert (window.moon_style, window.moon_scale) == ("marker", 1)

    window.toggle_moon_option()
    assert (window.moon_style, window.moon_scale) == ("marker", 3)
