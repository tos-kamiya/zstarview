from types import SimpleNamespace

from zstarview.gui.window_actions import SkyWindowActionsMixin


class _CheckedAction:
    def __init__(self) -> None:
        self.checked = True

    def setChecked(self, checked: bool) -> None:
        self.checked = checked


def test_toggle_twinkle_stops_and_restarts_without_changing_count() -> None:
    updates: list[str] = []
    owner = SimpleNamespace(
        twinkle_count=30,
        twinkle_enabled=True,
        state=SimpleNamespace(twinkle_targets=((4, 0.5),), twinkle_bucket=12),
        _action_toggle_twinkle=_CheckedAction(),
        request_client_update=lambda: updates.append("paint"),
        _update_twinkle=lambda: updates.append("sample"),
    )

    SkyWindowActionsMixin.toggle_twinkle(owner)

    assert owner.twinkle_enabled is False
    assert owner.twinkle_count == 30
    assert owner.state.twinkle_targets == ()
    assert owner.state.twinkle_bucket is None
    assert owner._action_toggle_twinkle.checked is False
    assert updates == ["paint"]

    SkyWindowActionsMixin.toggle_twinkle(owner)

    assert owner.twinkle_enabled is True
    assert owner._action_toggle_twinkle.checked is True
    assert updates == ["paint", "sample"]


def test_toggle_twinkle_is_locked_out_when_count_is_zero() -> None:
    owner = SimpleNamespace(twinkle_count=0, twinkle_enabled=False)

    SkyWindowActionsMixin.toggle_twinkle(owner)

    assert owner.twinkle_enabled is False
