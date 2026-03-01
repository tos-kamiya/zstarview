from __future__ import annotations

from zstarview.render.draw import effective_star_draw_vmag_limit


def test_effective_star_draw_vmag_limit_uses_base_when_not_interacting() -> None:
    assert effective_star_draw_vmag_limit(10.0, False, 7.0) == 10.0
    assert effective_star_draw_vmag_limit(6.0, False, 7.0) == 6.0


def test_effective_star_draw_vmag_limit_caps_when_interacting() -> None:
    assert effective_star_draw_vmag_limit(10.0, True, 7.0) == 7.0
    assert effective_star_draw_vmag_limit(8.2, True, 7.0) == 7.0
    assert effective_star_draw_vmag_limit(6.5, True, 7.0) == 6.5
