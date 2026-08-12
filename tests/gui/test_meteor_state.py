from datetime import datetime, timezone

from zstarview.gui.meteor_state import MeteorState
from zstarview.meteors.types import MeteorWindowResult


def test_meteor_state_result_clears_banner() -> None:
    state = MeteorState(banner_text="GMN meteors: loading...")
    now = datetime(2026, 8, 12, tzinfo=timezone.utc)
    result = MeteorWindowResult(
        trails=(),
        window_start_utc=now,
        window_end_utc=now,
        source_files=(),
        unavailable_files=(),
    )
    state.set_result(result)
    assert state.result is result
    assert state.banner_text is None
