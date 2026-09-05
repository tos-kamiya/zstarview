from __future__ import annotations

from datetime import datetime, timezone

from zstarview.gui.application_services import ApplicationServices
from zstarview.gui.meteor_controller import MeteorController
from zstarview.meteors.types import MeteorWindowResult


def test_meteor_loader_runs_under_native_work_lock() -> None:
    services = ApplicationServices(max_workers=1)
    calls: list[bool] = []
    now = datetime(2026, 1, 1, tzinfo=timezone.utc)

    def loader(*args, **kwargs) -> MeteorWindowResult:
        calls.append(services.native_work_lock.locked())
        return MeteorWindowResult(
            trails=(),
            display_time_utc=now,
            window_start_utc=now,
            window_end_utc=now,
            source_files=(),
            unavailable_files=(),
        )

    controller = MeteorController(loader=loader, services=services)
    try:
        assert controller.update(
            display_time_utc=now,
            observer_lat=35.0,
            observer_lon=135.0,
            observer_height_m=0.0,
        )
        controller.shutdown(wait_timeout_s=1.0)
        assert calls == [True]
    finally:
        services.shutdown(wait=True)
