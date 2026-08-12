from __future__ import annotations

from datetime import datetime, timedelta, timezone

from zstarview.meteors.service import load_celestial_meteor_trails
from zstarview.meteors.types import GmnLoadResult, MeteorObservation


class _StubRepository:
    def __init__(self, observations: tuple[MeteorObservation, ...]) -> None:
        self.observations = observations
        self.requested_display_time = None

    def load_latest_window(self, display_time_utc, *, now_utc=None):
        self.requested_display_time = (display_time_utc, now_utc)
        window_end = max(item.beginning_utc for item in self.observations)
        return GmnLoadResult(
            observations=self.observations,
            source_files=("sample.txt",),
            unavailable_files=("missing.txt",),
            used_stale_index=True,
            used_stale_files=True,
            window_end_utc=window_end,
        )


def test_service_anchors_24_hour_window_to_latest_observation() -> None:
    display_time = datetime(2026, 8, 10, 12, tzinfo=timezone.utc)
    observation = MeteorObservation(
        trajectory_id="overhead",
        beginning_utc=display_time - timedelta(hours=1),
        begin_lat_deg=35.0,
        begin_lon_deg=139.0,
        begin_height_km=100.0,
        end_lat_deg=35.0,
        end_lon_deg=139.0,
        end_height_km=80.0,
    )
    repository = _StubRepository((observation,))

    result = load_celestial_meteor_trails(
        display_time,
        observer_lat=35.0,
        observer_lon=139.0,
        observer_height_m=0.0,
        repository=repository,
        now_utc=display_time,
    )

    assert result.display_time_utc == display_time
    assert result.window_start_utc == display_time - timedelta(hours=25)
    assert result.window_end_utc == display_time - timedelta(hours=1)
    assert len(result.trails) == 1
    assert result.source_files == ("sample.txt",)
    assert result.unavailable_files == ("missing.txt",)
    assert result.used_stale_index
    assert result.used_stale_files
    assert repository.requested_display_time == (display_time, display_time)
