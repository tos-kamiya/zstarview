"""Bounded live cloud refresh acceleration based on observation age."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone

from .paths import CLOUD_UPDATE_INTERVAL

CLOUD_OLD_DATA_AGE = timedelta(minutes=30)
CLOUD_EARLY_REFRESH_SECONDS = 5 * 60


@dataclass
class CloudRefreshSchedule:
    last_observation: tuple[str, datetime] | None = None

    def next_delay_seconds(
        self, *, observed_at: datetime | None, completed_at: datetime, source: str
    ) -> int:
        """Shorten only for an old, newly seen observation from this source."""
        if observed_at is None:
            return CLOUD_UPDATE_INTERVAL
        observation = (
            observed_at.replace(tzinfo=timezone.utc)
            if observed_at.tzinfo is None
            else observed_at.astimezone(timezone.utc)
        )
        completed = (
            completed_at.replace(tzinfo=timezone.utc)
            if completed_at.tzinfo is None
            else completed_at.astimezone(timezone.utc)
        )
        if observation > completed:
            return CLOUD_UPDATE_INTERVAL
        key = (source, observation)
        repeated = key == self.last_observation
        self.last_observation = key
        if not repeated and completed - observation >= CLOUD_OLD_DATA_AGE:
            return CLOUD_EARLY_REFRESH_SECONDS
        return CLOUD_UPDATE_INTERVAL
