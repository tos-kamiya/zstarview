"""Qt-independent history records and minute selection contract."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path

UTC = timezone.utc
MINUTE = timedelta(minutes=1)


def parse_utc(value: object) -> datetime:
    result = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    if result.tzinfo is None:
        raise ValueError("Time must include a timezone")
    return result.astimezone(UTC)


def iso(value: datetime) -> str:
    return value.astimezone(UTC).isoformat().replace("+00:00", "Z")


def ceil_minute(value: datetime) -> datetime:
    rounded = value.replace(second=0, microsecond=0)
    return rounded if rounded == value else rounded + MINUTE


@dataclass(frozen=True)
class HistoryRecord:
    layer: str
    path: Path
    source_time: datetime
    fetched_at: datetime
    expires_at: datetime
    # Aircraft observations may have different contact times; do not merge
    # observations from different snapshots.
    contact_times: tuple[datetime, ...] = ()
    empty: bool = False
    digest: str = ""

    @property
    def starts_at(self) -> datetime:
        return (
            self.fetched_at
            if self.layer == "aircraft"
            else max(self.source_time, self.fetched_at)
        )

    def available(self, target: datetime) -> bool:
        if not self.starts_at <= target <= self.expires_at:
            return False
        if self.layer == "aircraft" and not self.empty:
            return any(
                contact <= target <= contact + timedelta(minutes=10)
                for contact in self.contact_times
            )
        return True
