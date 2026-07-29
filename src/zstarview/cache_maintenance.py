from __future__ import annotations

import json
import logging
import shutil
from datetime import datetime, timedelta, timezone
from pathlib import Path

from .paths import (
    CACHE_PATH,
    OVERTURE_DERIVED_ROOT_DIR,
    OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR,
    PLATEAU_DERIVED_ROOT_DIR,
)

logger = logging.getLogger(__name__)

CLEAR_LONG_LIVED_CACHE_METADATA = Path(CACHE_PATH) / "clear_long_lived_cache_meta.json"
CLEAR_LONG_LIVED_CACHE_COOLDOWN = timedelta(days=3)


class LongLivedCacheClearCooldownError(RuntimeError):
    def __init__(
        self, *, last_cleared_at_utc: datetime, retry_at_utc: datetime, message: str
    ) -> None:
        super().__init__(message)
        self.last_cleared_at_utc = last_cleared_at_utc
        self.retry_at_utc = retry_at_utc


def long_lived_cache_targets() -> tuple[Path, ...]:
    return (
        Path(CACHE_PATH) / "copernicus-dem",
        Path(PLATEAU_DERIVED_ROOT_DIR),
        Path(OVERTURE_DERIVED_ROOT_DIR),
        Path(OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR),
    )


def format_manual_long_lived_cache_clear_hint() -> str:
    targets = ", ".join(str(path) for path in long_lived_cache_targets())
    return f"If you really need to clear it again, remove these directories manually: {targets}"


def _normalize_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _parse_optional_utc(text: object) -> datetime | None:
    if not isinstance(text, str) or not text.strip():
        return None
    return _normalize_utc(datetime.fromisoformat(text.replace("Z", "+00:00")))


def _read_last_cleared_at_utc(metadata_path: Path) -> datetime | None:
    if not metadata_path.exists():
        return None
    try:
        payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if not isinstance(payload, dict):
        return None
    try:
        return _parse_optional_utc(payload.get("last_cleared_at_utc"))
    except Exception:
        return None


def _write_last_cleared_at_utc(
    metadata_path: Path, *, cleared_at_utc: datetime
) -> None:
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(
        json.dumps(
            {
                "last_cleared_at_utc": _normalize_utc(cleared_at_utc).isoformat(),
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )


def enforce_long_lived_cache_clear_cooldown(
    *,
    now_utc: datetime | None = None,
    metadata_path: Path | None = None,
) -> None:
    if metadata_path is None:
        metadata_path = CLEAR_LONG_LIVED_CACHE_METADATA
    now = _normalize_utc(now_utc or datetime.now(timezone.utc))
    last_cleared_at_utc = _read_last_cleared_at_utc(metadata_path)
    if last_cleared_at_utc is None:
        return
    retry_at_utc = last_cleared_at_utc + CLEAR_LONG_LIVED_CACHE_COOLDOWN
    if now >= retry_at_utc:
        return
    message = (
        "Long-lived cache was already cleared on "
        f"{last_cleared_at_utc.strftime('%Y-%m-%d %H:%M UTC')}. "
        f"Retry is allowed after {retry_at_utc.strftime('%Y-%m-%d %H:%M UTC')}. "
        "This option is for troubleshooting, not normal startup. "
        + format_manual_long_lived_cache_clear_hint()
    )
    raise LongLivedCacheClearCooldownError(
        last_cleared_at_utc=last_cleared_at_utc,
        retry_at_utc=retry_at_utc,
        message=message,
    )


def clear_long_lived_cache(
    *,
    now_utc: datetime | None = None,
    metadata_path: Path | None = None,
) -> tuple[Path, ...]:
    if metadata_path is None:
        metadata_path = CLEAR_LONG_LIVED_CACHE_METADATA
    now = _normalize_utc(now_utc or datetime.now(timezone.utc))
    enforce_long_lived_cache_clear_cooldown(now_utc=now, metadata_path=metadata_path)
    removed: list[Path] = []
    for path in long_lived_cache_targets():
        if not path.exists():
            logger.info("Long-lived cache already absent: %s", path)
            continue
        shutil.rmtree(path)
        removed.append(path)
        logger.info("Cleared long-lived cache: %s", path)
    _write_last_cleared_at_utc(metadata_path, cleared_at_utc=now)
    return tuple(removed)
