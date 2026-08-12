"""GMN daily-file discovery, fetching, caching, and window loading."""

from __future__ import annotations

import hashlib
import json
import tempfile
from collections.abc import Callable
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

from ..paths import GMN_METEOR_CACHE_ROOT_DIR
from .client import fetch_gmn_text, gmn_daily_file_url
from .constants import (
    GMN_CACHE_SCHEMA,
    GMN_DAILY_INDEX_URL,
    GMN_INDEX_FRESH_TTL,
    GMN_LATEST_LOOKBACK,
    GMN_RECENT_FILE_FRESH_TTL,
    GMN_WINDOW,
)
from .parser import parse_gmn_daily_index, parse_gmn_trajectory_summary
from .types import GmnDailyFile, GmnLoadResult, MeteorObservation

TextFetcher = Callable[..., str]


class GmnMeteorRepository:
    def __init__(
        self,
        *,
        cache_root: str | Path = GMN_METEOR_CACHE_ROOT_DIR,
        index_url: str = GMN_DAILY_INDEX_URL,
        fetcher: TextFetcher = fetch_gmn_text,
        timeout_s: float = 30.0,
    ) -> None:
        self.cache_root = Path(cache_root)
        self.index_url = str(index_url)
        self.fetcher = fetcher
        self.timeout_s = float(timeout_s)

    def load_window(
        self,
        window_start_utc: datetime,
        window_end_utc: datetime,
        *,
        now_utc: datetime | None = None,
    ) -> GmnLoadResult:
        start = _normalize_utc(window_start_utc)
        end = _normalize_utc(window_end_utc)
        if start > end:
            raise ValueError("GMN window start must not be after end")
        now = _normalize_utc(now_utc or datetime.now(timezone.utc))
        daily_files, stale_index = self._load_index(now)
        candidates = _select_candidate_files(daily_files, start.date(), end.date())

        observations: dict[str, MeteorObservation] = {}
        source_files: list[str] = []
        unavailable_files: list[str] = []
        used_stale_files = False
        for daily_file in candidates:
            try:
                text, stale = self._load_daily_file(daily_file, now)
            except Exception:
                unavailable_files.append(daily_file.filename)
                continue
            source_files.append(daily_file.filename)
            used_stale_files = used_stale_files or stale
            for observation in parse_gmn_trajectory_summary(text):
                if start <= observation.beginning_utc <= end:
                    observations[observation.trajectory_id] = observation

        ordered = tuple(
            sorted(
                observations.values(),
                key=lambda item: (item.beginning_utc, item.trajectory_id),
            )
        )
        return GmnLoadResult(
            observations=ordered,
            source_files=tuple(source_files),
            unavailable_files=tuple(unavailable_files),
            used_stale_index=stale_index,
            used_stale_files=used_stale_files,
        )

    def load_latest_window(
        self,
        display_time_utc: datetime,
        *,
        now_utc: datetime | None = None,
        window: timedelta = GMN_WINDOW,
        lookback: timedelta = GMN_LATEST_LOOKBACK,
    ) -> GmnLoadResult:
        display_time = _normalize_utc(display_time_utc)
        now = _normalize_utc(now_utc or datetime.now(timezone.utc))
        if window <= timedelta(0):
            raise ValueError("GMN window must be positive")
        if lookback < window:
            raise ValueError("GMN latest-window lookback must cover the window")

        daily_files, stale_index = self._load_index(now)
        search_start = display_time - lookback
        search_candidates = _select_candidate_files(
            daily_files,
            search_start.date(),
            display_time.date(),
        )
        loaded_text: dict[str, tuple[str, bool] | None] = {}
        observations: dict[str, MeteorObservation] = {}
        source_files: list[str] = []
        unavailable_files: list[str] = []
        used_stale_files = False

        def load_candidate(daily_file: GmnDailyFile) -> None:
            nonlocal used_stale_files
            if daily_file.filename in loaded_text:
                return
            try:
                text, stale = self._load_daily_file(daily_file, now)
            except Exception:
                loaded_text[daily_file.filename] = None
                unavailable_files.append(daily_file.filename)
                return
            loaded_text[daily_file.filename] = (text, stale)
            source_files.append(daily_file.filename)
            used_stale_files = used_stale_files or stale
            for observation in parse_gmn_trajectory_summary(text):
                if search_start <= observation.beginning_utc <= display_time:
                    observations[observation.trajectory_id] = observation

        latest: datetime | None = None
        for daily_file in reversed(search_candidates):
            load_candidate(daily_file)
            if observations:
                latest = max(item.beginning_utc for item in observations.values())
                break

        if latest is not None:
            window_start = latest - window
            for daily_file in _select_candidate_files(
                daily_files,
                window_start.date(),
                latest.date(),
            ):
                load_candidate(daily_file)
            observations = {
                key: item
                for key, item in observations.items()
                if window_start <= item.beginning_utc <= latest
            }

        ordered = tuple(
            sorted(
                observations.values(),
                key=lambda item: (item.beginning_utc, item.trajectory_id),
            )
        )
        return GmnLoadResult(
            observations=ordered,
            source_files=tuple(source_files),
            unavailable_files=tuple(unavailable_files),
            used_stale_index=stale_index,
            used_stale_files=used_stale_files,
            window_end_utc=latest,
        )

    def _load_index(self, now_utc: datetime) -> tuple[tuple[GmnDailyFile, ...], bool]:
        path = self.cache_root / "daily-index.json"
        cached = _read_json(path)
        if cached is not None and _cache_is_fresh(cached, now_utc, GMN_INDEX_FRESH_TTL):
            return _daily_files_from_cache(cached), False
        try:
            text = self.fetcher(self.index_url, timeout_s=self.timeout_s)
            files = parse_gmn_daily_index(text)
            if not files:
                raise ValueError("GMN daily index contains no trajectory files")
            payload: dict[str, object] = {
                "schema": GMN_CACHE_SCHEMA,
                "fetched_at_utc": now_utc.isoformat(),
                "files": [
                    {"filename": item.filename, "nominal_date": item.nominal_date.isoformat()}
                    for item in files
                ],
            }
            _write_json(path, payload)
            return files, False
        except Exception:
            if cached is None:
                raise
            return _daily_files_from_cache(cached), True

    def _load_daily_file(
        self,
        daily_file: GmnDailyFile,
        now_utc: datetime,
    ) -> tuple[str, bool]:
        data_path = self.cache_root / "daily" / daily_file.filename
        meta_path = data_path.with_suffix(data_path.suffix + ".json")
        cached_meta = _read_json(meta_path)
        cached_text = _read_cached_daily_text(
            data_path,
            cached_meta,
            expected_filename=daily_file.filename,
        )
        has_cache = cached_text is not None
        historical = daily_file.nominal_date < now_utc.date() - timedelta(days=2)
        if cached_text is not None and cached_meta is not None and (
            historical
            or _cache_is_fresh(cached_meta, now_utc, GMN_RECENT_FILE_FRESH_TTL)
        ):
            return cached_text, False
        try:
            text = self.fetcher(
                gmn_daily_file_url(daily_file.filename, index_url=self.index_url),
                timeout_s=self.timeout_s,
            )
            # Parse before replacing a usable cache, so an HTML error page or a
            # truncated response cannot become the authoritative local copy.
            parse_result = parse_gmn_trajectory_summary(text)
            if not parse_result and not any(line.startswith("#") for line in text.splitlines()):
                raise ValueError("GMN daily file is not a trajectory summary")
            _write_text(data_path, text)
            _write_json(
                meta_path,
                {
                    "schema": GMN_CACHE_SCHEMA,
                    "filename": daily_file.filename,
                    "nominal_date": daily_file.nominal_date.isoformat(),
                    "fetched_at_utc": now_utc.isoformat(),
                    "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
                },
            )
            return text, False
        except Exception:
            if not has_cache or cached_text is None:
                raise
            return cached_text, True


def _select_candidate_files(
    files: tuple[GmnDailyFile, ...],
    start_date: date,
    end_date: date,
) -> tuple[GmnDailyFile, ...]:
    # GMN daily summaries cover one degree of solar longitude rather than a
    # midnight-to-midnight UTC day. Adjacent nominal dates are therefore loaded
    # and the record timestamps provide the exact window filter.
    lower = start_date - timedelta(days=1)
    upper = end_date + timedelta(days=1)
    return tuple(item for item in files if lower <= item.nominal_date <= upper)


def _normalize_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _cache_is_fresh(payload: dict[str, object], now_utc: datetime, ttl: timedelta) -> bool:
    try:
        fetched_at = _normalize_utc(datetime.fromisoformat(str(payload["fetched_at_utc"])))
    except (KeyError, TypeError, ValueError):
        return False
    age = now_utc - fetched_at
    return timedelta(0) <= age <= ttl


def _daily_files_from_cache(payload: dict[str, object]) -> tuple[GmnDailyFile, ...]:
    if payload.get("schema") != GMN_CACHE_SCHEMA:
        raise ValueError("unsupported GMN cache schema")
    raw_files = payload.get("files")
    if not isinstance(raw_files, list):
        raise ValueError("GMN index cache is missing files")
    files: list[GmnDailyFile] = []
    for item in raw_files:
        if not isinstance(item, dict):
            raise ValueError("invalid GMN index cache entry")
        files.append(
            GmnDailyFile(
                filename=str(item["filename"]),
                nominal_date=date.fromisoformat(str(item["nominal_date"])),
            )
        )
    return tuple(files)


def _read_json(path: Path) -> dict[str, object] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def _read_cached_daily_text(
    path: Path,
    metadata: dict[str, object] | None,
    *,
    expected_filename: str,
) -> str | None:
    if not path.is_file() or metadata is None:
        return None
    if metadata.get("schema") != GMN_CACHE_SCHEMA:
        return None
    if metadata.get("filename") != expected_filename:
        return None
    expected_hash = metadata.get("sha256")
    if not isinstance(expected_hash, str):
        return None
    try:
        # GMN summaries currently use CRLF. Preserve the original line endings
        # because the metadata hash covers the downloaded UTF-8 text exactly.
        with path.open("r", encoding="utf-8", newline="") as handle:
            text = handle.read()
    except OSError:
        return None
    actual_hash = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return text if actual_hash == expected_hash else None


def _write_json(path: Path, payload: dict[str, object]) -> None:
    _write_text(path, json.dumps(payload, separators=(",", ":"), sort_keys=True))


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=path.name + ".",
        suffix=".tmp",
        dir=path.parent,
        delete=False,
    ) as handle:
        handle.write(text)
        temporary = Path(handle.name)
    temporary.replace(path)
