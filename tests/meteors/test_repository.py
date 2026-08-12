from __future__ import annotations

from datetime import datetime, timedelta, timezone

from zstarview.meteors.repository import GmnMeteorRepository


def _index_html(*filenames: str) -> str:
    return "\n".join(f'<a href="{name}">{name}</a>' for name in filenames)


FILENAMES = (
    "traj_summary_20260809_solrange_137.0-138.0.txt",
    "traj_summary_20260810_solrange_138.0-139.0.txt",
    "traj_summary_20260811_solrange_139.0-140.0.txt",
    "traj_summary_20260812_solrange_140.0-141.0.txt",
)


def test_repository_finds_latest_observation_then_loads_its_24_hour_window(
    tmp_path,
    summary_row_factory,
) -> None:
    index_url = "https://example.test/daily/"
    payloads = {
        index_url: _index_html(*FILENAMES),
        index_url + FILENAMES[0]: summary_row_factory(
            trajectory_id="inside-old",
            beginning_utc="2026-08-09 18:30:00.000000",
        ),
        index_url + FILENAMES[1]: summary_row_factory(
            trajectory_id="latest",
            beginning_utc="2026-08-10 18:00:00.000000",
        ),
        index_url + FILENAMES[2]: summary_row_factory(
            trajectory_id="after-display-time",
            beginning_utc="2026-08-11 18:00:00.000000",
        ),
        index_url + FILENAMES[3]: "# no completed trajectories yet\n",
    }
    calls: list[str] = []

    def fetcher(url: str, *, timeout_s: float) -> str:
        calls.append(url)
        return payloads[url]

    repository = GmnMeteorRepository(
        cache_root=tmp_path,
        index_url=index_url,
        fetcher=fetcher,
    )
    display_time = datetime(2026, 8, 11, 12, tzinfo=timezone.utc)
    result = repository.load_latest_window(
        display_time,
        now_utc=display_time,
    )

    assert result.window_end_utc == datetime(
        2026, 8, 10, 18, tzinfo=timezone.utc
    )
    assert [item.trajectory_id for item in result.observations] == [
        "inside-old",
        "latest",
    ]
    assert index_url + FILENAMES[3] in calls
    assert index_url + FILENAMES[2] in calls
    assert index_url + FILENAMES[1] in calls
    assert index_url + FILENAMES[0] in calls


def test_repository_discovers_adjacent_solar_days_and_filters_exact_window(
    tmp_path,
    summary_row_factory,
) -> None:
    index_url = "https://example.test/daily/"
    payloads = {
        index_url: _index_html(*FILENAMES),
        **{
            index_url + filename: summary_row_factory(
                trajectory_id=f"meteor-{offset}",
                beginning_utc=f"2026-08-{9 + offset:02d} 18:00:00.000000",
            )
            for offset, filename in enumerate(FILENAMES)
        },
    }
    calls: list[str] = []

    def fetcher(url: str, *, timeout_s: float) -> str:
        calls.append(url)
        return payloads[url]

    repository = GmnMeteorRepository(
        cache_root=tmp_path,
        index_url=index_url,
        fetcher=fetcher,
    )
    result = repository.load_window(
        datetime(2026, 8, 10, 12, tzinfo=timezone.utc),
        datetime(2026, 8, 11, 12, tzinfo=timezone.utc),
        now_utc=datetime(2026, 8, 12, 12, tzinfo=timezone.utc),
    )

    assert [item.trajectory_id for item in result.observations] == ["meteor-1"]
    assert result.source_files == FILENAMES
    assert not result.unavailable_files
    assert calls[0] == index_url


def test_repository_reuses_fresh_cache_without_network(
    tmp_path,
    summary_row_factory,
) -> None:
    index_url = "https://example.test/daily/"
    payloads = {
        index_url: _index_html(*FILENAMES),
        **{
            index_url + filename: summary_row_factory() + "\r\n"
            for filename in FILENAMES
        },
    }

    def fetcher(url: str, *, timeout_s: float) -> str:
        return payloads[url]

    repository = GmnMeteorRepository(
        cache_root=tmp_path,
        index_url=index_url,
        fetcher=fetcher,
    )
    now = datetime(2026, 8, 12, 12, tzinfo=timezone.utc)
    repository.load_window(now - timedelta(hours=24), now, now_utc=now)

    def fail_fetcher(url: str, *, timeout_s: float) -> str:
        raise AssertionError(f"network should not be used: {url}")

    cached_repository = GmnMeteorRepository(
        cache_root=tmp_path,
        index_url=index_url,
        fetcher=fail_fetcher,
    )
    result = cached_repository.load_window(
        now - timedelta(hours=24),
        now,
        now_utc=now + timedelta(hours=1),
    )

    assert result.source_files == FILENAMES[1:]
    assert not result.used_stale_index
    assert not result.used_stale_files


def test_repository_uses_stale_cache_when_refresh_fails(
    tmp_path,
    summary_row_factory,
) -> None:
    index_url = "https://example.test/daily/"
    payloads = {
        index_url: _index_html(*FILENAMES),
        **{index_url + filename: summary_row_factory() for filename in FILENAMES},
    }

    def fetcher(url: str, *, timeout_s: float) -> str:
        return payloads[url]

    now = datetime(2026, 8, 12, 12, tzinfo=timezone.utc)
    repository = GmnMeteorRepository(
        cache_root=tmp_path,
        index_url=index_url,
        fetcher=fetcher,
    )
    repository.load_window(now - timedelta(hours=24), now, now_utc=now)

    def failing_fetcher(url: str, *, timeout_s: float) -> str:
        raise OSError("offline")

    offline = GmnMeteorRepository(
        cache_root=tmp_path,
        index_url=index_url,
        fetcher=failing_fetcher,
    )
    result = offline.load_window(
        now - timedelta(hours=24),
        now,
        now_utc=now + timedelta(hours=7),
    )

    assert result.used_stale_index
    assert result.used_stale_files
    assert not result.unavailable_files


def test_repository_reports_missing_files_but_returns_partial_data(
    tmp_path,
    summary_row_factory,
) -> None:
    index_url = "https://example.test/daily/"

    def fetcher(url: str, *, timeout_s: float) -> str:
        if url == index_url:
            return _index_html(*FILENAMES)
        if url.endswith(FILENAMES[1]):
            raise OSError("missing")
        return summary_row_factory(trajectory_id=url.rsplit("/", 1)[-1])

    repository = GmnMeteorRepository(
        cache_root=tmp_path,
        index_url=index_url,
        fetcher=fetcher,
    )
    result = repository.load_window(
        datetime(2026, 8, 10, 0, tzinfo=timezone.utc),
        datetime(2026, 8, 11, 0, tzinfo=timezone.utc),
        now_utc=datetime(2026, 8, 12, 0, tzinfo=timezone.utc),
    )

    assert result.unavailable_files == (FILENAMES[1],)
    assert len(result.source_files) == 3


def test_repository_replaces_cache_when_content_hash_is_invalid(
    tmp_path,
    summary_row_factory,
) -> None:
    index_url = "https://example.test/daily/"
    payloads = {
        index_url: _index_html(*FILENAMES),
        **{index_url + filename: summary_row_factory() for filename in FILENAMES},
    }
    calls: list[str] = []

    def fetcher(url: str, *, timeout_s: float) -> str:
        calls.append(url)
        return payloads[url]

    now = datetime(2026, 8, 12, 12, tzinfo=timezone.utc)
    repository = GmnMeteorRepository(
        cache_root=tmp_path,
        index_url=index_url,
        fetcher=fetcher,
    )
    repository.load_window(now - timedelta(hours=24), now, now_utc=now)
    damaged_path = tmp_path / "daily" / FILENAMES[1]
    damaged_path.write_text("damaged", encoding="utf-8")
    calls.clear()

    repository.load_window(
        now - timedelta(hours=24),
        now,
        now_utc=now + timedelta(hours=1),
    )

    assert index_url + FILENAMES[1] in calls
    assert damaged_path.read_text(encoding="utf-8") != "damaged"
