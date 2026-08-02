from __future__ import annotations

import json
from datetime import datetime, timezone

from zstarview.location_resolver import place_cache


def _result(name: str, latitude: float = 35.0) -> tuple[dict[str, object], ...]:
    return (
        {
            "category": "railway",
            "display_name": name,
            "importance": 0.8,
            "lat": latitude,
            "lon": 133.0,
            "name": name,
            "type": "station",
        },
    )


def test_cache_key_normalizes_unicode_whitespace_and_countrycode() -> None:
    first = place_cache.build_place_cache_key("  松江駅　北口 ", " JP ", "ja")
    second = place_cache.build_place_cache_key("松江駅 北口", "jp", "ja")

    assert first == second
    assert first.query == "松江駅 北口"


def test_cache_filename_is_fixed_ascii_for_unsafe_query(tmp_path) -> None:
    key = place_cache.build_place_cache_key("../CON\\駅:*?🚉", None, "ja")
    path = place_cache.place_cache_path(key, cache_root=tmp_path)

    assert path.parent == tmp_path
    assert len(path.stem) == 64
    assert set(path.stem) <= set("0123456789abcdef")


def test_cache_round_trip_and_matching_entry_update(tmp_path) -> None:
    key = place_cache.build_place_cache_key("松江駅", "jp", "ja")
    first_time = datetime(2026, 8, 1, 1, 2, 3, tzinfo=timezone.utc)
    second_time = datetime(2026, 8, 2, 4, 5, 6, tzinfo=timezone.utc)

    place_cache.save_place_cache(
        key,
        _result("old"),
        original_query="松江駅",
        fetched_at_utc=first_time,
        cache_root=tmp_path,
    )
    place_cache.save_place_cache(
        key,
        _result("new"),
        original_query=" 松江駅 ",
        fetched_at_utc=second_time,
        cache_root=tmp_path,
    )

    cached = place_cache.load_place_cache(key, cache_root=tmp_path)
    assert cached is not None
    assert cached.fetched_at_utc == second_time
    assert cached.results[0]["name"] == "new"
    payload = json.loads(place_cache.place_cache_path(key, cache_root=tmp_path).read_text())
    assert len(payload["entries"]) == 1


def test_digest_collision_keeps_both_complete_keys(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(place_cache, "place_cache_digest", lambda _key: "a" * 64)
    first = place_cache.build_place_cache_key("First", None, "en")
    second = place_cache.build_place_cache_key("Second", None, "en")
    fetched = datetime(2026, 8, 2, tzinfo=timezone.utc)

    place_cache.save_place_cache(
        first, _result("First"), original_query="First", fetched_at_utc=fetched, cache_root=tmp_path
    )
    place_cache.save_place_cache(
        second,
        _result("Second"),
        original_query="Second",
        fetched_at_utc=fetched,
        cache_root=tmp_path,
    )

    first_cached = place_cache.load_place_cache(first, cache_root=tmp_path)
    second_cached = place_cache.load_place_cache(second, cache_root=tmp_path)
    assert first_cached is not None
    assert second_cached is not None
    assert first_cached.results[0]["name"] == "First"
    assert second_cached.results[0]["name"] == "Second"
    payload = json.loads(place_cache.place_cache_path(first, cache_root=tmp_path).read_text())
    assert len(payload["entries"]) == 2


def test_invalid_bucket_is_not_overwritten(tmp_path) -> None:
    key = place_cache.build_place_cache_key("Matsue", None, "en")
    path = place_cache.place_cache_path(key, cache_root=tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("not-json", encoding="utf-8")

    saved = place_cache.save_place_cache(
        key,
        _result("Matsue"),
        original_query="Matsue",
        fetched_at_utc=datetime.now(timezone.utc),
        cache_root=tmp_path,
    )

    assert saved is None
    assert path.read_text(encoding="utf-8") == "not-json"
