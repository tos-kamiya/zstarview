from __future__ import annotations

import urllib.error
from datetime import datetime, timezone
from email.message import Message

import pytest

from zstarview.location_resolver.place_cache import (
    build_place_cache_key,
    load_place_cache,
    save_place_cache,
)
from zstarview.location_resolver.place_search import (
    PlaceSearchCandidate,
    normalize_place_search_candidates,
    place_search_candidate_from_nominatim,
    search_place_candidates,
)


def test_place_search_candidate_from_nominatim_uses_namedetails_name() -> None:
    candidate = place_search_candidate_from_nominatim(
        {
            "lat": "35.681236",
            "lon": "139.767125",
            "display_name": "Tokyo Station, Chiyoda, Tokyo, Japan",
            "namedetails": {"name": "Tokyo Station"},
            "category": "railway",
            "type": "station",
            "importance": 0.9,
        }
    )

    assert candidate == PlaceSearchCandidate(
        name="Tokyo Station",
        display_name="Tokyo Station, Chiyoda, Tokyo, Japan",
        latitude_deg=35.681236,
        longitude_deg=139.767125,
        category="railway",
        type_name="station",
        importance=0.9,
    )


def test_place_search_candidate_from_nominatim_falls_back_to_display_name_prefix() -> None:
    candidate = place_search_candidate_from_nominatim(
        {
            "lat": "48.8582602",
            "lon": "2.2944991",
            "display_name": "Eiffel Tower, 5, Avenue Anatole France, Paris, France",
            "class": "man_made",
            "type": "tower",
            "importance": "0.75",
        }
    )

    assert candidate is not None
    assert candidate.name == "Eiffel Tower"
    assert candidate.category == "man_made"
    assert candidate.type_name == "tower"
    assert candidate.importance == 0.75


def test_normalize_place_search_candidates_filters_invalid_rows_and_sorts() -> None:
    candidates = normalize_place_search_candidates(
        [
            {
                "lat": "35.0",
                "lon": "139.0",
                "display_name": "Lower Importance, Tokyo, Japan",
                "name": "Lower Importance",
                "category": "place",
                "type": "city",
                "importance": 0.2,
            },
            {
                "lat": "35.1",
                "lon": "139.1",
                "display_name": "Higher Importance, Tokyo, Japan",
                "name": "Higher Importance",
                "category": "place",
                "type": "city",
                "importance": 0.8,
            },
            {
                "lat": "not-a-number",
                "lon": "139.2",
                "display_name": "Broken",
            },
        ]
    )

    assert [candidate.name for candidate in candidates] == ["Higher Importance", "Lower Importance"]


def test_search_place_candidates_uses_nominatim_helpers(monkeypatch, tmp_path) -> None:
    seen: dict[str, object] = {}

    def _build_url(query: str, *, limit: int, countrycode: str | None) -> str:
        seen["url"] = (query, limit, countrycode)
        return "https://example.invalid/search"

    def _fetch(url: str, *, language: str, user_agent: str):
        seen["fetch"] = (url, language, user_agent)
        return [
            {
                "lat": "35.681236",
                "lon": "139.767125",
                "display_name": "Tokyo Station, Chiyoda, Tokyo, Japan",
                "name": "Tokyo Station",
                "category": "railway",
                "type": "station",
                "importance": 0.9,
            }
        ]

    monkeypatch.setattr(
        "zstarview.location_resolver.place_search.nominatim._build_url",
        _build_url,
    )
    monkeypatch.setattr(
        "zstarview.location_resolver.place_search.nominatim._fetch",
        _fetch,
    )

    candidates = search_place_candidates(
        "Tokyo Station",
        limit=7,
        countrycode="jp",
        language="ja",
        user_agent="zstarview-test/1.0",
        cache_root=tmp_path,
    )

    assert seen["url"] == ("Tokyo Station", 7, "jp")
    assert seen["fetch"] == ("https://example.invalid/search", "ja", "zstarview-test/1.0")
    assert candidates[0].name == "Tokyo Station"


def test_search_place_candidates_reports_cache_miss_for_transport_errors(
    monkeypatch, tmp_path
) -> None:
    monkeypatch.setattr(
        "zstarview.location_resolver.place_search.nominatim._build_url",
        lambda query, *, limit, countrycode: "https://example.invalid/search",
    )
    monkeypatch.setattr(
        "zstarview.location_resolver.place_search.nominatim._fetch",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(urllib.error.URLError("offline")),
    )

    with pytest.raises(RuntimeError, match="no cache is available"):
        search_place_candidates("Tokyo Station", cache_root=tmp_path)


def test_search_place_candidates_uses_cache_only_after_network_failure(
    monkeypatch, tmp_path
) -> None:
    fetched_at = datetime(2026, 8, 1, 2, 3, 4, tzinfo=timezone.utc)
    key = build_place_cache_key("Matsue Station", "jp", "en")
    save_place_cache(
        key,
        (
            {
                "category": "railway",
                "display_name": "Matsue Station, Shimane, Japan",
                "importance": 0.9,
                "lat": 35.464,
                "lon": 133.063,
                "name": "Matsue Station",
                "type": "station",
            },
        ),
        original_query="Matsue Station",
        fetched_at_utc=fetched_at,
        cache_root=tmp_path,
    )
    monkeypatch.setattr(
        "zstarview.location_resolver.place_search.nominatim._fetch",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            urllib.error.URLError("offline")
        ),
    )

    candidates = search_place_candidates(
        "Matsue Station", countrycode="jp", cache_root=tmp_path
    )

    assert candidates[0].name == "Matsue Station"
    assert candidates[0].cache_fetched_at_utc == fetched_at


def test_online_success_replaces_cache_without_reading_it_first(
    monkeypatch, tmp_path
) -> None:
    key = build_place_cache_key("Matsue Station", None, "en")
    monkeypatch.setattr(
        "zstarview.location_resolver.place_search.load_place_cache",
        lambda *_args, **_kwargs: pytest.fail("online success must not read cache"),
    )
    monkeypatch.setattr(
        "zstarview.location_resolver.place_search.nominatim._fetch",
        lambda *_args, **_kwargs: [
            {
                "lat": "35.464",
                "lon": "133.063",
                "display_name": "Fresh Matsue Station",
                "name": "Fresh Matsue Station",
                "category": "railway",
                "type": "station",
                "importance": 1.0,
            }
        ],
    )

    candidates = search_place_candidates("Matsue Station", cache_root=tmp_path)

    assert candidates[0].name == "Fresh Matsue Station"
    cached = load_place_cache(key, cache_root=tmp_path)
    assert cached is not None
    assert cached.results[0]["name"] == "Fresh Matsue Station"


def test_empty_online_result_does_not_replace_existing_cache(monkeypatch, tmp_path) -> None:
    fetched_at = datetime(2026, 8, 1, tzinfo=timezone.utc)
    key = build_place_cache_key("Matsue Station", None, "en")
    original = (
        {
            "category": "railway",
            "display_name": "Cached Matsue Station",
            "importance": 0.9,
            "lat": 35.464,
            "lon": 133.063,
            "name": "Cached Matsue Station",
            "type": "station",
        },
    )
    save_place_cache(
        key,
        original,
        original_query="Matsue Station",
        fetched_at_utc=fetched_at,
        cache_root=tmp_path,
    )
    monkeypatch.setattr(
        "zstarview.location_resolver.place_search.nominatim._fetch",
        lambda *_args, **_kwargs: [],
    )

    assert search_place_candidates("Matsue Station", cache_root=tmp_path) == ()
    cached = load_place_cache(key, cache_root=tmp_path)
    assert cached is not None
    assert cached.fetched_at_utc == fetched_at


def test_http_error_does_not_fall_back_to_cache(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(
        "zstarview.location_resolver.place_search.nominatim._fetch",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            urllib.error.HTTPError(
                "https://example.invalid", 429, "rate limited", Message(), None
            )
        ),
    )

    with pytest.raises(urllib.error.HTTPError):
        search_place_candidates("Matsue Station", cache_root=tmp_path)


def test_request_slot_waits_until_one_second_after_previous_start(tmp_path) -> None:
    marker = tmp_path / "last_request.json"
    marker.write_text('{"last_request_started_at": 100.0}', encoding="ascii")
    times = iter((100.2, 101.0))
    sleeps: list[float] = []

    from zstarview.location_resolver.place_search import _nominatim_request_slot

    with _nominatim_request_slot(
        tmp_path,
        now_func=lambda: next(times),
        sleep_func=sleeps.append,
    ):
        pass

    assert sleeps == pytest.approx([0.8])
