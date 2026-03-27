from __future__ import annotations

import urllib.error

import pytest

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


def test_search_place_candidates_uses_nominatim_helpers(monkeypatch) -> None:
    seen = {}

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
    )

    assert seen["url"] == ("Tokyo Station", 7, "jp")
    assert seen["fetch"] == ("https://example.invalid/search", "ja", "zstarview-test/1.0")
    assert candidates[0].name == "Tokyo Station"


def test_search_place_candidates_propagates_transport_errors(monkeypatch) -> None:
    monkeypatch.setattr(
        "zstarview.location_resolver.place_search.nominatim._build_url",
        lambda query, *, limit, countrycode: "https://example.invalid/search",
    )
    monkeypatch.setattr(
        "zstarview.location_resolver.place_search.nominatim._fetch",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(urllib.error.URLError("offline")),
    )

    with pytest.raises(urllib.error.URLError):
        search_place_candidates("Tokyo Station")
