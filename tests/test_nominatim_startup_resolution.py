from __future__ import annotations

import pytest

from zstarview.startup import StartupAbortError, _startup_resolve_city


def test_startup_resolve_city_accepts_place_query(monkeypatch) -> None:
    saved_payloads: list[object] = []

    monkeypatch.setattr("zstarview.startup.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.startup.save_last_city", lambda value: saved_payloads.append(value))
    monkeypatch.setattr("zstarview.startup.load_admin1_names", lambda _path: {})
    monkeypatch.setattr(
        "zstarview.startup._resolve_nearest_city",
        lambda _lat, _lon, _admin1_map: type("City", (), {"tz": "Asia/Tokyo", "cc": "JP"})(),
    )
    monkeypatch.setattr(
        "zstarview.startup.search_nominatim",
        lambda query, **_kwargs: [
            {
                "name": "Matsue Station, Asahimachi, Matsue, Shimane, Japan",
                "lat": 35.4641778,
                "lon": 133.0628539,
                "category": "railway",
                "type": "station",
                "importance": 1.0,
            },
            {
                "name": "Another Matsue Station Candidate",
                "lat": 35.0,
                "lon": 133.0,
                "category": "place",
                "type": "village",
                "importance": 0.5,
            },
        ],
    )

    location = _startup_resolve_city(
        None,
        place_query="Matsue Station",
        place_countrycode="jp",
        place_lang="en",
    )

    assert location.kind == "place"
    assert location.display_name == "Matsue Station, Asahimachi, Matsue, Shimane, Japan"
    assert location.persistence_value is not None
    assert saved_payloads[0] == location.persistence_value
    assert location.tz == "Asia/Tokyo"


def test_startup_resolve_city_restores_saved_nominatim_place(monkeypatch) -> None:
    saved_place = {
        "resolver": "nominatim",
        "query": "Matsue Station",
        "result": {
            "name": "Matsue Station, Asahimachi, Matsue, Shimane, Japan",
            "lat": 35.4641778,
            "lon": 133.0628539,
            "category": "railway",
            "type": "station",
            "importance": 1.0,
        },
    }

    monkeypatch.setattr("zstarview.startup.load_last_city", lambda: saved_place)
    monkeypatch.setattr(
        "zstarview.startup.save_last_city",
        lambda _value: pytest.fail("save_last_city should not be called when restoring a saved place"),
    )
    monkeypatch.setattr("zstarview.startup.load_admin1_names", lambda _path: {})
    monkeypatch.setattr(
        "zstarview.startup._resolve_nearest_city",
        lambda _lat, _lon, _admin1_map: type("City", (), {"tz": "Asia/Tokyo", "cc": "JP"})(),
    )
    monkeypatch.setattr(
        "zstarview.startup.search_nominatim",
        lambda *_args, **_kwargs: pytest.fail("search_nominatim should not be called"),
    )

    location = _startup_resolve_city(None)

    assert location.kind == "place"
    assert location.display_name == "Matsue Station, Asahimachi, Matsue, Shimane, Japan"


def test_startup_resolve_city_raises_for_empty_place_results(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.startup.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.startup.load_admin1_names", lambda _path: {})
    monkeypatch.setattr("zstarview.startup.search_nominatim", lambda *_args, **_kwargs: [])

    with pytest.raises(StartupAbortError):
        _startup_resolve_city(None, place_query="No Such Station")
