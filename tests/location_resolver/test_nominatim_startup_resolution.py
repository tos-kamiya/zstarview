from __future__ import annotations

import pytest

from zstarview.location_resolver import LocationResolveError, resolve_launch_location


@pytest.fixture(autouse=True)
def _stub_ground_elevation(monkeypatch) -> None:
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_ground_elevation_m",
        lambda **_kwargs: 0.0,
    )


def test_startup_resolve_city_accepts_place_query(monkeypatch) -> None:
    saved_payloads: list[object] = []

    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda value: saved_payloads.append(value))
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_admin1_names", lambda _path: {})
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_nearest_city",
        lambda _lat, _lon, _admin1_map: type("City", (), {"tz": "Asia/Tokyo", "cc": "JP"})(),
    )
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve.search_nominatim",
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

    location = resolve_launch_location(
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

    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: saved_place)
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve.save_last_city",
        lambda _value: pytest.fail("save_last_city should not be called when restoring a saved place"),
    )
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_admin1_names", lambda _path: {})
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_nearest_city",
        lambda _lat, _lon, _admin1_map: type("City", (), {"tz": "Asia/Tokyo", "cc": "JP"})(),
    )
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve.search_nominatim",
        lambda *_args, **_kwargs: pytest.fail("search_nominatim should not be called"),
    )

    location = resolve_launch_location(None)

    assert location.kind == "place"
    assert location.display_name == "Matsue Station, Asahimachi, Matsue, Shimane, Japan"


def test_startup_resolve_city_raises_for_empty_place_results(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_admin1_names", lambda _path: {})
    monkeypatch.setattr("zstarview.location_resolver.resolve.search_nominatim", lambda *_args, **_kwargs: [])

    with pytest.raises(LocationResolveError):
        resolve_launch_location(None, place_query="No Such Station")
