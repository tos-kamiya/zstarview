from __future__ import annotations

from typing import Any

import pytest

from zstarview.location_resolver import LocationResolveError, resolve_launch_location
from zstarview.location_resolver import resolve as resolver

_IP_DATA = {
    "lat": 35.468,
    "lon": 133.0486,
    "city": "Matsue",
    "countryCode": "JP",
    "timezone": "Asia/Tokyo",
}


@pytest.fixture(autouse=True)
def _stub_ground_elevation(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resolver, "_resolve_ground_elevation_m", lambda **_kwargs: 0.0)


def _online(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resolver, "fetch_location_by_ip", lambda: dict(_IP_DATA))


def _resolve(city: str | None, last: str | dict[str, Any] | None) -> Any:
    return resolve_launch_location(
        city,
        load_last_city_func=lambda: last,
        save_last_city_func=lambda _value: None,
    )


def test_omitted_location_uses_auto_first(monkeypatch: pytest.MonkeyPatch) -> None:
    _online(monkeypatch)

    location = _resolve(None, "JP/Osaka")

    assert location.kind == "auto"
    assert location.display_name == "JP/Matsue"


def test_omitted_location_falls_back_to_last_when_auto_fails() -> None:
    location = _resolve("", "JP/Osaka")

    assert "Osaka" in location.display_name


def test_omitted_location_falls_back_to_tokyo_without_saved_location() -> None:
    location = _resolve(None, None)

    assert "Tokyo" in location.display_name


def test_last_keyword_uses_saved_location_even_when_auto_works(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _online(monkeypatch)

    location = _resolve(" LAST ", "JP/Osaka")

    assert "Osaka" in location.display_name


def test_last_keyword_restores_structured_auto_location() -> None:
    saved = {
        "resolver": "auto",
        "lat": 35.468,
        "lon": 133.0486,
        "display_name": "JP/Matsue",
        "tz": "Asia/Tokyo",
        "cc": "JP",
    }

    location = _resolve("last", saved)

    assert location.kind == "auto"
    assert location.display_name == "JP/Matsue"


def test_last_keyword_falls_back_to_tokyo_without_saved_location() -> None:
    location = _resolve("last", None)

    assert "Tokyo" in location.display_name


def test_explicit_auto_failure_is_still_an_error() -> None:
    with pytest.raises(LocationResolveError):
        _resolve("auto", "JP/Osaka")
