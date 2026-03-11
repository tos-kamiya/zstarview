from __future__ import annotations

from zstarview.startup import _startup_resolve_city
from zstarview.viewpoints import Viewpoint


def test_startup_resolve_city_accepts_tower_name(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.startup.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.startup.save_last_city", lambda _value: None)
    location = _startup_resolve_city("Tokyo Skytree")
    assert location.kind == "tower"
    assert location.display_name == "t/Tokyo Skytree"
    assert location.persistence_key == "t/Tokyo Skytree"
    assert location.observer_height_m == 635.7
    assert abs(location.lat - 35.710055555) < 1e-6
    assert abs(location.lon - 139.810722222) < 1e-6
    assert location.tz == "Asia/Tokyo"


def test_startup_resolve_city_accepts_mountain_name(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.startup.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.startup.save_last_city", lambda _value: None)
    location = _startup_resolve_city("Mount Fuji")
    assert location.kind == "mountain"
    assert location.display_name == "m/Mount Fuji"
    assert location.persistence_key == "m/Mount Fuji"
    assert location.observer_height_m == 1.7
    assert abs(location.lat - 35.360555555) < 1e-6
    assert abs(location.lon - 138.7275) < 1e-6
    assert location.tz == "Asia/Tokyo"


def test_startup_resolve_city_accepts_mountain_wikidata_key(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.startup.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.startup.save_last_city", lambda _value: None)
    location = _startup_resolve_city("wikidata:Q39231")
    assert location.kind == "mountain"
    assert location.display_name == "m/Mount Fuji"


def test_startup_resolve_city_accepts_explicit_tower_prefix(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.startup.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.startup.save_last_city", lambda _value: None)
    location = _startup_resolve_city("t/Tokyo Skytree")
    assert location.kind == "tower"
    assert location.display_name == "t/Tokyo Skytree"


def test_startup_resolve_city_accepts_explicit_mountain_prefix(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.startup.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.startup.save_last_city", lambda _value: None)
    location = _startup_resolve_city("m/Mount Hermon")
    assert location.kind == "mountain"
    assert location.display_name == "m/Mount Hermon"


def test_startup_resolve_city_uses_observer_height_when_present(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.startup.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.startup.save_last_city", lambda _value: None)
    monkeypatch.setattr(
        "zstarview.startup.resolve_tower_viewpoint",
        lambda _name: Viewpoint(
            id="wikidata:Q1",
            qid="Q1",
            kind="tower",
            name="Example Deck",
            labels={},
            names=("Example Deck",),
            latitude_deg=35.0,
            longitude_deg=139.0,
            height_m=300.0,
            observer_height_m=240.0,
            meta={},
        ),
    )
    monkeypatch.setattr(
        "zstarview.startup._resolve_nearest_city",
        lambda _lat, _lon, _admin1_map: type("City", (), {"tz": "Asia/Tokyo", "cc": "JP"})(),
    )
    monkeypatch.setattr("zstarview.startup.load_admin1_names", lambda _path: {})

    location = _startup_resolve_city("Example Deck")

    assert location.kind == "tower"
    assert location.observer_height_m == 241.7


def test_startup_resolve_city_formats_city_display_name_with_country(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.startup.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.startup.save_last_city", lambda _value: None)
    location = _startup_resolve_city("Tokyo")
    assert location.kind == "city"
    assert location.display_name == "JP/Tokyo"
