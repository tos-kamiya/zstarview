from __future__ import annotations

from zstarview.startup import _startup_resolve_city


def test_startup_resolve_city_accepts_tower_name(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.startup.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.startup.save_last_city", lambda _value: None)
    location = _startup_resolve_city("Tokyo Skytree")
    assert location.kind == "tower"
    assert location.display_name == "t/Tokyo Skytree"
    assert location.persistence_key == "t/Tokyo Skytree"
    assert location.observer_height_m == 634.0
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


def test_startup_resolve_city_formats_city_display_name_with_country(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.startup.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.startup.save_last_city", lambda _value: None)
    location = _startup_resolve_city("Tokyo")
    assert location.kind == "city"
    assert location.display_name == "JP/Tokyo"
