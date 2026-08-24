from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from zstarview.data.urban_outline_common import BuildingFootprint
from zstarview.location_resolver import LocationResolveError, resolve_launch_location
from zstarview.location_resolver.viewpoints import Viewpoint


def _write_overture_building_cache(
    root,
    *,
    lat_deg: float,
    lon_deg: float,
    radius_km: float,
    feature_type: str,
    height_m: float | None,
):
    from zstarview.data.import_overture_buildings import derive_dataset_name

    dataset_name = derive_dataset_name(
        lat_deg,
        lon_deg,
        radius_km,
        feature_type,
        0.0,
    )
    derived_dir = root / dataset_name / "bldg"
    derived_dir.mkdir(parents=True)
    buildings = []
    if height_m is not None:
        buildings.append(
            {
                "id": f"{feature_type}-1",
                "height_m": height_m,
                "rings": [
                    [
                        [lon_deg - 0.0001, lat_deg - 0.0001],
                        [lon_deg + 0.0001, lat_deg - 0.0001],
                        [lon_deg + 0.0001, lat_deg + 0.0001],
                        [lon_deg - 0.0001, lat_deg + 0.0001],
                        [lon_deg - 0.0001, lat_deg - 0.0001],
                    ]
                ],
            }
        )
    (derived_dir / f"{dataset_name}.json").write_text(
        json.dumps(
            {
                "tile": {
                    "bbox": {
                        "min_lat": lat_deg - 0.01,
                        "min_lon": lon_deg - 0.01,
                        "max_lat": lat_deg + 0.01,
                        "max_lon": lon_deg + 0.01,
                    }
                },
                "buildings": buildings,
            }
        ),
        encoding="utf-8",
    )
    (derived_dir.parent / "cache_meta.json").write_text(
        json.dumps(
            {
                "fetched_at_utc": "2099-01-01T00:00:00+00:00",
                "overture_release": "test-release",
                "query_lat_deg": lat_deg,
                "query_lon_deg": lon_deg,
                "query_radius_km": radius_km,
            }
        ),
        encoding="utf-8",
    )
    return derived_dir


@pytest.fixture(autouse=True)
def _stub_ground_elevation(monkeypatch) -> None:
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_ground_elevation_m",
        lambda **_kwargs: 0.0,
    )


def test_startup_resolve_city_accepts_tower_name(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    location = resolve_launch_location("Tokyo Skytree")
    assert location.kind == "tower"
    assert location.display_name == "t/Tokyo Skytree"
    assert location.persistence_key == "t/Tokyo Skytree"
    assert location.observer_height_m == 635.7
    assert abs(location.lat - 35.710055555) < 1e-6
    assert abs(location.lon - 139.810722222) < 1e-6
    assert location.tz == "Asia/Tokyo"


def test_startup_resolve_city_accepts_added_famous_tower_name(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    location = resolve_launch_location("N Seoul Tower")
    assert location.kind == "tower"
    assert location.display_name == "t/N Seoul Tower"
    assert abs(location.lat - 37.551216) < 1e-6
    assert abs(location.lon - 126.988276) < 1e-6


def test_startup_resolve_city_accepts_mountain_name(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    location = resolve_launch_location("Mount Fuji")
    assert location.kind == "mountain"
    assert location.display_name == "m/Mount Fuji"
    assert location.persistence_key == "m/Mount Fuji"
    assert location.observer_height_m == 1.7
    assert abs(location.lat - 35.360555555) < 1e-6
    assert abs(location.lon - 138.7275) < 1e-6
    assert location.tz == "Asia/Tokyo"


def test_startup_resolve_city_accepts_added_famous_mountain_name(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    location = resolve_launch_location("Kilimanjaro")

    assert location.kind == "mountain"
    assert location.display_name == "m/Mount Kilimanjaro"
    assert abs(location.lat - (-3.0666666666667)) < 1e-6
    assert abs(location.lon - 37.359166666667) < 1e-6


def test_startup_resolve_city_accepts_mountain_wikidata_key(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    location = resolve_launch_location("wikidata:Q39231")
    assert location.kind == "mountain"
    assert location.display_name == "m/Mount Fuji"


def test_startup_resolve_city_accepts_explicit_tower_prefix(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    location = resolve_launch_location("t/Tokyo Skytree")
    assert location.kind == "tower"
    assert location.display_name == "t/Tokyo Skytree"


def test_startup_resolve_city_accepts_explicit_mountain_prefix(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    location = resolve_launch_location("m/Mount Hermon")
    assert location.kind == "mountain"
    assert location.display_name == "m/Mount Hermon"


def test_startup_resolve_city_uses_viewpoint_height_when_present(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve.resolve_tower_viewpoint",
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
            viewpoint_height_m=240.0,
            meta={},
        ),
    )
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_nearest_city",
        lambda _lat, _lon, _admin1_map: type("City", (), {"tz": "Asia/Tokyo", "cc": "JP"})(),
    )
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_admin1_names", lambda _path: {})

    location = resolve_launch_location("Example Deck")

    assert location.kind == "tower"
    assert location.observer_height_m == 241.7


def test_startup_resolve_city_formats_city_display_name_with_country(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    location = resolve_launch_location("Tokyo")
    assert location.kind == "city"
    assert location.display_name == "JP/Tokyo"


def test_startup_resolve_city_accepts_at_lat_lon(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_admin1_names", lambda _path: {})
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_nearest_city",
        lambda _lat, _lon, _admin1_map: type("City", (), {"tz": "Asia/Tokyo", "cc": "JP"})(),
    )

    location = resolve_launch_location("@35.4824704,133.0683567")

    assert location.kind == "coords"
    assert location.display_name == "Lat: 35.4824704, Lon: 133.0683567"
    assert location.persistence_key == "35.482470;133.068357"
    assert location.tz == "Asia/Tokyo"
    assert location.cc == "JP"


def test_startup_resolve_city_accepts_plain_lat_lon_with_space(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_admin1_names", lambda _path: {})
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_nearest_city",
        lambda _lat, _lon, _admin1_map: type("City", (), {"tz": "Asia/Tokyo", "cc": "JP"})(),
    )

    location = resolve_launch_location("35.4824704, 133.0683567")

    assert location.kind == "coords"
    assert location.display_name == "Lat: 35.4824704, Lon: 133.0683567"
    assert location.persistence_key == "35.482470;133.068357"
    assert location.tz == "Asia/Tokyo"
    assert location.cc == "JP"


def test_startup_resolve_city_accepts_at_lat_lon_with_space(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_admin1_names", lambda _path: {})
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_nearest_city",
        lambda _lat, _lon, _admin1_map: type("City", (), {"tz": "Asia/Tokyo", "cc": "JP"})(),
    )

    location = resolve_launch_location("@35.4824704, 133.0683567")

    assert location.kind == "coords"
    assert location.display_name == "Lat: 35.4824704, Lon: 133.0683567"
    assert location.persistence_key == "35.482470;133.068357"
    assert location.tz == "Asia/Tokyo"
    assert location.cc == "JP"


def test_startup_resolve_city_accepts_google_maps_url(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_admin1_names", lambda _path: {})
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_nearest_city",
        lambda _lat, _lon, _admin1_map: type("City", (), {"tz": "Asia/Tokyo", "cc": "JP"})(),
    )

    location = resolve_launch_location(
        "https://www.google.com/maps/place/%E7%A5%9E%E6%88%B8%E4%B8%89%E5%AE%AE%E9%A7%85/"
        "@34.6938393,135.1914038,15.58z/data=!4m6!3m5!"
        "1s0x60008efaac8a2383:0xf35028084a4989de!8m2!3d34.6938491!4d135.1960454"
        "!16s%2Fg%2F1tdx_rsh?entry=ttu"
    )

    assert location.kind == "coords"
    assert location.display_name == "Lat: 34.6938491, Lon: 135.1960454"
    assert location.persistence_key == "34.693849;135.196045"
    assert location.tz == "Asia/Tokyo"


def test_startup_resolve_city_uses_building_top_when_enabled(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_admin1_names", lambda _path: {})
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_nearest_city",
        lambda _lat, _lon, _admin1_map: type("City", (), {"tz": "Asia/Tokyo", "cc": "JP"})(),
    )
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_building_top_height_m",
        lambda **_kwargs: 42.0,
    )

    location = resolve_launch_location("@35.4824704,133.0683567", use_building_top=True)

    assert location.kind == "coords"
    assert location.observer_height_m == 43.7
    assert location.location_height_label == "Building"
    assert location.location_height_m == 42.0


def test_startup_resolve_city_keeps_default_height_when_no_building_found(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_admin1_names", lambda _path: {})
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_nearest_city",
        lambda _lat, _lon, _admin1_map: type("City", (), {"tz": "Asia/Tokyo", "cc": "JP"})(),
    )
    monkeypatch.setattr(
        "zstarview.location_resolver.resolve._resolve_building_top_height_m",
        lambda **_kwargs: None,
    )

    location = resolve_launch_location("@35.4824704,133.0683567", use_building_top=True)

    assert location.kind == "coords"
    assert location.observer_height_m == 1.7
    assert location.location_height_label is None
    assert location.location_height_m == 0.0


def test_find_building_top_height_m_accepts_nearby_building_within_5m() -> None:
    from zstarview.location_resolver.resolve import _find_building_top_height_m

    buildings = (
        BuildingFootprint(
            building_id="near-lower",
            height_m=20.0,
            rings_lonlat=(
                (
                    (139.000010, 35.000000),
                    (139.000060, 35.000000),
                    (139.000060, 35.000050),
                    (139.000010, 35.000050),
                    (139.000010, 35.000000),
                ),
            ),
        ),
        BuildingFootprint(
            building_id="near-higher",
            height_m=45.0,
            rings_lonlat=(
                (
                    (139.000020, 35.000000),
                    (139.000070, 35.000000),
                    (139.000070, 35.000050),
                    (139.000020, 35.000050),
                    (139.000020, 35.000000),
                ),
            ),
        ),
    )

    got = _find_building_top_height_m(
        buildings,
        lon_deg=139.0,
        lat_deg=35.000025,
    )

    assert got == 45.0


def test_building_top_reuses_normal_overture_cache(monkeypatch, tmp_path) -> None:
    from zstarview.location_resolver import resolve as resolve_module

    lat_deg = 35.65305
    lon_deg = 139.67368
    for feature_type, height_m in (("building", 42.0), ("building_part", None)):
        _write_overture_building_cache(
            tmp_path,
            lat_deg=lat_deg,
            lon_deg=lon_deg,
            radius_km=2.5,
            feature_type=feature_type,
            height_m=height_m,
        )
    monkeypatch.setattr(resolve_module, "OVERTURE_DERIVED_ROOT_DIR", tmp_path)
    monkeypatch.setattr(
        resolve_module,
        "select_prepared_building_source",
        lambda **_kwargs: SimpleNamespace(source="overture", derived_dirs=()),
    )
    monkeypatch.setattr(
        "zstarview.data.import_overture_buildings.resolve_overture_release_for_cache_root",
        lambda **_kwargs: "test-release",
    )
    monkeypatch.setattr(
        "zstarview.data.import_overture_buildings.import_overture_buildings",
        lambda **_kwargs: pytest.fail("fresh 2.5km cache should avoid download"),
    )

    assert (
        resolve_module._resolve_building_top_height_m(
            lat_deg=lat_deg,
            lon_deg=lon_deg,
        )
        == 42.0
    )


def test_building_top_reuses_small_overture_cache(monkeypatch, tmp_path) -> None:
    from zstarview.location_resolver import resolve as resolve_module

    lat_deg = 35.65305
    lon_deg = 139.67368
    for feature_type, height_m in (("building", 30.0), ("building_part", 48.0)):
        _write_overture_building_cache(
            tmp_path,
            lat_deg=lat_deg,
            lon_deg=lon_deg,
            radius_km=0.15,
            feature_type=feature_type,
            height_m=height_m,
        )
    monkeypatch.setattr(resolve_module, "OVERTURE_DERIVED_ROOT_DIR", tmp_path)
    monkeypatch.setattr(
        resolve_module,
        "select_prepared_building_source",
        lambda **_kwargs: SimpleNamespace(source="overture", derived_dirs=()),
    )
    monkeypatch.setattr(
        "zstarview.data.import_overture_buildings.resolve_overture_release_for_cache_root",
        lambda **_kwargs: "test-release",
    )
    monkeypatch.setattr(
        "zstarview.data.import_overture_buildings.import_overture_buildings",
        lambda **_kwargs: pytest.fail("fresh 0.15km cache should avoid download"),
    )

    assert (
        resolve_module._resolve_building_top_height_m(
            lat_deg=lat_deg,
            lon_deg=lon_deg,
        )
        == 48.0
    )


def test_building_top_fetches_only_missing_small_cache(monkeypatch, tmp_path) -> None:
    from zstarview.location_resolver import resolve as resolve_module

    lat_deg = 35.65305
    lon_deg = 139.67368
    _write_overture_building_cache(
        tmp_path,
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        radius_km=0.15,
        feature_type="building",
        height_m=35.0,
    )
    fetched_types: list[str] = []

    def _fetch(**kwargs):
        feature_type = str(kwargs["feature_type"])
        fetched_types.append(feature_type)
        return _write_overture_building_cache(
            tmp_path,
            lat_deg=lat_deg,
            lon_deg=lon_deg,
            radius_km=0.15,
            feature_type=feature_type,
            height_m=None,
        )

    monkeypatch.setattr(resolve_module, "OVERTURE_DERIVED_ROOT_DIR", tmp_path)
    monkeypatch.setattr(
        resolve_module,
        "select_prepared_building_source",
        lambda **_kwargs: SimpleNamespace(source="overture", derived_dirs=()),
    )
    monkeypatch.setattr(
        "zstarview.data.import_overture_buildings.resolve_overture_release_for_cache_root",
        lambda **_kwargs: "test-release",
    )
    monkeypatch.setattr(
        "zstarview.data.import_overture_buildings.import_overture_buildings",
        _fetch,
    )

    assert (
        resolve_module._resolve_building_top_height_m(
            lat_deg=lat_deg,
            lon_deg=lon_deg,
        )
        == 35.0
    )
    assert fetched_types == ["building_part"]


def test_startup_resolve_city_rejects_google_maps_url_without_coordinates(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)

    with pytest.raises(LocationResolveError):
        resolve_launch_location("https://www.google.com/maps/place/Tokyo+Tower")


def test_startup_resolve_city_rejects_google_com_maps_host(monkeypatch) -> None:
    monkeypatch.setattr("zstarview.location_resolver.resolve.load_last_city", lambda: None)
    monkeypatch.setattr("zstarview.location_resolver.resolve.save_last_city", lambda _value: None)

    with pytest.raises(LocationResolveError):
        resolve_launch_location("https://google.com/maps/@35.0,139.0,17z")
