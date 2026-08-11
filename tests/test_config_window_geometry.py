import json

from zstarview import config


def test_city_and_window_geometry_are_persisted_together(tmp_path, monkeypatch) -> None:
    cfg = tmp_path / "config.json"
    monkeypatch.setattr(config, "_config_file", cfg)

    config.save_last_city("JP/Tokyo")
    config.save_last_window_geometry(100, 200, 800, 600)

    data = json.loads(cfg.read_text(encoding="utf-8"))
    assert data["city"] == "JP/Tokyo"
    assert data["window_geometry"] == {"x": 100, "y": 200, "width": 800, "height": 600}


def test_save_last_city_keeps_existing_window_geometry(tmp_path, monkeypatch) -> None:
    cfg = tmp_path / "config.json"
    monkeypatch.setattr(config, "_config_file", cfg)
    cfg.write_text(
        json.dumps(
            {
                "city": "Old",
                "window_geometry": {"x": 1, "y": 2, "width": 3, "height": 4},
            },
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    config.save_last_city("JP/Osaka")
    assert config.load_last_city() == "JP/Osaka"
    assert config.load_last_window_geometry() == (1, 2, 3, 4)


def test_load_last_window_geometry_returns_none_for_invalid_data(tmp_path, monkeypatch) -> None:
    cfg = tmp_path / "config.json"
    monkeypatch.setattr(config, "_config_file", cfg)
    cfg.write_text(json.dumps({"window_geometry": {"x": 1, "y": 2, "width": "bad"}}), encoding="utf-8")

    assert config.load_last_window_geometry() is None


def test_save_last_city_accepts_structured_location_payload(tmp_path, monkeypatch) -> None:
    cfg = tmp_path / "config.json"
    monkeypatch.setattr(config, "_config_file", cfg)

    payload = {
        "resolver": "nominatim",
        "query": "Matsue Station",
        "result": {"name": "Matsue Station", "lat": 35.4, "lon": 133.0},
    }
    config.save_last_city(payload)

    assert config.load_last_city() == payload


def test_save_last_city_ignores_write_permission_errors(tmp_path, monkeypatch) -> None:
    cfg = tmp_path / "config.json"
    monkeypatch.setattr(config, "_config_file", cfg)

    def raise_permission_error(*_args, **_kwargs):
        raise PermissionError("blocked")

    monkeypatch.setattr(config.Path, "write_text", raise_permission_error)

    config.save_last_city("JP/Tokyo")


def test_save_last_window_geometry_ignores_write_permission_errors(tmp_path, monkeypatch) -> None:
    cfg = tmp_path / "config.json"
    monkeypatch.setattr(config, "_config_file", cfg)

    def raise_permission_error(*_args, **_kwargs):
        raise PermissionError("blocked")

    monkeypatch.setattr(config.Path, "write_text", raise_permission_error)

    config.save_last_window_geometry(1, 2, 3, 4)


def test_open_meteo_terms_confirmation_is_versioned_and_preserves_config(
    tmp_path, monkeypatch
) -> None:
    cfg = tmp_path / "config.json"
    monkeypatch.setattr(config, "_config_file", cfg)
    config.save_last_city("JP/Tokyo")

    assert config.open_meteo_noncommercial_terms_accepted() is False
    config.accept_open_meteo_noncommercial_terms()

    assert config.open_meteo_noncommercial_terms_accepted() is True
    data = json.loads(cfg.read_text(encoding="utf-8"))
    assert data["city"] == "JP/Tokyo"
    assert data["open_meteo_noncommercial_terms_version"] == 1
