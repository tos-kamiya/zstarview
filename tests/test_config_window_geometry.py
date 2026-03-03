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

