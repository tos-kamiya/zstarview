from __future__ import annotations

from pathlib import Path

from zstarview import cache_maintenance


def test_clear_long_lived_cache_removes_known_roots(monkeypatch, tmp_path: Path) -> None:
    cache_root = tmp_path / "cache"
    dem_root = cache_root / "copernicus-dem"
    overture_root = cache_root / "overture_buildings"
    skyscraper_root = cache_root / "overture_skyscrapers"
    for path in (dem_root, overture_root, skyscraper_root):
        path.mkdir(parents=True)
        (path / "marker.txt").write_text("x", encoding="utf-8")

    monkeypatch.setattr(cache_maintenance, "CACHE_PATH", str(cache_root))
    monkeypatch.setattr(cache_maintenance, "OVERTURE_DERIVED_ROOT_DIR", str(overture_root))
    monkeypatch.setattr(cache_maintenance, "OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR", str(skyscraper_root))

    removed = cache_maintenance.clear_long_lived_cache()

    assert removed == (dem_root, overture_root, skyscraper_root)
    assert not dem_root.exists()
    assert not overture_root.exists()
    assert not skyscraper_root.exists()
