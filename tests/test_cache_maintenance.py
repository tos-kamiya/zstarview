from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

from zstarview import cache_maintenance


def test_clear_long_lived_cache_removes_known_roots(monkeypatch, tmp_path: Path) -> None:
    cache_root = tmp_path / "cache"
    metadata_path = cache_root / "clear_long_lived_cache_meta.json"
    dem_root = cache_root / "copernicus-dem"
    overture_root = cache_root / "overture_buildings"
    skyscraper_root = cache_root / "overture_skyscrapers"
    for path in (dem_root, overture_root, skyscraper_root):
        path.mkdir(parents=True)
        (path / "marker.txt").write_text("x", encoding="utf-8")

    monkeypatch.setattr(cache_maintenance, "CACHE_PATH", str(cache_root))
    monkeypatch.setattr(cache_maintenance, "OVERTURE_DERIVED_ROOT_DIR", str(overture_root))
    monkeypatch.setattr(cache_maintenance, "OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR", str(skyscraper_root))

    removed = cache_maintenance.clear_long_lived_cache(metadata_path=metadata_path)

    assert removed == (dem_root, overture_root, skyscraper_root)
    assert not dem_root.exists()
    assert not overture_root.exists()
    assert not skyscraper_root.exists()
    metadata_payload = metadata_path.read_text(encoding="utf-8")
    assert "last_cleared_at_utc" in metadata_payload


def test_clear_long_lived_cache_enforces_cooldown(monkeypatch, tmp_path: Path) -> None:
    cache_root = tmp_path / "cache"
    metadata_path = cache_root / "clear_long_lived_cache_meta.json"
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    now = datetime(2026, 3, 27, 3, 0, tzinfo=timezone.utc)
    metadata_path.write_text(
        '{\n  "last_cleared_at_utc": "%s"\n}' % (now - timedelta(days=1)).isoformat(),
        encoding="utf-8",
    )
    monkeypatch.setattr(cache_maintenance, "CACHE_PATH", str(cache_root))
    monkeypatch.setattr(cache_maintenance, "OVERTURE_DERIVED_ROOT_DIR", str(cache_root / "overture_buildings"))
    monkeypatch.setattr(cache_maintenance, "OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR", str(cache_root / "overture_skyscrapers"))

    try:
        cache_maintenance.clear_long_lived_cache(now_utc=now, metadata_path=metadata_path)
    except cache_maintenance.LongLivedCacheClearCooldownError as exc:
        assert "Long-lived cache was already cleared on" in str(exc)
        assert "Retry is allowed after" in str(exc)
        assert "remove these directories manually" in str(exc)
    else:
        raise AssertionError("Expected cooldown error")
