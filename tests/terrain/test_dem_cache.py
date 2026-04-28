from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from zstarview.terrain.dem import (
    COPERNICUS_DEM_BUCKET,
    dem_tile_metadata_path,
    fetch_copernicus_dem,
    read_dem_tile_fetched_at_utc,
)


class _FakeS3Client:
    class exceptions:
        class NoSuchKey(Exception):
            pass

    def __init__(self, *, should_fail: bool) -> None:
        self.should_fail = should_fail
        self.calls: list[tuple[str, str]] = []

    def download_fileobj(self, bucket: str, key: str, handle) -> None:
        self.calls.append((bucket, key))
        if self.should_fail:
            raise RuntimeError("network down")
        handle.write(b"fake-dem")


class _FakeInvalidTileS3Client(_FakeS3Client):
    def download_fileobj(self, bucket: str, key: str, handle) -> None:
        self.calls.append((bucket, key))
        handle.write(b"not-a-geotiff")


def test_read_dem_tile_fetched_at_utc_migrates_legacy_cache(monkeypatch, tmp_path: Path) -> None:
    tile_path = tmp_path / "tile.tif"
    tile_path.write_bytes(b"legacy")
    now = datetime(2026, 3, 27, 2, 0, tzinfo=timezone.utc)
    monkeypatch.setattr("zstarview.terrain.dem._is_valid_dem_tile", lambda path: path.read_bytes() == b"legacy")

    fetched_at_utc = read_dem_tile_fetched_at_utc(tile_path, now_utc=now)

    assert fetched_at_utc == now
    assert dem_tile_metadata_path(tile_path).exists()


def test_read_dem_tile_fetched_at_utc_discards_invalid_cache(tmp_path: Path) -> None:
    tile_path = tmp_path / "tile.tif"
    tile_path.write_bytes(b"not-a-geotiff")
    dem_tile_metadata_path(tile_path).write_text(
        '{"fetched_at_utc": "2026-03-01T00:00:00+00:00"}',
        encoding="utf-8",
    )

    fetched_at_utc = read_dem_tile_fetched_at_utc(tile_path, now_utc=datetime(2026, 3, 27, 2, 0, tzinfo=timezone.utc))

    assert fetched_at_utc is None
    assert not tile_path.exists()
    assert not dem_tile_metadata_path(tile_path).exists()


def test_fetch_copernicus_dem_uses_stale_tile_when_refresh_fails(monkeypatch, tmp_path: Path) -> None:
    tile_path = tmp_path / "Copernicus_DSM_COG_30_N35_00_E139_00_DEM" / "Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"
    tile_path.parent.mkdir(parents=True, exist_ok=True)
    tile_path.write_bytes(b"stale")
    now = datetime(2026, 3, 27, 2, 0, tzinfo=timezone.utc)
    dem_tile_metadata_path(tile_path).write_text(
        "{\"fetched_at_utc\": \"%s\"}" % (now - timedelta(days=120)).isoformat(),
        encoding="utf-8",
    )
    fake_s3 = _FakeS3Client(should_fail=True)
    monkeypatch.setattr("zstarview.terrain.dem.anonymous_s3_client", lambda: fake_s3)
    monkeypatch.setattr("zstarview.terrain.dem.build_download_bbox", lambda **_kwargs: (0.0, 0.0, 1.0, 1.0))
    monkeypatch.setattr(
        "zstarview.terrain.dem.collect_copernicus_tile_keys",
        lambda _bbox: ["Copernicus_DSM_COG_30_N35_00_E139_00_DEM/Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"],
    )
    monkeypatch.setattr(
        "zstarview.terrain.dem._is_valid_dem_tile",
        lambda path: path.read_bytes() == b"stale",
    )

    got = fetch_copernicus_dem(
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        max_distance_km=1.0,
        margin_km=0.0,
        cache_dir=tmp_path,
        now_utc=now,
    )

    assert got.paths == (tile_path,)
    assert got.source == "cache-stale"
    assert fake_s3.calls == [
        (
            COPERNICUS_DEM_BUCKET,
            "Copernicus_DSM_COG_30_N35_00_E139_00_DEM/Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif",
        )
    ]


def test_fetch_copernicus_dem_discards_invalid_existing_tile(monkeypatch, tmp_path: Path) -> None:
    tile_relpath = "Copernicus_DSM_COG_30_N35_00_E139_00_DEM/Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"
    tile_path = tmp_path / tile_relpath
    tile_path.parent.mkdir(parents=True, exist_ok=True)
    tile_path.write_bytes(b"invalid")
    dem_tile_metadata_path(tile_path).write_text(
        '{"fetched_at_utc": "2026-03-01T00:00:00+00:00"}',
        encoding="utf-8",
    )
    fake_s3 = _FakeS3Client(should_fail=False)
    now = datetime(2026, 3, 27, 2, 0, tzinfo=timezone.utc)
    monkeypatch.setattr("zstarview.terrain.dem.anonymous_s3_client", lambda: fake_s3)
    monkeypatch.setattr("zstarview.terrain.dem.build_download_bbox", lambda **_kwargs: (0.0, 0.0, 1.0, 1.0))
    monkeypatch.setattr("zstarview.terrain.dem.collect_copernicus_tile_keys", lambda _bbox: [tile_relpath])
    monkeypatch.setattr(
        "zstarview.terrain.dem._is_valid_dem_tile",
        lambda path: path.suffix == ".tmp" or path.read_bytes() == b"fake-dem",
    )

    got = fetch_copernicus_dem(
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        max_distance_km=1.0,
        margin_km=0.0,
        cache_dir=tmp_path,
        now_utc=now,
    )

    assert got.paths == (tile_path,)
    assert got.source == "download"
    assert fake_s3.calls == [(COPERNICUS_DEM_BUCKET, tile_relpath)]
    assert tile_path.read_bytes() == b"fake-dem"


def test_fetch_copernicus_dem_skips_invalid_downloaded_tile(monkeypatch, tmp_path: Path) -> None:
    tile_relpath = "Copernicus_DSM_COG_30_N35_00_E139_00_DEM/Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"
    fake_s3 = _FakeInvalidTileS3Client(should_fail=False)
    monkeypatch.setattr("zstarview.terrain.dem.anonymous_s3_client", lambda: fake_s3)
    monkeypatch.setattr("zstarview.terrain.dem.build_download_bbox", lambda **_kwargs: (0.0, 0.0, 1.0, 1.0))
    monkeypatch.setattr("zstarview.terrain.dem.collect_copernicus_tile_keys", lambda _bbox: [tile_relpath])
    monkeypatch.setattr("zstarview.terrain.dem._is_valid_dem_tile", lambda path: path.read_bytes() == b"fake-dem")

    with pytest.raises(RuntimeError, match="No Copernicus DEM tiles were downloaded for the requested area."):
        fetch_copernicus_dem(
            observer_lat_deg=35.0,
            observer_lon_deg=139.0,
            max_distance_km=1.0,
            margin_km=0.0,
            cache_dir=tmp_path,
            now_utc=datetime(2026, 3, 27, 2, 0, tzinfo=timezone.utc),
        )

    assert fake_s3.calls == [(COPERNICUS_DEM_BUCKET, tile_relpath)]
    assert not (tmp_path / tile_relpath).exists()


def test_fetch_copernicus_dem_treats_legacy_tile_as_fresh_after_migration(monkeypatch, tmp_path: Path) -> None:
    tile_path = tmp_path / "Copernicus_DSM_COG_30_N35_00_E139_00_DEM" / "Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"
    tile_path.parent.mkdir(parents=True, exist_ok=True)
    tile_path.write_bytes(b"legacy")
    fake_s3 = _FakeS3Client(should_fail=False)
    now = datetime(2026, 3, 27, 2, 0, tzinfo=timezone.utc)
    monkeypatch.setattr("zstarview.terrain.dem.anonymous_s3_client", lambda: fake_s3)
    monkeypatch.setattr("zstarview.terrain.dem.build_download_bbox", lambda **_kwargs: (0.0, 0.0, 1.0, 1.0))
    monkeypatch.setattr(
        "zstarview.terrain.dem.collect_copernicus_tile_keys",
        lambda _bbox: ["Copernicus_DSM_COG_30_N35_00_E139_00_DEM/Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"],
    )
    monkeypatch.setattr(
        "zstarview.terrain.dem._is_valid_dem_tile",
        lambda path: path.read_bytes() == b"legacy",
    )

    got = fetch_copernicus_dem(
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        max_distance_km=1.0,
        margin_km=0.0,
        cache_dir=tmp_path,
        now_utc=now,
    )

    assert got.paths == (tile_path,)
    assert got.source == "cache"
    assert fake_s3.calls == []
    assert read_dem_tile_fetched_at_utc(tile_path, now_utc=now) == now
