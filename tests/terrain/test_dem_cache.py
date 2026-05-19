from __future__ import annotations

import io
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError

import pytest

from zstarview.terrain.dem import (
    DownloadCancelledError,
    dem_tile_metadata_path,
    fetch_copernicus_dem,
    read_dem_tile_fetched_at_utc,
)


def _tile_bytes(payload: bytes = b"fake-dem") -> bytes:
    return payload


class _FakeResponse:
    def __init__(self, payload: bytes) -> None:
        self._buffer = io.BytesIO(payload)

    def __enter__(self) -> "_FakeResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        return None

    def read(self, size: int = -1) -> bytes:
        return self._buffer.read(size)


class _FakeUrlopen:
    def __init__(self, responses: list[bytes | Exception]) -> None:
        self._responses = list(responses)
        self.calls: list[str] = []

    def __call__(self, request, timeout=None):  # noqa: ANN001
        self.calls.append(request.full_url)
        if not self._responses:
            raise AssertionError("Unexpected urlopen call")
        response = self._responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return _FakeResponse(response)


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
    fake_urlopen = _FakeUrlopen([URLError("network down")])
    monkeypatch.setattr("zstarview.terrain.dem.urlopen", fake_urlopen)
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
    assert len(fake_urlopen.calls) == 1


def test_fetch_copernicus_dem_discards_invalid_existing_tile(monkeypatch, tmp_path: Path) -> None:
    tile_relpath = "Copernicus_DSM_COG_30_N35_00_E139_00_DEM/Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"
    tile_path = tmp_path / tile_relpath
    tile_path.parent.mkdir(parents=True, exist_ok=True)
    tile_path.write_bytes(b"invalid")
    dem_tile_metadata_path(tile_path).write_text(
        '{"fetched_at_utc": "2026-03-01T00:00:00+00:00"}',
        encoding="utf-8",
    )
    fake_urlopen = _FakeUrlopen([_tile_bytes()])
    now = datetime(2026, 3, 27, 2, 0, tzinfo=timezone.utc)
    monkeypatch.setattr("zstarview.terrain.dem.urlopen", fake_urlopen)
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
    assert len(fake_urlopen.calls) == 1
    assert fake_urlopen.calls == [
        "https://copernicus-dem-90m.s3.eu-central-1.amazonaws.com/Copernicus_DSM_COG_30_N35_00_E139_00_DEM/Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"
    ]
    assert tile_path.read_bytes() == b"fake-dem"


def test_fetch_copernicus_dem_skips_invalid_downloaded_tile(monkeypatch, tmp_path: Path) -> None:
    tile_relpath = "Copernicus_DSM_COG_30_N35_00_E139_00_DEM/Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"
    fake_urlopen = _FakeUrlopen([b"not-a-geotiff"])
    monkeypatch.setattr("zstarview.terrain.dem.urlopen", fake_urlopen)
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

    assert fake_urlopen.calls == [
        "https://copernicus-dem-90m.s3.eu-central-1.amazonaws.com/Copernicus_DSM_COG_30_N35_00_E139_00_DEM/Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"
    ]
    assert not (tmp_path / tile_relpath).exists()


def test_fetch_copernicus_dem_404_without_cache_reports_no_tiles(monkeypatch, tmp_path: Path) -> None:
    tile_relpath = "Copernicus_DSM_COG_30_N35_00_E139_00_DEM/Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"
    fake_urlopen = _FakeUrlopen([
        HTTPError(
            url="https://copernicus-dem-90m.s3.eu-central-1.amazonaws.com/",
            code=404,
            msg="Not Found",
            hdrs=None,
            fp=None,
        )
    ])
    monkeypatch.setattr("zstarview.terrain.dem.urlopen", fake_urlopen)
    monkeypatch.setattr("zstarview.terrain.dem.build_download_bbox", lambda **_kwargs: (0.0, 0.0, 1.0, 1.0))
    monkeypatch.setattr("zstarview.terrain.dem.collect_copernicus_tile_keys", lambda _bbox: [tile_relpath])

    with pytest.raises(RuntimeError, match="No Copernicus DEM tiles were downloaded for the requested area."):
        fetch_copernicus_dem(
            observer_lat_deg=35.0,
            observer_lon_deg=139.0,
            max_distance_km=1.0,
            margin_km=0.0,
            cache_dir=tmp_path,
            now_utc=datetime(2026, 3, 27, 2, 0, tzinfo=timezone.utc),
        )

    assert len(fake_urlopen.calls) == 1


def test_fetch_copernicus_dem_treats_legacy_tile_as_fresh_after_migration(monkeypatch, tmp_path: Path) -> None:
    tile_path = tmp_path / "Copernicus_DSM_COG_30_N35_00_E139_00_DEM" / "Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"
    tile_path.parent.mkdir(parents=True, exist_ok=True)
    tile_path.write_bytes(b"legacy")
    fake_urlopen = _FakeUrlopen([AssertionError("must not be called")])
    now = datetime(2026, 3, 27, 2, 0, tzinfo=timezone.utc)
    monkeypatch.setattr("zstarview.terrain.dem.urlopen", fake_urlopen)
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
    assert fake_urlopen.calls == []
    assert read_dem_tile_fetched_at_utc(tile_path, now_utc=now) == now


def test_fetch_copernicus_dem_can_be_cancelled(monkeypatch, tmp_path: Path) -> None:
    abort_event = threading.Event()

    class _StreamingResponse:
        def __init__(self) -> None:
            self._chunks = [b"part1", b"part2"]

        def __enter__(self) -> "_StreamingResponse":
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

        def read(self, size: int = -1) -> bytes:
            if self._chunks:
                chunk = self._chunks.pop(0)
                if chunk == b"part1":
                    abort_event.set()
                return chunk
            return b""

    class _AbortAfterFirstCall:
        def __init__(self) -> None:
            self.calls = 0

        def __call__(self, request, timeout=None):  # noqa: ANN001
            self.calls += 1
            return _StreamingResponse()

    fake_urlopen = _AbortAfterFirstCall()
    monkeypatch.setattr("zstarview.terrain.dem.urlopen", fake_urlopen)
    monkeypatch.setattr("zstarview.terrain.dem.build_download_bbox", lambda **_kwargs: (0.0, 0.0, 1.0, 1.0))
    monkeypatch.setattr(
        "zstarview.terrain.dem.collect_copernicus_tile_keys",
        lambda _bbox: ["Copernicus_DSM_COG_30_N35_00_E139_00_DEM/Copernicus_DSM_COG_30_N35_00_E139_00_DEM.tif"],
    )

    with pytest.raises(DownloadCancelledError):
        fetch_copernicus_dem(
            observer_lat_deg=35.0,
            observer_lon_deg=139.0,
            max_distance_km=1.0,
            margin_km=0.0,
            cache_dir=tmp_path,
            now_utc=datetime(2026, 3, 27, 2, 0, tzinfo=timezone.utc),
            abort_event=abort_event,
        )

    assert fake_urlopen.calls == 1
