from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr
from pyproj import CRS

from zstarview.clouddisc.geo_area import GeoArea
from zstarview.clouddisc.types import CloudSourceData
from zstarview.clouddisc.types import SourceKey
from zstarview.clouddisc.workers import cloud_source as cloud_source_worker
from zstarview.clouddisc.workers import cloud_source_worker as cloud_source_process_worker
from zstarview.clouddisc.workers.cloud_source import (
    CloudSourceFetchRequest,
    build_cloud_source_fetch_request,
    fetch_cloud_source,
)
from zstarview.clouddisc.workers.cloud_source_worker import load_cloud_source_worker_result


def test_cloud_source_fetch_request_normalizes_inputs() -> None:
    request = build_cloud_source_fetch_request(
        lat="35.5",
        lon="139.75",
        when_utc=datetime(2026, 6, 2, 1, 23, 45, tzinfo=timezone.utc),
        cloud_shells_km=(6374.0, 6376.0),
    )

    assert isinstance(request, CloudSourceFetchRequest)
    assert request.lat == pytest.approx(35.5)
    assert request.lon == pytest.approx(139.75)
    assert request.when_utc == datetime(2026, 6, 2, 1, 20, tzinfo=timezone.utc)
    assert request.cloud_shells_km == (6374.0, 6376.0)


def test_fetch_cloud_source_uses_himawari_provider(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    request = build_cloud_source_fetch_request(lat=35.0, lon=139.0)
    calls: list[tuple[str, object]] = []

    class _FakeHima:
        def fetch_bt_c13(self, **kwargs):  # noqa: ANN001
            calls.append(("hima", kwargs))
            da = SimpleNamespace(attrs={})
            return da, datetime(2026, 6, 2, 1, 20, tzinfo=timezone.utc), [tmp_path / "hima.nc"]

    class _FakeGoes:
        def fetch_bt_c13_with_failover(self, **kwargs):  # noqa: ANN001
            calls.append(("goes", kwargs))
            da = SimpleNamespace(attrs={})
            return (da, datetime(2026, 6, 2, 1, 20, tzinfo=timezone.utc), [tmp_path / "goes.nc"]), "G19"

    class _FakeContext:
        goes = _FakeGoes()
        hima = _FakeHima()

        def make_source_key(self, *, lat: float, lon: float, when_utc=None):  # noqa: ANN001
            del lat, lon, when_utc
            return SourceKey(
                satellite="HIMAWARI",
                provider="HIMAWARI",
                timeslot_utc=datetime(2026, 6, 2, 1, 20, tzinfo=timezone.utc),
            )

    result = fetch_cloud_source(_FakeContext(), request)

    assert result.satellite == "HIMAWARI"
    assert result.product == "ISatSS-B13"
    assert result.src_paths == [tmp_path / "hima.nc"]
    assert calls and calls[0][0] == "hima"
    assert calls[0][1]["observer_lat"] == pytest.approx(35.0)
    assert calls[0][1]["observer_lon"] == pytest.approx(139.0)


def test_fetch_cloud_source_uses_goes_provider(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    request = build_cloud_source_fetch_request(lat=35.0, lon=-90.0)
    calls: list[tuple[str, object]] = []

    class _FakeHima:
        def fetch_bt_c13(self, **kwargs):  # noqa: ANN001
            calls.append(("hima", kwargs))
            da = SimpleNamespace(attrs={})
            return da, datetime(2026, 6, 2, 1, 20, tzinfo=timezone.utc), [tmp_path / "hima.nc"]

    class _FakeGoes:
        def fetch_bt_c13_with_failover(self, **kwargs):  # noqa: ANN001
            calls.append(("goes", kwargs))
            da = SimpleNamespace(attrs={"source_expected_count": 2, "source_available_count": 2, "source_completeness_ratio": 1.0})
            return (da, datetime(2026, 6, 2, 1, 20, tzinfo=timezone.utc), [tmp_path / "goes.nc"]), "G19"

    class _FakeContext:
        goes = _FakeGoes()
        hima = _FakeHima()

        def make_source_key(self, *, lat: float, lon: float, when_utc=None):  # noqa: ANN001
            del lat, lon, when_utc
            return SourceKey(
                satellite="G19",
                provider="GOES",
                timeslot_utc=datetime(2026, 6, 2, 1, 20, tzinfo=timezone.utc),
            )

    monkeypatch.setattr(cloud_source_worker, "visible_satellites", lambda *_args, **_kwargs: ["G19", "G18"])

    result = fetch_cloud_source(_FakeContext(), request)

    assert result.satellite == "G19"
    assert result.product == "CMIPF-C13"
    assert result.source_expected_count == 2
    assert result.source_available_count == 2
    assert result.source_completeness_ratio == pytest.approx(1.0)
    assert result.src_paths == [tmp_path / "goes.nc"]
    assert calls and calls[0][0] == "goes"


def test_one_shot_worker_writes_manifest_and_artifact(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    area = GeoArea(
        area_id="test-area",
        description="test",
        proj_id="geos",
        projection=CRS.from_epsg(4326),
        width=2,
        height=2,
        area_extent=(0.0, 0.0, 1.0, 1.0),
    )
    source = CloudSourceData(
        source_key=SourceKey(
            satellite="HIMAWARI",
            provider="HIMAWARI",
            timeslot_utc=datetime(2026, 6, 2, 1, 20, tzinfo=timezone.utc),
        ),
        data_array=xr.DataArray(np.zeros((2, 2), dtype=np.float32), attrs={"area": area}),
        satellite="HIMAWARI",
        product="ISatSS-B13",
        time_utc=datetime(2026, 6, 2, 1, 20, tzinfo=timezone.utc),
        src_paths=[tmp_path / "source.nc"],
        source_expected_count=88,
        source_available_count=82,
        source_completeness_ratio=82.0 / 88.0,
    )

    class _FakeCloudDisc:
        def __init__(self, cfg):
            self.cfg = cfg

    monkeypatch.setattr(cloud_source_process_worker, "CloudDisc", _FakeCloudDisc)
    monkeypatch.setattr(cloud_source_process_worker, "fetch_cloud_source", lambda _context, _request: source)

    args = argparse.Namespace(
        work_dir=tmp_path,
        request_id=17,
        lat=35.0,
        lon=139.0,
        when_utc=None,
        cache_dir=tmp_path / "cache",
        sat_priority=["AUTO"],
        bt_warm_k=310.0,
        bt_cold_k=190.0,
        alt_min_deg=0.0,
        search_back_minutes=120,
        connect_timeout=5.0,
        read_timeout=30.0,
        cloud_shells_km=(6374.0, 6376.0, 6378.0),
    )

    assert cloud_source_process_worker._run_one_shot_worker(args=args) == 0

    loaded = load_cloud_source_worker_result(tmp_path / cloud_source_process_worker.WORKER_RESULT_FILENAME)
    assert loaded.satellite == "HIMAWARI"
    assert loaded.product == "ISatSS-B13"
    assert loaded.source_expected_count == 88
    assert loaded.source_available_count == 82
    assert loaded.source_completeness_ratio == pytest.approx(82.0 / 88.0)
    assert loaded.src_paths == [tmp_path / "source.nc"]
