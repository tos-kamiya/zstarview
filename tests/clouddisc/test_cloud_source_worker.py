from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import pytest

from zstarview.clouddisc.types import SourceKey
from zstarview.clouddisc.workers import cloud_source as cloud_source_worker
from zstarview.clouddisc.workers.cloud_source import (
    CloudSourceFetchRequest,
    build_cloud_source_fetch_request,
    fetch_cloud_source,
)


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
