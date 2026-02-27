from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from zstarview.clouddisc.config import CloudDiscConfig
from zstarview.clouddisc.providers.meteosat import MeteosatFetchResult, MeteosatProvider
from zstarview.clouddisc.types import DataNotFoundError, DownloadError, RenderError


def _dummy_da() -> xr.DataArray:
    return xr.DataArray(np.zeros((2, 2), dtype=np.float32))


def test_fetch_bt_c13_returns_first_available_slot(tmp_path, monkeypatch) -> None:
    provider = MeteosatProvider(CloudDiscConfig(cache_dir=tmp_path, search_back_minutes=20))
    calls: list[dt.datetime] = []

    def fake_fetch_slot(search_time_utc: dt.datetime):
        calls.append(search_time_utc)
        if len(calls) < 2:
            return None
        return MeteosatFetchResult(
            da=_dummy_da(),
            used_time=search_time_utc,
            src_paths=[Path("/tmp/fake_meteosat.nc")],
            product="SEVIRI-IR",
        )

    monkeypatch.setattr(provider, "_fetch_slot", fake_fetch_slot)
    da, used_time, src_paths = provider.fetch_bt_c13(
        dt.datetime(2026, 2, 27, 12, 50, tzinfo=dt.timezone.utc)
    )

    assert isinstance(da, xr.DataArray)
    assert used_time.tzinfo is not None
    assert src_paths == [Path("/tmp/fake_meteosat.nc")]
    assert len(calls) == 2


def test_fetch_bt_c13_raises_data_not_found(tmp_path, monkeypatch) -> None:
    provider = MeteosatProvider(CloudDiscConfig(cache_dir=tmp_path, search_back_minutes=20))
    monkeypatch.setattr(provider, "_fetch_slot", lambda _t: None)

    with pytest.raises(DataNotFoundError):
        provider.fetch_bt_c13(dt.datetime(2026, 2, 27, 12, 50, tzinfo=dt.timezone.utc))


def test_fetch_bt_c13_propagates_download_error(tmp_path, monkeypatch) -> None:
    provider = MeteosatProvider(CloudDiscConfig(cache_dir=tmp_path, search_back_minutes=20))

    def raise_download(_t):
        raise DownloadError("mock download failed")

    monkeypatch.setattr(provider, "_fetch_slot", raise_download)

    with pytest.raises(DownloadError):
        provider.fetch_bt_c13(dt.datetime(2026, 2, 27, 12, 50, tzinfo=dt.timezone.utc))


def test_fetch_bt_c13_wraps_unexpected_exception_as_render_error(tmp_path, monkeypatch) -> None:
    provider = MeteosatProvider(CloudDiscConfig(cache_dir=tmp_path, search_back_minutes=20))

    def raise_unexpected(_t):
        raise ValueError("unexpected")

    monkeypatch.setattr(provider, "_fetch_slot", raise_unexpected)

    with pytest.raises(RenderError):
        provider.fetch_bt_c13(dt.datetime(2026, 2, 27, 12, 50, tzinfo=dt.timezone.utc))
