from __future__ import annotations

import datetime as dt

import numpy as np
import pytest
import xarray as xr

from zstarview.clouddisc.config import CloudDiscConfig
from zstarview.clouddisc.providers import goes as goes_module
from zstarview.clouddisc.providers.goes import GoesProvider
from zstarview.clouddisc.types import DataNotFoundError


def _dummy_result() -> tuple[xr.DataArray, dt.datetime, list]:
    da = xr.DataArray(np.zeros((2, 2), dtype=np.float32))
    used = dt.datetime(2026, 2, 27, 12, 50, tzinfo=dt.timezone.utc)
    return da, used, []


def test_goes_failover_respects_allowed_satellites(tmp_path, monkeypatch) -> None:
    provider = GoesProvider(CloudDiscConfig(cache_dir=tmp_path))
    calls: list[str] = []

    def fake_once(
        sat: str,
        when_utc: dt.datetime,
        search_back_minutes: int,
        abort_event=None,
    ):
        del when_utc, search_back_minutes, abort_event
        calls.append(sat)
        if sat == "G18":
            return _dummy_result()
        return None

    monkeypatch.setattr(provider, "_fetch_bt_c13_once", fake_once)

    with pytest.raises(DataNotFoundError):
        provider.fetch_bt_c13_with_failover(
            sat="G19",
            when_utc=dt.datetime(2026, 2, 27, 13, 0, tzinfo=dt.timezone.utc),
            allowed_sats=("G19",),
        )

    assert set(calls) == {"G19"}


def test_goes_does_not_cache_current_hour_listing(tmp_path, monkeypatch) -> None:
    provider = GoesProvider(CloudDiscConfig(cache_dir=tmp_path))
    current_hour = dt.datetime(2026, 3, 31, 12, 0, tzinfo=dt.timezone.utc)
    listed_keys = [
        ["ABI-L2-CMIPF/2026/090/12/first.nc"],
        ["ABI-L2-CMIPF/2026/090/12/first.nc", "ABI-L2-CMIPF/2026/090/12/second.nc"],
    ]
    list_calls: list[str] = []

    class FrozenDateTime(dt.datetime):
        @classmethod
        def now(cls, tz=None):
            now = current_hour + dt.timedelta(minutes=24)
            if tz is None:
                return now.replace(tzinfo=None)
            return now.astimezone(tz)

    monkeypatch.setattr(goes_module.dt, "datetime", FrozenDateTime)

    def fake_list_s3_keys(**kwargs):
        del kwargs
        list_calls.append("list")
        return listed_keys.pop(0)

    monkeypatch.setattr(goes_module, "list_s3_keys", fake_list_s3_keys)

    first = provider._list_hour("noaa-goes19", current_hour)
    second = provider._list_hour("noaa-goes19", current_hour)

    assert len(list_calls) == 2
    assert first == ["ABI-L2-CMIPF/2026/090/12/first.nc"]
    assert second == [
        "ABI-L2-CMIPF/2026/090/12/first.nc",
        "ABI-L2-CMIPF/2026/090/12/second.nc",
    ]
