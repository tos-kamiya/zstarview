from __future__ import annotations

import datetime as dt
import pickle

from zstarview.clouddisc.types import CloudBandData, CloudSourceData, SourceKey

_TIME = dt.datetime(2026, 11, 20, 14, 30, tzinfo=dt.timezone.utc)


def _source() -> CloudSourceData:
    return CloudSourceData(
        source_key=SourceKey(
            satellite="HIMAWARI", provider="HIMAWARI", timeslot_utc=_TIME
        ),
        data_array=None,
        satellite="HIMAWARI",
        product="ISatSS-B13",
        time_utc=_TIME,
        src_paths=[],
        source_completeness_ratio=1.0,
    )


def test_cloud_source_data_pickle_round_trip() -> None:
    source = _source()
    restored = pickle.loads(pickle.dumps(source, protocol=pickle.HIGHEST_PROTOCOL))
    assert restored == source


def test_cloud_source_data_ignores_legacy_auxiliary_bands() -> None:
    source = _source()
    state = {name: getattr(source, name) for name in source.__slots__}
    state["auxiliary_bands"] = {
        "B16": CloudBandData(
            band=16,
            data_array=None,
            product="ISatSS-B16",
            time_utc=_TIME,
            src_paths=[],
        )
    }
    restored = CloudSourceData.__new__(CloudSourceData)
    # Pickles written before B16 removal carry (None, slot_state).
    restored.__setstate__((None, state))
    assert restored == source
    assert not hasattr(restored, "auxiliary_bands")
