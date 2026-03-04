from __future__ import annotations

import datetime as dt
from pathlib import Path

from zstarview.clouddisc import CloudDisc, CloudDiscConfig
from zstarview.clouddisc.types import CloudMeta, CloudSourceData, SourceKey


def test_make_source_key_rounds_time_slot(tmp_path: Path) -> None:
    clouddisc = CloudDisc(CloudDiscConfig(cache_dir=tmp_path))
    key = clouddisc.make_source_key(
        lat=35.0,
        lon=139.0,
        when_utc=dt.datetime(2026, 3, 4, 12, 37, 12, tzinfo=dt.timezone.utc),
    )
    assert key.timeslot_utc == dt.datetime(2026, 3, 4, 12, 30, 0, tzinfo=dt.timezone.utc)


def test_render_now_delegates_to_fetch_and_render(monkeypatch, tmp_path: Path) -> None:
    clouddisc = CloudDisc(CloudDiscConfig(cache_dir=tmp_path))
    source = CloudSourceData(
        source_key=SourceKey(
            satellite="HIMAWARI",
            timeslot_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc),
            provider="HIMAWARI",
        ),
        data_array=object(),
        satellite="HIMAWARI",
        product="ISatSS-B13",
        time_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc),
        src_paths=[],
    )
    expected_meta = CloudMeta(
        satellite="HIMAWARI",
        product="ISatSS-B13",
        time_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc),
        src_paths=[],
    )
    calls = {"fetch": 0, "render": 0}

    def fake_fetch_source(*, lat: float, lon: float, when_utc=None) -> CloudSourceData:
        assert lat == 35.0
        assert lon == 139.0
        assert when_utc is None
        calls["fetch"] += 1
        return source

    def fake_render_from_source(**kwargs):
        calls["render"] += 1
        assert kwargs["source"] is source
        return "fake_img", expected_meta

    monkeypatch.setattr(clouddisc, "fetch_source", fake_fetch_source)
    monkeypatch.setattr(clouddisc, "render_from_source", fake_render_from_source)

    img, meta = clouddisc.render_now(
        lat=35.0,
        lon=139.0,
        alt=45.0,
        az=180.0,
        radius_px=256,
    )
    assert img == "fake_img"
    assert meta == expected_meta
    assert calls == {"fetch": 1, "render": 1}
