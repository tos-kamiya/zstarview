from __future__ import annotations

import datetime as dt
from pathlib import Path
import numpy as np

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


def test_fetch_source_and_render_from_source_compose_current_render_flow(monkeypatch, tmp_path: Path) -> None:
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

    fetched = clouddisc.fetch_source(lat=35.0, lon=139.0)
    img, meta = clouddisc.render_from_source(
        source=fetched,
        lat=35.0,
        lon=139.0,
        alt=45.0,
        az=180.0,
        radius_px=256,
    )
    assert img == "fake_img"
    assert meta == expected_meta
    assert calls == {"fetch": 1, "render": 1}


def test_render_from_source_reuses_sampler_for_same_source(monkeypatch, tmp_path: Path) -> None:
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

    build_calls = {"count": 0}

    def fake_build_bt_sampler(_data_array):
        build_calls["count"] += 1

        def _sampler(lon_grid, lat_grid):
            return np.full(lon_grid.shape, 250.0, dtype=np.float32)

        return _sampler

    monkeypatch.setattr("zstarview.clouddisc.core.build_bt_sampler", fake_build_bt_sampler)
    monkeypatch.setattr(
        "zstarview.clouddisc.core.az_project_lonlat_grid",
        lambda **_kwargs: (
            np.array([[139.0, 139.1], [139.2, 139.3]], dtype=np.float32),
            np.array([[35.0, 35.1], [35.2, 35.3]], dtype=np.float32),
            np.array([[True, True], [True, True]], dtype=bool),
        ),
    )
    monkeypatch.setattr(
        "zstarview.clouddisc.core.estimate_bt_warm_from_equator_band",
        lambda *_args, **_kwargs: (280.0, np.array([280.0], dtype=np.float32)),
    )
    monkeypatch.setattr(
        "zstarview.clouddisc.core.estimate_bt_cold_hybrid",
        lambda *_args, **_kwargs: 240.0,
    )

    clouddisc.render_from_source(
        source=source,
        lat=35.0,
        lon=139.0,
        alt=45.0,
        az=180.0,
        radius_px=1,
    )
    clouddisc.render_from_source(
        source=source,
        lat=35.0,
        lon=139.0,
        alt=50.0,
        az=190.0,
        radius_px=1,
    )

    assert build_calls["count"] == 1
    assert source.sampler is not None
