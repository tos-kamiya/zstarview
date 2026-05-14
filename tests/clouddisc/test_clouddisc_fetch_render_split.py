from __future__ import annotations

import datetime as dt
from pathlib import Path
import numpy as np
import xarray as xr

from zstarview.clouddisc import CloudDisc, CloudDiscConfig
from zstarview.clouddisc.core import (
    DEFAULT_CLOUD_SHELLS_KM,
    DEFAULT_CLOUD_SHELL_BLEND_HIGH_CLOUD_AMOUNT,
    DEFAULT_CLOUD_SHELL_BLEND_LOW_CLOUD_AMOUNT,
    DEFAULT_CLOUD_SHELL_WEIGHTS,
    _blend_cloud_shell_weights,
)
from zstarview.clouddisc.types import CloudMeta, CloudSourceData, SourceKey
from zstarview.paths import CLOUD_SHELLS_KM


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

    def fake_fetch_source(*, lat: float, lon: float, when_utc=None, abort_event=None) -> CloudSourceData:
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


def test_default_cloud_shells_use_weighted_blend() -> None:
    assert CLOUD_SHELLS_KM == (6374.0, 6376.0, 6378.0)
    assert DEFAULT_CLOUD_SHELLS_KM == CLOUD_SHELLS_KM
    assert DEFAULT_CLOUD_SHELL_WEIGHTS == (0.20, 0.60, 0.20)


def test_cloud_shell_weights_transition_from_middle_only_to_default_mix() -> None:
    low = _blend_cloud_shell_weights(0.0, DEFAULT_CLOUD_SHELLS_KM)
    mid = _blend_cloud_shell_weights(
        0.5,
        DEFAULT_CLOUD_SHELLS_KM,
    )
    high = _blend_cloud_shell_weights(1.0, DEFAULT_CLOUD_SHELLS_KM)

    assert low == (0.0, 1.0, 0.0)
    assert high == DEFAULT_CLOUD_SHELL_WEIGHTS
    assert low[1] > mid[1] > high[1]
    assert low[0] < mid[0] < high[0]
    assert low[2] < mid[2] < high[2]
    assert DEFAULT_CLOUD_SHELL_BLEND_LOW_CLOUD_AMOUNT == 0.25
    assert DEFAULT_CLOUD_SHELL_BLEND_HIGH_CLOUD_AMOUNT == 0.65


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


def test_fetch_source_uses_max_shell_for_himawari_selection(monkeypatch, tmp_path: Path) -> None:
    clouddisc = CloudDisc(CloudDiscConfig(cache_dir=tmp_path))
    calls: list[float] = []

    monkeypatch.setattr(clouddisc, "_select_satellite", lambda _lat, _lon: "HIMAWARI")

    def fake_fetch_bt_c13(*, when_utc, observer_lat, observer_lon, cloud_shell_km, abort_event=None):
        calls.append(float(cloud_shell_km))
        return object(), when_utc, []

    monkeypatch.setattr(clouddisc.hima, "fetch_bt_c13", fake_fetch_bt_c13)

    source = clouddisc.fetch_source(
        lat=35.0,
        lon=139.0,
        when_utc=dt.datetime(2026, 3, 4, 12, 37, tzinfo=dt.timezone.utc),
        cloud_shells_km=(6374.0, 6379.0),
    )

    assert source.satellite == "HIMAWARI"
    assert calls == [6379.0]


def test_render_from_source_with_coverage_blends_multiple_shells(monkeypatch, tmp_path: Path) -> None:
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

    def fake_build_bt_sampler(_data_array):
        def _sampler(lon_grid, lat_grid):
            return lon_grid.astype(np.float32)

        return _sampler

    shell_masks = {
        6374.0: np.array([[True, True], [True, True]], dtype=bool),
        6376.0: np.array([[True, True], [True, True]], dtype=bool),
        6378.0: np.array([[True, True], [True, True]], dtype=bool),
    }
    shell_values = {6374.0: 6310.0, 6376.0: 6320.0, 6378.0: 6400.0}

    def fake_project(*, cloud_shell_km, **_kwargs):
        mask = shell_masks[float(cloud_shell_km)]
        lon_grid = np.full(mask.shape, shell_values[float(cloud_shell_km)], dtype=np.float32)
        lat_grid = np.zeros(mask.shape, dtype=np.float32)
        return lon_grid, lat_grid, mask

    def fake_convert(bt, mask_inside, bt_warm, bt_cold):
        out = np.zeros(mask_inside.shape + (4,), dtype=np.uint8)
        out[..., :3][mask_inside] = int(round(float(bt[mask_inside][0]) - 6300.0))
        out[..., 3][mask_inside] = 200
        return out

    monkeypatch.setattr("zstarview.clouddisc.core.build_bt_sampler", fake_build_bt_sampler)
    monkeypatch.setattr("zstarview.clouddisc.core.az_project_lonlat_grid", fake_project)
    monkeypatch.setattr(
        "zstarview.clouddisc.core.estimate_bt_warm_from_equator_band",
        lambda *_args, **_kwargs: (280.0, np.array([280.0], dtype=np.float32)),
    )
    monkeypatch.setattr(
        "zstarview.clouddisc.core.estimate_bt_cold_hybrid",
        lambda *_args, **_kwargs: 240.0,
    )
    monkeypatch.setattr("zstarview.clouddisc.core._estimate_scene_cloud_amount", lambda *_args, **_kwargs: 1.0)
    monkeypatch.setattr("zstarview.clouddisc.core.convert_bt_to_rgba_image", fake_convert)

    img, _meta, missing_mask, coverage_ratio = clouddisc.render_from_source_with_coverage(
        source=source,
        lat=35.0,
        lon=139.0,
        alt=45.0,
        az=180.0,
        radius_px=1,
        cloud_shells_km=(6374.0, 6376.0, 6378.0),
    )

    assert float(coverage_ratio) == 1.0
    assert int(np.count_nonzero(missing_mask)) == 0
    assert int(img[0, 0, 0]) == 34
    assert int(img[0, 1, 0]) == 34
    assert int(img[0, 0, 3]) == 200
    assert int(img[0, 1, 3]) == 200


def test_render_from_source_reuses_previous_bt_warm_when_equator_band_missing(monkeypatch, tmp_path: Path) -> None:
    clouddisc = CloudDisc(CloudDiscConfig(cache_dir=tmp_path))
    clouddisc._last_bt_warm = 287.5
    source = CloudSourceData(
        source_key=SourceKey(
            satellite="HIMAWARI",
            timeslot_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc),
            provider="HIMAWARI",
        ),
        data_array=xr.DataArray(
            np.zeros((2, 2), dtype=np.float32),
            dims=("y", "x"),
            attrs={"equator_band_missing": True},
        ),
        satellite="HIMAWARI",
        product="ISatSS-B13",
        time_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc),
        src_paths=[],
    )

    def fake_build_bt_sampler(_data_array):
        def _sampler(lon_grid, lat_grid):
            return np.full(lon_grid.shape, 250.0, dtype=np.float32)

        return _sampler

    captured_bt_warm: list[float] = []

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
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("equator-band estimator should not be called")),
    )
    monkeypatch.setattr(
        "zstarview.clouddisc.core.estimate_bt_cold_hybrid",
        lambda *args, **_kwargs: captured_bt_warm.append(float(args[3])) or 240.0,
    )

    clouddisc.render_from_source(
        source=source,
        lat=35.0,
        lon=139.0,
        alt=45.0,
        az=180.0,
        radius_px=1,
    )

    assert captured_bt_warm == [287.5]
    assert clouddisc._last_bt_warm == 287.5
