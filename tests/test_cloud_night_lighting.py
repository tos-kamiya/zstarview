from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from zstarview.cloud_night_lighting import (
    cloud_illumination_key,
    cloud_light_input_key,
    cloud_moon_input,
    moon_cloud_strength,
)
from zstarview.clouddisc.workers import cloud_voxel_worker
from zstarview.geosatellite.types import GeoSatelliteVoxelSource
from zstarview.render import cloud_ground_light
from zstarview.render.cloud_voxels import _night_lighting_inputs, _render


def test_lunar_phase_horizon_and_offscreen_input() -> None:
    assert moon_cloud_strength(60.0, 180.0) > moon_cloud_strength(60.0, 90.0) > 0
    assert moon_cloud_strength(60.0, 0.0) == 0
    assert moon_cloud_strength(-1.0, 180.0) == 0
    assert moon_cloud_strength(0.0, 180.0) == 0
    assert moon_cloud_strength(0.001, 180.0) < 1e-8
    assert moon_cloud_strength(float("nan"), 180.0) == 0
    assert cloud_moon_input([SimpleNamespace(
        name="moon", alt=60.0, az=120.0, phase_angle=180.0, is_visible=False
    )]) == (60.0, 120.0, 180.0)
    assert cloud_moon_input([]) is None


def test_cloud_light_keys_match_worker_inputs_and_quantize_small_changes() -> None:
    moon = SimpleNamespace(name="moon", alt=45.0, az=120.0, phase_angle=180.0)
    sun = SimpleNamespace(name="sun", alt=-20.0, az=220.0)
    assert cloud_illumination_key([sun, moon]) == cloud_light_input_key(
        (-20.0, 220.0), (45.0, 120.0, 180.0)
    )
    first = cloud_illumination_key([sun, moon])
    moon.alt = 45.01
    assert cloud_illumination_key([sun, moon]) == first
    moon.alt = 45.2
    assert cloud_illumination_key([sun, moon]) != first


def test_stale_lighting_completion_keeps_a_refresh_pending() -> None:
    from typing import Any

    from zstarview.gui.window_update_cloud import SkyWindowCloudUpdatesMixin

    moon = SimpleNamespace(name="moon", alt=45.0, az=120.0, phase_angle=180.0)
    owner: Any = SimpleNamespace(
        _is_shutting_down=False,
        state=SimpleNamespace(
            celestial_data=SimpleNamespace(planets=[moon]), cloud_projection_next_refresh_utc=None
        ),
    )
    stale = SkyWindowCloudUpdatesMixin._cloud_lighting_result_is_stale(
        owner, {"lighting_key": cloud_light_input_key(None, None)}
    )
    assert stale
    assert owner.state.cloud_projection_next_refresh_utc is not None
    owner.state.cloud_projection_next_refresh_utc = None
    assert not SkyWindowCloudUpdatesMixin._cloud_lighting_result_is_stale(
        owner, {"lighting_key": cloud_illumination_key([moon])}
    )
    assert owner.state.cloud_projection_next_refresh_utc is None


def test_first_celestial_lighting_schedules_cloud_refresh() -> None:
    from typing import Any

    from tests._window_render_support import (
        _DummyCompositor,
        _DummySignal,
        _DummyTimer,
        _WindowStub,
    )
    from zstarview.gui.window import SkyWindow
    from zstarview.gui.window_state import SkyWindowState
    from zstarview.render.geometry import get_screen_geometry
    from zstarview.types import ViewerData

    owner: Any = _WindowStub()
    owner.viewer_data = ViewerData(
        location=(35, 135), timezone_name="UTC", city_name="Test", view_center=(20, 30)
    )
    owner.state = SkyWindowState(render_view_center=(20, 30))
    owner.runtime_options = SimpleNamespace(cloud_voxel=True)
    owner._compositor = _DummyCompositor()
    owner._sky_data_update_timer = _DummyTimer(active=True)
    owner._cloud_update_timer = _DummyTimer(active=False)
    owner._clouddisc = None
    owner.cloud_disc_alpha = 0.2
    owner.sky_update_interval = 60
    owner.initial_data_loaded = _DummySignal()
    owner._is_shutting_down = False
    owner._disc_generation = 0
    owner.width = lambda: 640
    owner.height = lambda: 480
    owner.request_sky_data_update = lambda *_args, **_kwargs: None
    owner._safe_request_cloud_repaint = lambda: None
    owner.update = lambda: None
    SkyWindow._on_sky_data_calculated(owner, {
        "celestial": SimpleNamespace(planets=[SimpleNamespace(
            name="moon", alt=45.0, az=120.0, phase_angle=180.0
        )]),
        "sky_disc": object(), "view_center": (20, 30), "render_generation": 0,
        "geometry": get_screen_geometry(640, 480, owner.viewer_data.view_alt_deg),
    })
    assert owner.state.cloud_projection_next_refresh_utc is not None


def _scene(blocker=0.0, *, moon_strength=0.0, ground_strength=0.0, day_mix=0.0):
    density = np.zeros((2, 1, 2))
    density[0, 0, 1] = 0.4
    density[1, 0, 0] = blocker
    bases = np.repeat(np.eye(2)[None], 2, axis=0)
    # The top cell at x=0 lies directly above lower cell x=1, not x=0.
    centers = np.array([[-1.0, 0.0], [0.0, 0.0]])
    moon = np.array([1.0, 0.0, 0.05])
    moon /= np.linalg.norm(moon)
    return _render(
        density, np.array([0.5, -0.5, 1.5]), np.array([[0.0, 1.0, 0.0]]),
        moon, day_mix, np.array([0.36, 0.385, 0.44]), np.zeros((1, 3)), 1.0,
        np.zeros(2), bases, centers, bases, moon, moon_strength,
        np.full(density.shape, ground_strength),
    )


def test_vertical_ground_shadow_uses_physical_position_across_shifted_layers() -> None:
    clear_rgb, clear_t = _scene(ground_strength=1.0)
    blocked_rgb, blocked_t = _scene(4.0, ground_strength=1.0)
    ambient_rgb, _ = _scene(4.0)
    assert np.all(clear_rgb > blocked_rgb)
    assert np.all(blocked_rgb >= ambient_rgb)
    np.testing.assert_array_equal(clear_t, blocked_t)


def test_ground_and_moon_add_light_without_changing_view_extinction() -> None:
    ambient, transmission = _scene()
    lunar, lunar_t = _scene(moon_strength=0.25)
    ground, ground_t = _scene(ground_strength=1.0)
    combined, combined_t = _scene(moon_strength=0.25, ground_strength=1.0)
    assert np.all(lunar > ambient)
    assert np.all(ground > ambient)
    assert np.all(combined > lunar)
    assert np.all(combined > ground)
    for actual in (lunar_t, ground_t, combined_t):
        np.testing.assert_array_equal(actual, transmission)


def test_new_night_lighting_blends_continuously_and_does_not_affect_day() -> None:
    day, transmission = _scene(day_mix=1.0)
    lit_day, _ = _scene(day_mix=1.0, moon_strength=0.25, ground_strength=1.0)
    night, _ = _scene(moon_strength=0.25, ground_strength=1.0)
    twilight, twilight_t = _scene(day_mix=0.5, moon_strength=0.25, ground_strength=1.0)
    np.testing.assert_array_equal(day, lit_day)
    np.testing.assert_allclose(twilight, 0.5 * (day + night))
    np.testing.assert_array_equal(twilight_t, transmission)


def test_ground_sampler_skipped_in_day_and_inputs_validated() -> None:
    bases = np.eye(2)[None]
    def forbidden(*args):
        raise AssertionError("daytime must not read ground data")
    _, _, ground = _night_lighting_inputs(
        None, forbidden, 1.0, 0, 0, (1, 1, 1), np.zeros(2), np.zeros((1, 2)), bases
    )
    assert not ground.any()
    with pytest.raises(ValueError, match="finite"):
        _night_lighting_inputs(
            (float("nan"), 0, 180), None, 0, 0, 0, (1, 1, 1),
            np.zeros(2), np.zeros((1, 2)), bases,
        )


def test_original_raster_sampling_handles_nodata_edges_and_separate_blocks(tmp_path) -> None:
    path = tmp_path / "lights.tif"
    data = np.zeros((2, 1025), dtype=np.float32)
    data[0, 0] = 20.0
    data[0, 512] = 40.0
    data[0, 1024] = -9999.0
    with rasterio.open(
        path, "w", driver="GTiff", height=2, width=1025, count=1,
        dtype="float32", crs="EPSG:4326", transform=from_origin(0, 2, 0.001, 1),
        nodata=-9999.0,
    ) as dataset:
        dataset.write(data, 1)
    with rasterio.open(path) as dataset:
        values = cloud_ground_light._sample_raster(
            dataset, np.array([0.0005, 0.5125, 1.0245, -0.001]), np.full(4, 1.5)
        )
    np.testing.assert_array_equal(values, [20, 40, 0, 0])


def test_direct_ground_sampling_follows_layer_centers_and_tile_failure(monkeypatch) -> None:
    from contextlib import nullcontext
    positions = []
    def sample(dataset, longitude, latitude):
        positions.append(longitude.copy())
        return np.full(len(longitude), 20.0)
    monkeypatch.setattr(cloud_ground_light, "_ensure_night_light_tiles",
                        lambda **kwargs: {next(iter(kwargs["tile_names"])): "fake"})
    monkeypatch.setattr(cloud_ground_light.rasterio, "open", lambda _: nullcontext(None))
    monkeypatch.setattr(cloud_ground_light, "_sample_raster", sample)
    args = (35.0, 135.0, (1, 1, 2), np.array([0.5, 0.5]),
            np.array([[0.0, 0.0], [10.0, 0.0]]), np.repeat(np.eye(2)[None], 2, axis=0))
    ground = cloud_ground_light.sample_cloud_ground_light(*args)
    np.testing.assert_allclose(ground, 1 - np.exp(-1))
    assert positions[1][0] > positions[0][0]
    calls = []
    def fail(**kwargs):
        calls.append(kwargs)
        raise OSError("tile missing")
    monkeypatch.setattr(cloud_ground_light, "_ensure_night_light_tiles", fail)
    assert not cloud_ground_light.sample_cloud_ground_light(*args).any()
    assert len(calls) == 1  # No repeated failed acquisition in each layer.


def test_missing_ground_tile_does_not_discard_neighboring_tile(monkeypatch) -> None:
    from contextlib import nullcontext

    def fetch(**kwargs):
        tile = next(iter(kwargs["tile_names"]))
        if tile == "B1":
            raise OSError("west tile missing")
        return {tile: "fake"}

    monkeypatch.setattr(cloud_ground_light, "_ensure_night_light_tiles", fetch)
    monkeypatch.setattr(cloud_ground_light.rasterio, "open", lambda _: nullcontext(None))
    monkeypatch.setattr(cloud_ground_light, "_sample_raster",
                        lambda dataset, longitude, latitude: np.full(len(longitude), 20.0))
    ground = cloud_ground_light.sample_cloud_ground_light(
        35.0, 0.0, (2, 1, 1), np.array([1.0, 0.5]), np.zeros((1, 2)),
        (np.eye(2) * 10)[None],
    )
    assert ground[0, 0, 0] == 0
    assert ground[1, 0, 0] > 0


@pytest.mark.parametrize("geo", [False, True])
def test_worker_protocol_passes_lunar_and_vertical_ground_inputs(monkeypatch, geo, tmp_path) -> None:
    source = (
        GeoSatelliteVoxelSource(
            np.zeros((2, 2)), np.ones((2, 2), dtype=bool), tmp_path / "grid.npz",
            "infrared", datetime.now(timezone.utc), "lighting-source",
        ) if geo else SimpleNamespace(source_key="lighting-source")
    )
    calls = []
    def shade(source, lat, lon, alt, az, sun_alt, sun_az, base, **kwargs):
        calls.append(kwargs)
        return np.full_like(base, 0.1), {}, np.full(len(base), 0.5)
    def run(command, **kwargs):
        code = cloud_voxel_worker._worker_main(*(Path(value) for value in command[-3:]))
        return SimpleNamespace(returncode=code, stderr="")
    monkeypatch.setattr(cloud_voxel_worker, "shade_native_voxels", shade)
    monkeypatch.setattr(cloud_voxel_worker, "shade_geo_satellite_voxels", shade)
    monkeypatch.setattr(cloud_voxel_worker.subprocess, "run", run)
    rgba = cloud_voxel_worker.render_cloud_voxels_in_subprocess(
        source, lat=35, lon=135, view_center=(90, 0), edge_fov_deg=90,
        content_fov_deg=90, radius_px=2, sun_alt_deg=-20, sun_az_deg=0,
        sunlight_mix=0, night_color_rgb=(0.36, 0.385, 0.44), request_id=7,
        moon_altaz_phase=(45, 120, 180),
    )
    assert rgba.shape == (5, 5, 4)
    assert calls[0]["moon_altaz_phase"] == (45, 120, 180)
    assert calls[0]["ground_light_sampler"] is cloud_ground_light.sample_cloud_ground_light
    assert calls[0]["height_layer_transform"] is not geo
