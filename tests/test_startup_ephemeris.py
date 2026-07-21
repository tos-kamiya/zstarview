from __future__ import annotations

import io
import sys
from datetime import timedelta
from types import SimpleNamespace

import numpy as np
import pytest

from zstarview.astro import _starfield_load, load_ephemeris
from zstarview.gui.viewer import _verify_ephemeris_for_launch
from zstarview.launch_time import LaunchSetupError
from zstarview.render import geometry as render_geometry
from zstarview.types import ViewerData


def test_de442s_uses_naif_planets_url() -> None:
    from zstarview.paths import EPHEMERIS_URL

    assert _starfield_load.build_url("de442s.bsp") == EPHEMERIS_URL


def test_startup_verify_ephemeris_passes_when_loader_succeeds(monkeypatch) -> None:
    calls: list[str] = []

    def fake_load_ephemeris() -> object:
        calls.append("de442s.bsp")
        return object()

    monkeypatch.setattr("zstarview.gui.viewer.load_ephemeris", fake_load_ephemeris)

    _verify_ephemeris_for_launch()

    assert calls == ["de442s.bsp"]


def test_startup_verify_ephemeris_aborts_on_oserror(monkeypatch) -> None:
    def fake_load_ephemeris() -> object:
        raise OSError("network blocked")

    monkeypatch.setattr("zstarview.gui.viewer.load_ephemeris", fake_load_ephemeris)

    with pytest.raises(LaunchSetupError):
        _verify_ephemeris_for_launch()


def test_load_ephemeris_provides_dummy_standard_streams_when_missing(monkeypatch) -> None:
    captured: list[tuple[object | None, object | None, str]] = []

    def fake_loader(filename: str) -> object:
        captured.append((sys.stdout, sys.stderr, filename))
        return object()

    monkeypatch.setattr("zstarview.astro._starfield_load", fake_loader)
    monkeypatch.setattr("zstarview.astro._cached_ephemeris", None)
    monkeypatch.setattr("sys.stdout", None)
    monkeypatch.setattr("sys.stderr", None)

    load_ephemeris()

    assert len(captured) == 1
    stdout_obj, stderr_obj, filename = captured[0]
    assert isinstance(stdout_obj, io.StringIO)
    assert isinstance(stderr_obj, io.StringIO)
    assert filename == "de442s.bsp"
    assert sys.stdout is None
    assert sys.stderr is None


def test_load_ephemeris_caches_loaded_kernel(monkeypatch) -> None:
    calls: list[tuple[object | None, object | None, str]] = []
    sentinel = object()

    def fake_loader(filename: str) -> object:
        calls.append((sys.stdout, sys.stderr, filename))
        return sentinel

    monkeypatch.setattr("zstarview.astro._starfield_load", fake_loader)
    monkeypatch.setattr("zstarview.astro._cached_ephemeris", None)

    first = load_ephemeris()
    second = load_ephemeris()

    assert first is sentinel
    assert second is sentinel
    assert len(calls) == 1


def test_compute_sky_snapshot_uses_provided_ephemeris(monkeypatch) -> None:
    from PySide6.QtGui import QImage

    from zstarview.gui import sky_worker
    from zstarview.paths import THEME_STYLES_BY_PRESET

    ephemeris = object()
    captured: dict[str, object] = {}

    def fake_calculate_visible_stars(*_args, **_kwargs):
        return ({"star_index": []}, object())

    def fake_calculate_visible_deep_sky_objects(*_args, **_kwargs):
        return {"id": [], "name": [], "type": [], "alt": [], "az": [], "vmag": [], "major_arcmin": [], "minor_arcmin": [], "pa_deg": []}

    def fake_calculate_planets(*args, **_kwargs):
        captured["planets"] = args[5]
        return [SimpleNamespace(name="moon", alt=0.0, az=0.0)]

    monkeypatch.setattr(sky_worker, "calculate_visible_stars", fake_calculate_visible_stars)
    monkeypatch.setattr(
        sky_worker,
        "calculate_visible_deep_sky_objects",
        fake_calculate_visible_deep_sky_objects,
    )
    monkeypatch.setattr(sky_worker, "calculate_planets", fake_calculate_planets)
    monkeypatch.setattr(
        sky_worker,
        "calculate_celestial_equator_points",
        lambda *_args, **_kwargs: [],
    )
    monkeypatch.setattr(sky_worker, "calculate_ecliptic_points", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(sky_worker, "calculate_horizon_points", lambda: [])
    monkeypatch.setattr(sky_worker, "compute_night_light_glow_profile", lambda **_kwargs: None)
    monkeypatch.setattr(
        sky_worker.sky_disc,
        "draw_sky_color_disc",
        lambda *args, **kwargs: QImage(),
    )
    monkeypatch.setattr(sky_worker.sky_disc, "draw_uniform_sky_color_disc", lambda *args, **kwargs: QImage())

    viewer_data = ViewerData(
        location=(0.0, 0.0),
        timezone_name="UTC",
        city_name="Test",
        view_center=(0.0, 0.0),
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
        observer_height_m=0.0,
    )
    payload = sky_worker.compute_sky_snapshot(
        ephemeris=ephemeris,
        viewer_data=viewer_data,
        geometry=render_geometry.get_screen_geometry(16, 16, viewer_data.view_alt_deg),
        star_catalog={"catalog_index": []},
        dso_catalog=None,
        star_vmag_limit=None,
        star_subset_indices=None,
        delta_t=timedelta(0),
        sky_disc_alpha=0.0,
        theme=THEME_STYLES_BY_PRESET["night"],
        image_size=(16, 16),
        render_generation=0,
    )

    assert captured["planets"] is ephemeris
    assert "celestial" in payload


def test_compute_sky_snapshot_skips_night_light_without_terrain(monkeypatch) -> None:
    from PySide6.QtGui import QImage

    from zstarview.gui import sky_worker
    from zstarview.paths import THEME_STYLES_BY_PRESET

    ephemeris = object()

    def fake_calculate_visible_stars(*_args, **_kwargs):
        return ({"star_index": []}, object())

    def fake_calculate_visible_deep_sky_objects(*_args, **_kwargs):
        return {"id": [], "name": [], "type": [], "alt": [], "az": [], "vmag": [], "major_arcmin": [], "minor_arcmin": [], "pa_deg": []}

    def fake_calculate_planets(*args, **_kwargs):
        captured_ephemeris = args[5]
        assert captured_ephemeris is ephemeris
        return [SimpleNamespace(name="sun", alt=-5.0, az=0.0, solar_eclipse_info=None)]

    def fake_compute_night_light_glow_profile(**_kwargs):
        raise AssertionError("night light should not compute without terrain")

    monkeypatch.setattr(sky_worker, "calculate_visible_stars", fake_calculate_visible_stars)
    monkeypatch.setattr(
        sky_worker,
        "calculate_visible_deep_sky_objects",
        fake_calculate_visible_deep_sky_objects,
    )
    monkeypatch.setattr(sky_worker, "calculate_planets", fake_calculate_planets)
    monkeypatch.setattr(sky_worker, "calculate_celestial_equator_points", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(sky_worker, "calculate_ecliptic_points", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(sky_worker, "calculate_horizon_points", lambda: [])
    monkeypatch.setattr(sky_worker, "compute_night_light_glow_profile", fake_compute_night_light_glow_profile)
    monkeypatch.setattr(sky_worker.sky_disc, "draw_sky_color_disc", lambda *args, **kwargs: QImage())
    monkeypatch.setattr(sky_worker.sky_disc, "draw_uniform_sky_color_disc", lambda *args, **kwargs: QImage())

    viewer_data = ViewerData(
        location=(0.0, 0.0),
        timezone_name="UTC",
        city_name="Test",
        view_center=(0.0, 0.0),
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
        observer_height_m=0.0,
    )
    payload = sky_worker.compute_sky_snapshot(
        ephemeris=ephemeris,
        viewer_data=viewer_data,
        geometry=render_geometry.get_screen_geometry(16, 16, viewer_data.view_alt_deg),
        star_catalog={"catalog_index": []},
        dso_catalog=None,
        star_vmag_limit=None,
        star_subset_indices=None,
        delta_t=timedelta(0),
        sky_disc_alpha=0.0,
        theme=THEME_STYLES_BY_PRESET["night"],
        image_size=(16, 16),
        render_generation=0,
    )

    assert payload["night_light_glow_profile"] is None


def test_compute_sky_snapshot_passes_absolute_elevation_to_night_light(monkeypatch) -> None:
    from PySide6.QtGui import QImage

    from zstarview.gui import sky_worker
    from zstarview.paths import THEME_STYLES_BY_PRESET

    captured: dict[str, object] = {}

    monkeypatch.setattr(
        sky_worker,
        "calculate_visible_stars",
        lambda *_args, **_kwargs: ({"star_index": []}, object()),
    )
    monkeypatch.setattr(
        sky_worker,
        "calculate_visible_deep_sky_objects",
        lambda *_args, **_kwargs: {
            "id": [],
            "name": [],
            "type": [],
            "alt": [],
            "az": [],
            "vmag": [],
            "major_arcmin": [],
            "minor_arcmin": [],
            "pa_deg": [],
        },
    )
    monkeypatch.setattr(
        sky_worker,
        "calculate_planets",
        lambda *_args, **_kwargs: [
            SimpleNamespace(name="sun", alt=-10.0, az=0.0, solar_eclipse_info=None)
        ],
    )
    monkeypatch.setattr(sky_worker, "calculate_celestial_equator_points", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(sky_worker, "calculate_ecliptic_points", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(sky_worker, "calculate_horizon_points", lambda: [])

    def fake_compute_night_light_glow_profile(**kwargs):
        captured.update(kwargs)
        return object()

    monkeypatch.setattr(
        sky_worker,
        "compute_night_light_glow_profile",
        fake_compute_night_light_glow_profile,
    )
    monkeypatch.setattr(sky_worker.sky_disc, "draw_sky_color_disc", lambda *args, **kwargs: QImage())
    monkeypatch.setattr(
        sky_worker.sky_disc,
        "draw_uniform_sky_color_disc",
        lambda *args, **kwargs: QImage(),
    )

    viewer_data = ViewerData(
        location=(40.0, -106.0),
        timezone_name="UTC",
        city_name="Test",
        view_center=(0.0, 330.0),
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
        observer_height_m=1.7,
        ground_elevation_m=2_885.32,
    )
    payload = sky_worker.compute_sky_snapshot(
        ephemeris=object(),
        viewer_data=viewer_data,
        geometry=render_geometry.get_screen_geometry(16, 16, viewer_data.view_alt_deg),
        star_catalog={"catalog_index": []},
        dso_catalog=None,
        star_vmag_limit=None,
        star_subset_indices=None,
        delta_t=timedelta(0),
        sky_disc_alpha=0.0,
        theme=THEME_STYLES_BY_PRESET["night"],
        image_size=(16, 16),
        terrain_horizon_profile_altaz=[(0.7, 330.0)],
        terrain_horizon_profile_distances_m=[15_000.0],
        terrain_sample_distances_m=np.asarray([15_000.0]),
        terrain_sample_terrain_elevation_m=np.asarray([[3_000.0]]),
        night_light_opacity=0.2,
        render_generation=0,
    )

    assert captured["observer_elevation_m"] == pytest.approx(2_887.02)
    assert payload["night_light_glow_profile"] is not None
