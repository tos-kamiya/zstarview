from __future__ import annotations

from dataclasses import dataclass
from datetime import timedelta
from types import SimpleNamespace

import zstarview.cli.export_image as mod
from zstarview.cli.args import SKY_OPACITY_DEFAULT
from zstarview.gui.window_inputs import SkyWindowRuntimeOptions


@dataclass
class _City:
    display_name: str = "Test City"
    lat: float = 35.0
    lon: float = 139.0
    tz: str = "UTC"
    observer_height_m: float = 1.7
    ground_elevation_m: float | None = 35.0
    location_height_label: str | None = None
    location_height_m: float | None = None


class _Args:
    city = "Test City"
    place = None
    place_countrycode = None
    place_lang = "en"
    timezone = None
    datetime = None
    days = 0
    hours = 0
    vmag_limit = 6.0
    view_center_alt = 90.0
    view_center_az = 180.0
    cloud_stripe = ("width", 50, 0.2)
    theme = "night"
    vmag_brightness_multiplier = 2.5
    content_fov_deg = 100.0
    edge_fov_deg = 95.0
    observer_height_m = None
    search = ""
    list = False
    view_center_alt_specified = False
    view_center_az_specified = False
    output = None
    image_size = (1280, 720)
    layer_timeout_seconds = 30.0
    allow_partial_data = False
    include_direction_grid = False
    window_frame = "frameless"
    sky_opacity = SKY_OPACITY_DEFAULT
    cloud_opacity = 0.15
    satellite_opacity = 0.5
    aircraft_opacity = 0.5
    terrain_horizon_opacity = 0.05
    urban_outline_opacity = 0.2
    ground_tint_opacity = 0.1
    enlarge_moon = False
    star_base_radius = 4.0
    show_dso_initial = None
    show_asterisms_initial = None
    urban_outline_radius_km = 2.5
    urban_outline_skyscraper_radius_km = 60.0
    urban_outline_min_height_m = 0.0
    urban_outline_feature_type = "both"
    urban_outline_skyscraper_only = False
    cloud_missing_tint_opacity = 0.0
    expected_render_width = 600


class _ArgsWithoutWindowFrame(_Args):
    @property
    def window_frame(self) -> str:
        raise AttributeError("window_frame")


def _patch_common(monkeypatch, *, delta_t: timedelta) -> None:
    monkeypatch.setattr(mod, "resolve_launch_location", lambda *args, **kwargs: _City())
    monkeypatch.setattr(mod, "parse_launch_time_arguments", lambda *args, **kwargs: delta_t)
    monkeypatch.setattr(mod, "_load_star_catalog_for_export", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(mod, "_load_dso_catalog_for_export", lambda: None)
    monkeypatch.setattr(mod, "_verify_ephemeris_for_export", lambda: None)
    monkeypatch.setattr(
        mod,
        "prepare_window_catalogs",
        lambda *args, **kwargs: SimpleNamespace(
            star_catalog_np=None,
            star_catalog_lod6_indices=None,
            star_catalog_meta=None,
            dso_catalog_np=None,
        ),
    )
    monkeypatch.setattr(mod, "prepare_window_viewer_data", lambda *args, **kwargs: SimpleNamespace())


def test_build_window_inputs_disables_all_realtime_overlays_for_past(monkeypatch) -> None:
    _patch_common(monkeypatch, delta_t=timedelta(days=-1))

    _catalogs, _viewer_data, user_options, _runtime_options, _search_overlay_target = mod._build_window_inputs_from_args(_Args())

    assert user_options.cloud_disc_alpha == 0.0
    assert user_options.aircraft_opacity == 0.0
    assert user_options.satellite_opacity == 0.0


def test_build_window_inputs_disables_all_realtime_overlays_for_future(monkeypatch) -> None:
    _patch_common(monkeypatch, delta_t=timedelta(days=1))

    _catalogs, _viewer_data, user_options, _runtime_options, _search_overlay_target = mod._build_window_inputs_from_args(_Args())

    assert user_options.cloud_disc_alpha == 0.0
    assert user_options.aircraft_opacity == 0.0
    assert user_options.satellite_opacity == 0.0


def test_build_window_inputs_propagates_cloud_stripe_mode(monkeypatch) -> None:
    _patch_common(monkeypatch, delta_t=timedelta(0))

    args = _Args()
    args.cloud_stripe = ("alpha", 50, 0.2)
    _catalogs, _viewer_data, _user_options, runtime_options, _search_overlay_target = mod._build_window_inputs_from_args(args)

    assert runtime_options.cloud_stripe_mode == "alpha"


def test_build_window_inputs_defaults_export_image_to_frameless(monkeypatch) -> None:
    _patch_common(monkeypatch, delta_t=timedelta(0))

    captured: dict[str, str] = {}
    real_prepare_window_runtime_options = mod.prepare_window_runtime_options

    def _capture_window_frame_mode(*args, **kwargs):
        captured["window_frame_mode"] = kwargs["window_frame_mode"]
        return real_prepare_window_runtime_options(*args, **kwargs)

    monkeypatch.setattr(mod, "prepare_window_runtime_options", _capture_window_frame_mode)

    mod._build_window_inputs_from_args(_ArgsWithoutWindowFrame())

    assert captured["window_frame_mode"] == "frameless"


def test_render_image_draws_direction_grid_when_requested(monkeypatch) -> None:
    scene = SimpleNamespace(
        viewer=SimpleNamespace(
            view_alt_deg=90.0,
            view_center=(0.0, 0.0),
            edge_fov_deg=95.0,
            content_fov_deg=100.0,
        )
    )
    style = SimpleNamespace()
    compositor = SimpleNamespace()
    calls: list[tuple[tuple[object, ...], dict[str, object]]] = []

    monkeypatch.setattr(mod, "render_base_scene_into_painter", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        mod.render_guides,
        "draw_direction_grid_overlay",
        lambda *args, **kwargs: calls.append((args, kwargs)),
    )

    image = mod._render_image(
        image_size=(64, 64),
        scene=scene,
        style=style,
        compositor=compositor,
        draw_direction_grid=True,
    )

    assert image.width() == 64
    assert image.height() == 64
    assert len(calls) == 1


def test_fetch_urban_outline_layer_skips_skyscraper_lookup_when_radius_zero(monkeypatch) -> None:
    viewer_data = SimpleNamespace(lat_deg=35.0, lon_deg=139.0, observer_height_m=1.7)
    runtime_options = SkyWindowRuntimeOptions(
        urban_outline_radius_km=2.5,
        urban_outline_skyscraper_radius_km=0.0,
        urban_outline_feature_type="both",
        urban_outline_min_height_m=0.0,
        urban_outline_skyscraper_only=False,
    )

    monkeypatch.setattr(mod, "_required_feature_types", lambda _mode: ())
    monkeypatch.setattr(
        mod,
        "resolve_overture_release_for_cache_root",
        lambda **_kwargs: None,
    )
    monkeypatch.setattr(
        mod,
        "select_skyscraper_seed_tiles_for_viewer",
        lambda **_kwargs: (_ for _ in ()).throw(AssertionError("skyscraper lookup should be skipped")),
    )

    got = mod._fetch_urban_outline_layer(
        viewer_data=viewer_data,
        runtime_options=runtime_options,
        deadline=None,
    )

    assert got is None
