from __future__ import annotations

import math
import threading
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import zstarview.cli.export_image as mod
from zstarview.cli.args import SKY_OPACITY_DEFAULT
from zstarview.gui.window_inputs import SkyWindowRuntimeOptions
from zstarview.water_mask_interface import WaterSurfaceBandStats
from zstarview.water_overlay import WaterOverlayPoint


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
    clear_long_lived_cache = False
    print_cache_dir = False
    sixel = False
    sky_opacity = SKY_OPACITY_DEFAULT
    sky_disc_style = "smooth"
    sky_disc_altaz_rings = "dimalt"
    sky_disc_altaz_rings_hover = "altaz"
    sky_update_interval = 60
    night_light_opacity = 0.07
    ridge_glow_opacity = 0.04
    cloud_opacity = 0.15
    geo_satellite = False
    satellite_opacity = 0.5
    aircraft_opacity = 0.5
    tropical_cyclone_opacity = 0.4
    precipitation_opacity = 0.6
    terrain_horizon_opacity = 0.05
    earth_guide_opacity = 0.028
    urban_outline_opacity = 0.2
    ground_tint_opacity = 0.1
    overlay_font_size = 11
    bright_bodies = "outline"
    star_base_radius = 4.0
    observation_info = "auto"
    show_dso_initial = None
    show_asterisms_initial = None
    show_guidelines_initial = None
    visibility_boost = 1.0
    urban_outline_radius_km = 2.5
    urban_outline_skyscraper_radius_km = 60.0
    urban_outline_min_height_m = 0.0
    urban_outline_max_candidates = 5000
    urban_outline_download_timeout_seconds = 120.0
    urban_outline_feature_type = "both"
    urban_outline_skyscraper_only = False
    cloud_missing_tint_opacity = 0.0
    water_surface_opacity = 0.4
    expected_render_width = 600


class _ArgsWithoutWindowFrame(_Args):
    @property
    def window_frame(self) -> str:
        raise AttributeError("window_frame")


def _patch_common(monkeypatch, *, delta_t: timedelta) -> None:
    monkeypatch.setattr(mod, "resolve_launch_location", lambda *args, **kwargs: _City())
    monkeypatch.setattr(
        mod, "parse_launch_time_arguments", lambda *args, **kwargs: delta_t
    )
    monkeypatch.setattr(
        mod, "_load_star_catalog_for_export", lambda *_args, **_kwargs: object()
    )
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
    monkeypatch.setattr(
        mod, "prepare_window_viewer_data", lambda *args, **kwargs: SimpleNamespace()
    )


def test_build_window_inputs_disables_all_realtime_overlays_for_past(
    monkeypatch,
) -> None:
    _patch_common(monkeypatch, delta_t=timedelta(days=-1))

    (
        _catalogs,
        _viewer_data,
        user_options,
        _runtime_options,
        _search_overlay_target,
        _place_location,
    ) = mod._build_window_inputs_from_args(_Args())

    assert user_options.cloud_disc_alpha == 0.0
    assert user_options.aircraft_opacity == 0.0
    assert user_options.satellite_opacity == 0.0
    assert user_options.tropical_cyclone_opacity == 0.0
    assert user_options.precipitation_opacity == 0.0
    assert user_options.overlay_font_size == 11


def test_build_window_inputs_disables_all_realtime_overlays_for_future(
    monkeypatch,
) -> None:
    _patch_common(monkeypatch, delta_t=timedelta(days=1))

    (
        _catalogs,
        _viewer_data,
        user_options,
        _runtime_options,
        _search_overlay_target,
        _place_location,
    ) = mod._build_window_inputs_from_args(_Args())

    assert user_options.cloud_disc_alpha == 0.0
    assert user_options.aircraft_opacity == 0.0
    assert user_options.satellite_opacity == 0.0
    assert user_options.tropical_cyclone_opacity == 0.0
    assert user_options.precipitation_opacity == 0.0


def test_build_window_inputs_propagates_cloud_stripe_mode(monkeypatch) -> None:
    _patch_common(monkeypatch, delta_t=timedelta(0))

    args = _Args()
    args.cloud_stripe = ("alpha", 50, 0.2)
    (
        _catalogs,
        _viewer_data,
        _user_options,
        runtime_options,
        _search_overlay_target,
        _place_location,
    ) = mod._build_window_inputs_from_args(args)

    assert runtime_options.cloud_stripe_mode == "alpha"


def test_build_window_inputs_defaults_export_image_to_frameless(monkeypatch) -> None:
    _patch_common(monkeypatch, delta_t=timedelta(0))

    captured: dict[str, str] = {}
    real_prepare_window_runtime_options = mod.prepare_window_runtime_options

    def _capture_window_frame_mode(*args, **kwargs):
        captured["window_frame_mode"] = kwargs["window_frame_mode"]
        return real_prepare_window_runtime_options(*args, **kwargs)

    monkeypatch.setattr(
        mod, "prepare_window_runtime_options", _capture_window_frame_mode
    )

    mod._build_window_inputs_from_args(_ArgsWithoutWindowFrame())

    assert captured["window_frame_mode"] == "frameless"


def test_render_image_draws_direction_grid_when_requested(monkeypatch) -> None:
    scene = SimpleNamespace(
        viewer=SimpleNamespace(
            view_alt_deg=90.0,
            view_center=(0.0, 0.0),
            edge_fov_deg=95.0,
            content_fov_deg=100.0,
        ),
        time_obj=None,
    )
    style = SimpleNamespace()
    compositor = SimpleNamespace()
    calls: list[tuple[tuple[object, ...], dict[str, object]]] = []

    monkeypatch.setattr(
        mod, "render_base_scene_into_painter", lambda *args, **kwargs: None
    )
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


def test_render_image_places_time_of_day_marker_in_top_left(monkeypatch) -> None:
    scene = SimpleNamespace(
        viewer=SimpleNamespace(
            view_alt_deg=90.0,
            view_center=(0.0, 0.0),
            edge_fov_deg=95.0,
            content_fov_deg=100.0,
        ),
        time_obj=None,
    )
    captured: dict[str, object] = {}

    def _capture_render(*args, **kwargs) -> None:
        captured["hud"] = kwargs["hud"]

    monkeypatch.setattr(mod, "render_base_scene_into_painter", _capture_render)

    mod._render_image(
        image_size=(64, 64),
        scene=scene,
        style=SimpleNamespace(),
        compositor=SimpleNamespace(),
        draw_direction_grid=False,
    )

    assert captured["hud"].time_of_day_marker_bottom_left is False


def test_fetch_urban_outline_layer_skips_skyscraper_lookup_when_radius_zero(
    monkeypatch,
) -> None:
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
        lambda **_kwargs: (_ for _ in ()).throw(
            AssertionError("skyscraper lookup should be skipped")
        ),
    )

    got = mod._fetch_urban_outline_layer(
        viewer_data=viewer_data,
        runtime_options=runtime_options,
        deadline=None,
    )

    assert got.outlines is None
    assert got.source == "Overture Maps"


def test_fetch_urban_outline_layer_uses_plateau_without_overture(monkeypatch) -> None:
    viewer_data = SimpleNamespace(
        lat_deg=35.455,
        lon_deg=133.055,
        view_center=(45.0, 180.0),
        content_fov_deg=90.0,
    )
    runtime_options = SkyWindowRuntimeOptions(
        urban_outline_radius_km=2.5,
        urban_outline_skyscraper_radius_km=60.0,
        urban_outline_feature_type="both",
        urban_outline_min_height_m=0.0,
        urban_outline_skyscraper_only=False,
    )

    monkeypatch.setattr(
        mod,
        "select_prepared_building_source",
        lambda **_kwargs: SimpleNamespace(
            source="plateau", derived_dirs=(Path("plateau") / "bldg",)
        ),
    )
    monkeypatch.setattr(
        mod,
        "resolve_urban_outline_layer_for_viewer",
        lambda *_args, **kwargs: kwargs["derived_dirs"],
    )
    monkeypatch.setattr(
        mod,
        "resolve_overture_release_for_cache_root",
        lambda **_kwargs: (_ for _ in ()).throw(
            AssertionError("Overture release lookup should be skipped")
        ),
    )

    result = mod._fetch_urban_outline_layer(
        viewer_data=viewer_data,
        runtime_options=runtime_options,
        deadline=None,
    )

    assert result.outlines == (Path("plateau") / "bldg",)
    assert result.source == "PLATEAU"


def test_fetch_water_overlay_layer_uses_observer_ground_and_eye_height(
    monkeypatch, caplog
) -> None:
    viewer_data = SimpleNamespace(
        lat_deg=35.0,
        lon_deg=139.0,
        observer_height_m=1.7,
        ground_elevation_m=42.0,
        view_center=(45.0, 180.0),
        content_fov_deg=110.0,
    )
    captured: dict[str, float] = {}

    def _sample_water_surface_interface_points(*_args, **kwargs):
        captured["observer_height_m"] = float(kwargs["observer_height_m"])
        captured["max_distance_km"] = float(kwargs["max_distance_km"])
        captured["azimuth_step_deg"] = float(kwargs["azimuth_step_deg"])
        return (
            (WaterOverlayPoint("water", 10.0, 20.0, 0.5, water_category="sea-500"),),
            (
                WaterSurfaceBandStats("125m", 0, 0, 0, 0),
                WaterSurfaceBandStats("250m", 0, 0, 0, 0),
                WaterSurfaceBandStats("500m", 1, 9, 1, 1),
            ),
        )

    monkeypatch.setattr(
        mod,
        "sample_water_surface_interface_points_with_stats",
        _sample_water_surface_interface_points,
    )
    monkeypatch.setattr(
        mod,
        "sample_water_overlay_points_for_observer",
        lambda **_kwargs: (),
    )
    monkeypatch.setattr(mod, "fetch_water_overlay_footprints", lambda **_kwargs: ())

    with caplog.at_level("INFO", logger="zstarview.cli.export_image"):
        got = mod._fetch_water_overlay_dots_layer(
            viewer_data=viewer_data,
            surface_size_px=(1280, 720),
            deadline=None,
        )

    assert got == [
        WaterOverlayPoint("water", 10.0, 20.0, 0.5, water_category="sea-500")
    ]
    assert captured["observer_height_m"] == 43.7
    assert captured["azimuth_step_deg"] == 2.0
    assert captured["max_distance_km"] == mod.resolve_water_scan_radius_km(
        43.7,
        minimum_distance_km=mod.DEFAULT_WATER_RADIUS_KM,
    )
    assert "Water band stats: 500m tiles=1 raw=9 collapsed=1 visible=1" in caplog.text
    assert (
        "Water mask dots: 1 visible, nearest sea dot 0.500 km, bands: 125m=0 250m=0 500m=1"
        in caplog.text
    )


def test_fetch_water_overlay_layer_passes_target_ground_sampler(monkeypatch) -> None:
    viewer_data = SimpleNamespace(
        lat_deg=35.0,
        lon_deg=139.0,
        observer_height_m=1.7,
        ground_elevation_m=42.0,
        view_center=(45.0, 180.0),
        content_fov_deg=110.0,
    )
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        mod,
        "sample_water_surface_interface_points_with_stats",
        lambda **_kwargs: (
            (WaterOverlayPoint("water", 10.0, 20.0, 0.5, water_category="sea-500"),),
            (
                WaterSurfaceBandStats("125m", 0, 0, 0, 0),
                WaterSurfaceBandStats("250m", 0, 0, 0, 0),
                WaterSurfaceBandStats("500m", 1, 9, 1, 1),
            ),
        ),
    )
    monkeypatch.setattr(
        mod, "_load_or_fetch_water_overlay_footprints", lambda **_kwargs: ("footprint",)
    )

    def _sample_water_overlay_points_for_observer(*_args, **kwargs):
        captured["target_ground_elevation_m_sampler"] = kwargs[
            "target_ground_elevation_m_sampler"
        ]
        sampler = kwargs["target_ground_elevation_m_sampler"]
        if sampler is None:
            raise AssertionError("expected a DEM sampler")
        captured["sampler_value"] = float(sampler(35.0, 139.0))
        return (WaterOverlayPoint("inland", 11.0, 314.0, 1.5, water_category="lake"),)

    monkeypatch.setattr(
        mod,
        "sample_water_overlay_points_for_observer",
        _sample_water_overlay_points_for_observer,
    )

    got = mod._fetch_water_overlay_dots_layer(
        viewer_data=viewer_data,
        surface_size_px=(1280, 720),
        deadline=None,
        target_ground_sampler=lambda *_args: 77.0,
    )

    assert got == [
        WaterOverlayPoint("water", 10.0, 20.0, 0.5, water_category="sea-500"),
        WaterOverlayPoint("inland", 11.0, 314.0, 1.5, water_category="lake"),
    ]
    assert captured["target_ground_elevation_m_sampler"] is not None
    assert captured["sampler_value"] == 77.0


def test_fetch_cloud_layer_skips_clouds_in_supported_band_without_geo_satellite(
    monkeypatch,
) -> None:
    viewer_data = SimpleNamespace(
        lat_deg=51.5,
        lon_deg=-0.1,
        view_alt_deg=12.0,
        view_az_deg=180.0,
        edge_fov_deg=60.0,
    )
    user_options = SimpleNamespace(
        cloud_disc_alpha=0.2,
        geo_satellite=False,
    )

    warnings: list[str] = []
    monkeypatch.setattr(mod, "is_within_europe_band", lambda *_args: True)
    monkeypatch.setattr(
        mod.logger,
        "warning",
        lambda message, *args, **kwargs: warnings.append(
            message % args if args else message
        ),
    )

    (
        cloud_rgba,
        missing_mask,
        cloud_amount_field,
        cloud_coverage_ratio,
        cloud_altaz_grid,
    ) = mod._fetch_cloud_layer(
        viewer_data=viewer_data,
        user_options=user_options,
        deadline=None,
    )

    assert cloud_rgba is None
    assert missing_mask is None
    assert cloud_amount_field is None
    assert cloud_coverage_ratio is None
    assert cloud_altaz_grid is None
    assert warnings == [
        "Cloud rendering is unavailable for this Europe-band location without --geo-satellite true; skipping the cloud layer."
    ]


def test_fetch_cloud_layer_uses_geo_satellite_branch_when_enabled(monkeypatch) -> None:
    viewer_data = SimpleNamespace(
        lat_deg=51.5,
        lon_deg=-0.1,
        view_alt_deg=12.0,
        view_az_deg=180.0,
        edge_fov_deg=60.0,
    )
    user_options = SimpleNamespace(
        cloud_disc_alpha=0.2,
        geo_satellite=True,
    )

    calls: dict[str, object] = {}
    timeout_checks = [False, True]
    monkeypatch.setattr(
        mod,
        "_timed_out",
        lambda _deadline: timeout_checks.pop(0) if timeout_checks else False,
    )
    monkeypatch.setattr(mod, "is_within_europe_band", lambda *_args: True)
    monkeypatch.setattr(
        mod,
        "run_geo_satellite_pipeline",
        lambda **kwargs: calls.setdefault(
            "pipeline",
            SimpleNamespace(
                download=SimpleNamespace(
                    fetched_at_utc=datetime(2026, 5, 26, tzinfo=timezone.utc)
                ),
                disc_gray=np.full((4, 4), 180, dtype=np.uint8),
                altaz_grid=SimpleNamespace(tag="grid"),
            ),
        ),
    )
    monkeypatch.setattr(
        mod,
        "render_gray_image_to_cloud_rgba",
        lambda gray: np.dstack(
            [
                np.full((*gray.shape, 3), 255, dtype=np.uint8),
                gray.astype(np.uint8)[..., None],
            ]
        ),
    )
    (
        cloud_rgba,
        missing_mask,
        cloud_amount_field,
        cloud_coverage_ratio,
        cloud_altaz_grid,
    ) = mod._fetch_cloud_layer(
        viewer_data=viewer_data,
        user_options=user_options,
        deadline=None,
    )

    assert "pipeline" in calls
    assert cloud_rgba.shape == (4, 4, 4)
    assert missing_mask is None
    assert cloud_amount_field is None
    assert cloud_coverage_ratio == pytest.approx(1.0)
    assert cloud_altaz_grid.tag == "grid"
    assert timeout_checks == [True]


def test_fetch_terrain_horizon_layer_uses_sea_level_fallback(monkeypatch) -> None:
    viewer_data = SimpleNamespace(
        lat_deg=35.0,
        lon_deg=139.0,
        observer_height_m=1.7,
    )

    def _raise_no_tiles(**_kwargs):
        raise RuntimeError(
            "No Copernicus DEM tiles were downloaded for the requested area."
        )

    monkeypatch.setattr(mod, "fetch_copernicus_dem", _raise_no_tiles)

    got = mod._fetch_terrain_horizon_layer(viewer_data=viewer_data, deadline=None)

    assert len(got["profile_altaz"]) == 720
    assert len(got["profile_distances_m"]) == 720
    assert len(got["secondary_ridges_altaz_layers"]) >= 1
    assert len(got["secondary_ridges_distances_m_layers"]) >= 1
    assert all(math.isfinite(alt) for alt, _az in got["profile_altaz"])
    assert min(got["profile_distances_m"]) > 0.0
    assert max(got["profile_distances_m"]) == pytest.approx(
        min(got["profile_distances_m"])
    )


def test_main_uses_independent_layer_deadlines(monkeypatch) -> None:
    deadline_calls = 0
    real_deadline_after = mod._deadline_after
    utc_timezone = timezone.utc

    def _counting_deadline_after(timeout_seconds: float) -> float | None:
        nonlocal deadline_calls
        deadline_calls += 1
        return real_deadline_after(timeout_seconds)

    class _DummyTime:
        def to_datetime(self, timezone=None):
            return datetime(2026, 5, 26, tzinfo=utc_timezone)

    catalogs = SimpleNamespace(
        star_catalog_np=object(),
        star_catalog_lod6_indices=object(),
        star_catalog_meta=None,
        dso_catalog_np=None,
    )
    viewer = SimpleNamespace(
        lat_deg=51.5,
        lon_deg=-0.1,
        observer_height_m=1.7,
        ground_elevation_m=0.0,
        timezone_name="UTC",
        city_name="London",
        view_center=(12.0, 180.0),
        view_alt_deg=12.0,
        view_az_deg=180.0,
        content_fov_deg=100.0,
        height_add_m=1.7,
        location_height_label=None,
        location_height_m=0.0,
        edge_fov_deg=60.0,
    )
    user_options = SimpleNamespace(
        vmag_limit=6.0,
        sky_disc_alpha=0.05,
        sky_disc_style="smooth",
        cloud_disc_alpha=0.2,
        geo_satellite=True,
        water_overlay_opacity=0.4,
        satellite_opacity=0.2,
        terrain_horizon_opacity=0.05,
        urban_outline_opacity=0.0,
        aircraft_opacity=0.2,
        overlay_font_size=11,
        visual_preset="night",
        cloud_missing_tint_opacity=0.0,
    )
    runtime_options = SimpleNamespace(
        delta_t=0.0,
        star_render_expected_width=600,
        urban_outline_skyscraper_only=False,
        urban_outline_feature_type="both",
        urban_outline_radius_km=2.5,
        urban_outline_skyscraper_radius_km=60.0,
        urban_outline_min_height_m=0.0,
    )

    monkeypatch.setattr(
        mod,
        "parse_export_image_args",
        lambda: SimpleNamespace(
            city="London",
            place=None,
            place_countrycode=None,
            place_lang="en",
            timezone=None,
            datetime=None,
            days=0,
            hours=0,
            vmag_limit=6.0,
            view_center_alt=90.0,
            view_center_az=180.0,
            cloud_stripe=("width", 50, 0.2),
            theme="night",
            vmag_brightness_multiplier=2.5,
            content_fov_deg=100.0,
            edge_fov_deg=95.0,
            observer_height_m=None,
            search="",
            list=False,
            view_center_alt_specified=False,
            view_center_az_specified=False,
            output="-",
            image_size=(64, 64),
            layer_timeout_seconds=30.0,
            allow_partial_data=False,
            include_direction_grid=False,
            window_frame="frameless",
            clear_long_lived_cache=False,
            print_cache_dir=False,
            sixel=False,
            sky_opacity=SKY_OPACITY_DEFAULT,
            sky_disc_style="smooth",
            sky_disc_altaz_rings="dimalt",
            sky_disc_altaz_rings_hover="altaz",
            sky_update_interval=60,
            night_light_opacity=0.0,
            ridge_glow_opacity=0.04,
            cloud_opacity=0.15,
            geo_satellite=True,
            satellite_opacity=0.2,
            aircraft_opacity=0.2,
            terrain_horizon_opacity=0.05,
            earth_guide_opacity=0.028,
            urban_outline_opacity=0.0,
            ground_tint_opacity=0.1,
            overlay_font_size=11,
            bright_bodies="outline",
            star_base_radius=4.0,
            observation_info="auto",
            show_dso_initial=None,
            show_asterisms_initial=None,
            show_guidelines_initial=None,
            visibility_boost=1.0,
            urban_outline_radius_km=2.5,
            urban_outline_skyscraper_radius_km=60.0,
            urban_outline_min_height_m=0.0,
            urban_outline_feature_type="both",
            urban_outline_skyscraper_only=False,
            cloud_missing_tint_opacity=0.0,
            water_surface_opacity=0.4,
            expected_render_width=600,
        ),
    )
    monkeypatch.setattr(mod, "setup_root_logger", lambda: None)
    monkeypatch.setattr(mod, "_deadline_after", _counting_deadline_after)
    monkeypatch.setattr(
        mod,
        "setup_app",
        lambda _name: SimpleNamespace(setQuitOnLastWindowClosed=lambda _flag: None),
    )
    monkeypatch.setattr(mod, "_load_fonts", lambda *_args: (object(), object()))
    monkeypatch.setattr(mod, "_build_compositor", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(
        mod,
        "_build_window_inputs_from_args",
        lambda _args: (catalogs, viewer, user_options, runtime_options, None, None),
    )
    monkeypatch.setattr(mod, "load_ephemeris", lambda: object())
    monkeypatch.setattr(
        mod,
        "compute_sky_snapshot",
        lambda **_kwargs: {
            "celestial": SimpleNamespace(
                time=_DummyTime(),
                planets=[SimpleNamespace(name="sun", alt=-10.0)],
            ),
            "sky_disc": None,
        },
    )
    monkeypatch.setattr(mod, "is_within_europe_band", lambda *_args: True)
    monkeypatch.setattr(
        mod,
        "_fetch_cloud_layer",
        lambda **_kwargs: (
            np.zeros((2, 2, 4), dtype=np.uint8),
            None,
            SimpleNamespace(),
            1.0,
        ),
    )
    monkeypatch.setattr(
        mod,
        "_fetch_terrain_horizon_layer",
        lambda **_kwargs: {
            "profile_altaz": [(0.0, 0.0)],
            "profile_distances_m": [0.0],
            "secondary_ridges_altaz_layers": [],
            "secondary_ridges_distances_m_layers": [],
        },
    )
    monkeypatch.setattr(
        mod, "_build_water_target_ground_sampler", lambda **_kwargs: lambda *_args: 0.0
    )
    monkeypatch.setattr(mod, "_fetch_water_overlay_dots_layer", lambda **_kwargs: [])
    monkeypatch.setattr(mod, "_fetch_aircraft_snapshots", lambda **_kwargs: [])
    monkeypatch.setattr(mod, "_fetch_satellite_records_by_group", lambda **_kwargs: {})
    monkeypatch.setattr(
        mod,
        "_build_render_style",
        lambda **_kwargs: SimpleNamespace(vmag_limit=6.0),
    )
    monkeypatch.setattr(
        mod, "_write_export_overlay_summary_to_stderr", lambda **_kwargs: None
    )
    monkeypatch.setattr(mod, "_render_image", lambda **_kwargs: SimpleNamespace())
    monkeypatch.setattr(mod, "_write_png_to_stdout", lambda _image, **_kwargs: True)

    mod.main()

    assert deadline_calls == 6


def test_main_parallelizes_independent_export_layers(monkeypatch) -> None:
    _patch_common(monkeypatch, delta_t=timedelta(0))

    phase1_release = threading.Event()
    phase2_release = threading.Event()
    started: list[str] = []
    started_events = {
        "aircraft": threading.Event(),
        "satellite": threading.Event(),
        "terrain": threading.Event(),
        "urban": threading.Event(),
        "nightlight": threading.Event(),
        "water": threading.Event(),
    }
    lock = threading.Lock()
    night_light_kwargs: dict[str, object] = {}

    def _record_start(name: str) -> None:
        with lock:
            started.append(name)
        started_events[name].set()

    def _make_blocking_task(
        name: str,
        result: object,
        *,
        release_event: threading.Event,
    ):
        def _task(**_kwargs):
            _record_start(name)
            if not release_event.wait(timeout=2.0):
                raise RuntimeError(f"{name} task did not receive release signal")
            return result

        return _task

    catalogs = SimpleNamespace(
        star_catalog_np=object(),
        star_catalog_lod6_indices=object(),
        star_catalog_meta=None,
        dso_catalog_np=None,
    )
    viewer = SimpleNamespace(
        lat_deg=51.5,
        lon_deg=-0.1,
        observer_height_m=1.7,
        ground_elevation_m=35.0,
        timezone_name="UTC",
        city_name="London",
        view_alt_deg=12.0,
        view_az_deg=180.0,
        edge_fov_deg=60.0,
        content_fov_deg=100.0,
        view_center=(0.0, 0.0),
        height_add_m=1.7,
        location_height_label=None,
        location_height_m=0.0,
    )
    user_options = SimpleNamespace(
        vmag_limit=6.0,
        sky_disc_alpha=0.05,
        sky_disc_style="smooth",
        cloud_disc_alpha=0.0,
        geo_satellite=False,
        satellite_opacity=0.2,
        water_overlay_opacity=0.4,
        aircraft_opacity=0.2,
        terrain_horizon_opacity=0.05,
        urban_outline_opacity=0.2,
        night_light_opacity=0.07,
        overlay_font_size=11,
        visual_preset="night",
        cloud_missing_tint_opacity=0.0,
    )
    runtime_options = SimpleNamespace(
        delta_t=0.0,
        star_render_expected_width=600,
        urban_outline_skyscraper_only=False,
        urban_outline_feature_type="both",
        urban_outline_radius_km=2.5,
        urban_outline_skyscraper_radius_km=60.0,
        urban_outline_min_height_m=0.0,
        urban_outline_max_candidates=5000,
    )

    class _DummyTime:
        def to_datetime(self, timezone=None):
            return datetime(2026, 5, 26, tzinfo=timezone)

    monkeypatch.setattr(
        mod,
        "parse_export_image_args",
        lambda: SimpleNamespace(
            city="London",
            place=None,
            place_countrycode=None,
            place_lang="en",
            timezone=None,
            datetime=None,
            days=0,
            hours=0,
            vmag_limit=6.0,
            view_center_alt=90.0,
            view_center_az=180.0,
            cloud_stripe=("width", 50, 0.2),
            theme="night",
            vmag_brightness_multiplier=2.5,
            content_fov_deg=100.0,
            edge_fov_deg=95.0,
            observer_height_m=None,
            search="",
            list=False,
            view_center_alt_specified=False,
            view_center_az_specified=False,
            output="-",
            image_size=(64, 64),
            layer_timeout_seconds=30.0,
            allow_partial_data=False,
            include_direction_grid=False,
            window_frame="frameless",
            clear_long_lived_cache=False,
            print_cache_dir=False,
            sixel=False,
            sky_opacity=SKY_OPACITY_DEFAULT,
            sky_disc_style="smooth",
            sky_disc_altaz_rings="dimalt",
            sky_disc_altaz_rings_hover="altaz",
            sky_update_interval=60,
            night_light_opacity=0.07,
            ridge_glow_opacity=0.04,
            cloud_opacity=0.0,
            geo_satellite=False,
            satellite_opacity=0.2,
            aircraft_opacity=0.2,
            terrain_horizon_opacity=0.05,
            earth_guide_opacity=0.028,
            urban_outline_opacity=0.2,
            ground_tint_opacity=0.1,
            overlay_font_size=11,
            bright_bodies="outline",
            star_base_radius=4.0,
            observation_info="auto",
            show_dso_initial=None,
            show_asterisms_initial=None,
            show_guidelines_initial=None,
            visibility_boost=1.0,
            urban_outline_radius_km=2.5,
            urban_outline_skyscraper_radius_km=60.0,
            urban_outline_min_height_m=0.0,
            urban_outline_feature_type="both",
            urban_outline_skyscraper_only=False,
            cloud_missing_tint_opacity=0.0,
            water_surface_opacity=0.4,
            expected_render_width=600,
        ),
    )
    monkeypatch.setattr(mod, "setup_root_logger", lambda: None)
    monkeypatch.setattr(
        mod,
        "setup_app",
        lambda _name: SimpleNamespace(setQuitOnLastWindowClosed=lambda _flag: None),
    )
    monkeypatch.setattr(mod, "_load_fonts", lambda *_args: (object(), object()))
    monkeypatch.setattr(mod, "_build_compositor", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(
        mod,
        "_build_window_inputs_from_args",
        lambda _args: (catalogs, viewer, user_options, runtime_options, None, None),
    )
    monkeypatch.setattr(mod, "load_ephemeris", lambda: object())
    monkeypatch.setattr(
        mod,
        "compute_sky_snapshot",
        lambda **_kwargs: {
            "celestial": SimpleNamespace(
                time=_DummyTime(),
                planets=[SimpleNamespace(name="sun", alt=-10.0)],
            ),
            "sky_disc": None,
        },
    )
    monkeypatch.setattr(
        mod,
        "_fetch_aircraft_snapshots",
        _make_blocking_task(
            "aircraft",
            [],
            release_event=phase1_release,
        ),
    )
    monkeypatch.setattr(
        mod,
        "_fetch_satellite_records_by_group",
        _make_blocking_task(
            "satellite",
            {},
            release_event=phase1_release,
        ),
    )
    monkeypatch.setattr(
        mod,
        "_fetch_terrain_horizon_layer",
        _make_blocking_task(
            "terrain",
            {
                "profile_altaz": [(0.0, 0.0)],
                "profile_distances_m": [0.0],
                "secondary_ridges_altaz_layers": [],
                "secondary_ridges_distances_m_layers": [],
                "sample_distances_m": None,
                "sample_terrain_elevation_m": None,
                "ground_elevation_m": 42.0,
            },
            release_event=phase1_release,
        ),
    )
    monkeypatch.setattr(
        mod,
        "_fetch_urban_outline_layer",
        _make_blocking_task(
            "urban",
            [],
            release_event=phase2_release,
        ),
    )
    def _compute_night_light(**kwargs):
        night_light_kwargs.update(kwargs)
        return _make_blocking_task(
            "nightlight",
            SimpleNamespace(),
            release_event=phase2_release,
        )()

    monkeypatch.setattr(mod, "compute_night_light_glow_profile", _compute_night_light)
    monkeypatch.setattr(
        mod,
        "_fetch_water_overlay_dots_layer",
        _make_blocking_task(
            "water",
            [],
            release_event=phase2_release,
        ),
    )
    monkeypatch.setattr(
        mod,
        "_build_render_style",
        lambda **_kwargs: SimpleNamespace(vmag_limit=6.0),
    )
    monkeypatch.setattr(
        mod, "_write_export_overlay_summary_to_stderr", lambda **_kwargs: None
    )
    monkeypatch.setattr(mod, "_render_image", lambda **_kwargs: SimpleNamespace())
    monkeypatch.setattr(mod, "_write_png_to_stdout", lambda _image, **_kwargs: True)

    main_thread = threading.Thread(target=mod.main)
    main_thread.start()

    assert started_events["aircraft"].wait(timeout=2.0)
    assert started_events["satellite"].wait(timeout=2.0)
    assert started_events["terrain"].wait(timeout=2.0)
    assert set(started) == {"aircraft", "satellite", "terrain"}

    phase1_release.set()

    assert started_events["urban"].wait(timeout=2.0)
    assert started_events["nightlight"].wait(timeout=2.0)
    assert started_events["water"].wait(timeout=2.0)
    assert started.index("terrain") < started.index("urban")
    assert started.index("terrain") < started.index("nightlight")
    assert started.index("terrain") < started.index("water")

    phase2_release.set()
    main_thread.join(timeout=2.0)
    assert not main_thread.is_alive()
    assert night_light_kwargs["observer_elevation_m"] == pytest.approx(43.7)
