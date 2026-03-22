from __future__ import annotations

from dataclasses import dataclass
from datetime import timedelta
from types import SimpleNamespace

import zstarview.cli.export_image as mod


@dataclass
class _City:
    display_name: str = "Test City"
    lat: float = 35.0
    lon: float = 139.0
    tz: str = "UTC"
    observer_height_m: float = 1.7
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
    cloud_stripe = (50, 0.2)
    theme = "night"
    vmag_brightness_multiplier = 2.5
    content_fov_deg = 100.0
    observer_height_m = None
    sky_opacity = 0.15
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
    urban_outline_min_height_m = 0.0
    urban_outline_feature_type = "both"
    urban_outline_skyscraper_only = False
    cloud_missing_tint_opacity = 0.0
    expected_render_width = 600


def _patch_common(monkeypatch, *, delta_t: timedelta) -> None:
    monkeypatch.setattr(mod, "resolve_launch_location", lambda *args, **kwargs: _City())
    monkeypatch.setattr(mod, "parse_launch_time_arguments", lambda *args, **kwargs: delta_t)
    monkeypatch.setattr(mod, "_load_star_catalog_for_export", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(mod, "_load_dso_catalog_for_export", lambda: None)
    monkeypatch.setattr(mod, "_verify_ephemeris_for_export", lambda: None)
    monkeypatch.setattr(
        mod,
        "prepare_window_catalogs",
        lambda *args, **kwargs: SimpleNamespace(star_catalog_full_np=None, star_catalog_lod6_np=None, dso_catalog_np=None),
    )
    monkeypatch.setattr(mod, "prepare_window_viewer_data", lambda *args, **kwargs: SimpleNamespace())


def test_build_window_inputs_disables_cloud_and_aircraft_for_past(monkeypatch) -> None:
    _patch_common(monkeypatch, delta_t=timedelta(days=-1))

    _catalogs, _viewer_data, user_options, _runtime_options = mod._build_window_inputs_from_args(_Args())

    assert user_options.cloud_disc_alpha == 0.0
    assert user_options.aircraft_opacity == 0.0
    assert user_options.satellite_opacity == 0.5


def test_build_window_inputs_disables_all_realtime_overlays_for_future(monkeypatch) -> None:
    _patch_common(monkeypatch, delta_t=timedelta(days=1))

    _catalogs, _viewer_data, user_options, _runtime_options = mod._build_window_inputs_from_args(_Args())

    assert user_options.cloud_disc_alpha == 0.0
    assert user_options.aircraft_opacity == 0.0
    assert user_options.satellite_opacity == 0.0
