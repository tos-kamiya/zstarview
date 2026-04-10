# -*- coding: utf-8 -*-
"""Prepared input data for SkyWindow construction."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import timedelta
from typing import Optional

import numpy as np
import polars as pl

from ..data.skyscraper_tiles import SKYSCRAPER_OUTER_RADIUS_KM
from ..paths import CLOUD_MISSING_TINT_RGBA
from ..astro import (
    DeepSkyCatalogArrays,
    StarCatalogArrays,
    prepare_star_catalog_meta,
    prepare_deep_sky_catalog_arrays,
    prepare_star_catalog_arrays,
)
from ..types import StarCatalogMeta, ViewerData
from .famous_star_shortcuts import (
    NamedStarShortcut,
    SearchJumpTarget,
    build_named_star_shortcuts,
    build_search_jump_targets,
)


@dataclass(frozen=True)
class PreparedWindowCatalogs:
    """Precomputed catalog inputs consumed by SkyWindow."""

    star_catalog_np: StarCatalogArrays
    star_catalog_lod6_indices: np.ndarray
    star_catalog_meta: StarCatalogMeta
    dso_catalog_np: Optional[DeepSkyCatalogArrays]
    named_stars_by_band: dict[str, list[NamedStarShortcut]]
    named_stars_search_all: list[SearchJumpTarget]


@dataclass(frozen=True)
class SkyWindowUserOptions:
    """User-facing window options that influence rendering and toggles."""

    sky_disc_alpha: float = 0.3
    cloud_disc_alpha: float = 0.6
    satellite_opacity: float = 0.5
    aircraft_opacity: float = 0.6
    terrain_horizon_opacity: float = 0.028
    urban_outline_opacity: float = 0.2
    ground_tint_opacity: float = 1.0
    enlarge_moon: bool = False
    star_base_radius: float = 4.0
    vmag_limit: float = 6.0
    visual_preset: str = "night"
    star_visibility_boost: float = 1.0
    show_dso_initial: Optional[bool] = None
    show_asterisms_initial: Optional[bool] = None
    show_guidelines_initial: Optional[bool] = None
    show_overlay_info_initial: Optional[bool] = None
    sky_disc_gui_allowed: bool = True
    cloud_gui_allowed: bool = True
    satellite_gui_allowed: bool = True
    aircraft_gui_allowed: bool = True
    terrain_horizon_gui_allowed: bool = True
    urban_outline_gui_allowed: bool = True


@dataclass(frozen=True)
class SkyWindowRuntimeOptions:
    """Runtime and window-hosting options for SkyWindow."""

    delta_t: timedelta = timedelta(0)
    sky_update_interval: int = 60
    urban_outline_radius_km: float = 2.5
    urban_outline_skyscraper_radius_km: float = SKYSCRAPER_OUTER_RADIUS_KM
    urban_outline_min_height_m: float = 0.0
    urban_outline_feature_type: str = "both"
    urban_outline_skyscraper_only: bool = False
    cloud_stripe_style: tuple[int, float] = (50, 0.85)
    cloud_stripe_mode: str = "width"
    cloud_missing_tint_opacity: float = float(CLOUD_MISSING_TINT_RGBA[3]) / 255.0
    star_render_expected_width: int = 600
    content_fov_deg: float = 100.0
    window_geometry_arg: Optional[str | tuple[int, int, int, int]] = None
    window_frame_mode: str = "frameless"


def prepare_window_viewer_data(
    city_name: str,
    city_data: tuple[float, float, str],
    view_center: tuple[float, float],
    *,
    content_fov_deg: float = 100.0,
    observer_height_m: float = 1.7,
    location_height_label: str | None = None,
    location_height_m: float | None = None,
    show_observer_height: bool = False,
) -> ViewerData:
    """Build the viewer input consumed by SkyWindow."""
    lat, lon, tz_name = city_data
    return ViewerData(
        location=(lat, lon),
        timezone_name=tz_name,
        city_name=city_name,
        view_center=view_center,
        content_fov_deg=float(content_fov_deg),
        observer_height_m=float(observer_height_m),
        location_height_label=location_height_label,
        location_height_m=None if location_height_m is None else float(location_height_m),
        show_observer_height=bool(show_observer_height),
    )


def prepare_window_catalogs(
    star_catalog: pl.DataFrame,
    dso_catalog: Optional[pl.DataFrame] = None,
    *,
    vmag_brightness_scale: float = -0.39,
) -> PreparedWindowCatalogs:
    """Build prepared SkyWindow catalog inputs from raw data frames."""
    star_catalog_np = prepare_star_catalog_arrays(
        star_catalog,
        vmag_brightness_scale=vmag_brightness_scale,
    )
    return PreparedWindowCatalogs(
        star_catalog_np=star_catalog_np,
        star_catalog_lod6_indices=(star_catalog_np["vmag"] <= 6.0).nonzero()[0].astype("int32", copy=False),
        star_catalog_meta=prepare_star_catalog_meta(star_catalog),
        dso_catalog_np=None if dso_catalog is None else prepare_deep_sky_catalog_arrays(dso_catalog),
        named_stars_by_band=build_named_star_shortcuts(star_catalog, max_vmag=2.0, include_satellites=True),
        named_stars_search_all=build_search_jump_targets(star_catalog, include_satellites=True),
    )


def prepare_window_user_options(
    *,
    sky_disc_alpha: float = 0.3,
    cloud_disc_alpha: float = 0.6,
    satellite_opacity: float = 0.5,
    aircraft_opacity: float = 0.6,
    terrain_horizon_opacity: float = 0.028,
    urban_outline_opacity: float = 0.2,
    ground_tint_opacity: float = 1.0,
    enlarge_moon: bool = False,
    star_base_radius: float = 4.0,
    vmag_limit: float = 6.0,
    visual_preset: str = "night",
    star_visibility_boost: float = 1.0,
    show_dso_initial: Optional[bool] = None,
    show_asterisms_initial: Optional[bool] = None,
    show_guidelines_initial: Optional[bool] = None,
    show_overlay_info_initial: Optional[bool] = None,
    sky_disc_gui_allowed: bool = True,
    cloud_gui_allowed: bool = True,
    satellite_gui_allowed: bool = True,
    aircraft_gui_allowed: bool = True,
    terrain_horizon_gui_allowed: bool = True,
    urban_outline_gui_allowed: bool = True,
) -> SkyWindowUserOptions:
    """Normalize user-facing options before constructing SkyWindow."""
    return SkyWindowUserOptions(
        sky_disc_alpha=min(1.0, max(0.0, sky_disc_alpha)),
        cloud_disc_alpha=min(1.0, max(0.0, cloud_disc_alpha)),
        satellite_opacity=min(1.0, max(0.0, satellite_opacity)),
        aircraft_opacity=min(1.0, max(0.0, aircraft_opacity)),
        terrain_horizon_opacity=min(1.0, max(0.0, terrain_horizon_opacity)),
        urban_outline_opacity=min(1.0, max(0.0, urban_outline_opacity)),
        ground_tint_opacity=min(1.0, max(0.0, ground_tint_opacity)),
        enlarge_moon=bool(enlarge_moon),
        star_base_radius=max(2.0, star_base_radius),
        vmag_limit=vmag_limit,
        visual_preset=visual_preset,
        star_visibility_boost=star_visibility_boost,
        show_dso_initial=show_dso_initial,
        show_asterisms_initial=show_asterisms_initial,
        show_guidelines_initial=show_guidelines_initial,
        show_overlay_info_initial=show_overlay_info_initial,
        sky_disc_gui_allowed=bool(sky_disc_gui_allowed),
        cloud_gui_allowed=bool(cloud_gui_allowed),
        satellite_gui_allowed=bool(satellite_gui_allowed),
        aircraft_gui_allowed=bool(aircraft_gui_allowed),
        terrain_horizon_gui_allowed=bool(terrain_horizon_gui_allowed),
        urban_outline_gui_allowed=bool(urban_outline_gui_allowed),
    )


def prepare_window_runtime_options(
    *,
    delta_t: timedelta = timedelta(0),
    sky_update_interval: int = 60,
    urban_outline_radius_km: float = 2.5,
    urban_outline_skyscraper_radius_km: float = SKYSCRAPER_OUTER_RADIUS_KM,
    urban_outline_min_height_m: float = 0.0,
    urban_outline_feature_type: str = "both",
    urban_outline_skyscraper_only: bool = False,
    cloud_stripe_style: tuple[int, float] = (50, 0.85),
    cloud_stripe_mode: str = "width",
    cloud_missing_tint_opacity: float = float(CLOUD_MISSING_TINT_RGBA[3]) / 255.0,
    star_render_expected_width: int = 600,
    content_fov_deg: float = 100.0,
    window_geometry_arg: Optional[str | tuple[int, int, int, int]] = None,
    window_frame_mode: str = "frameless",
) -> SkyWindowRuntimeOptions:
    """Normalize runtime options before constructing SkyWindow."""
    return SkyWindowRuntimeOptions(
        delta_t=delta_t,
        sky_update_interval=max(1, int(sky_update_interval)),
        urban_outline_radius_km=max(0.0, float(urban_outline_radius_km)),
        urban_outline_skyscraper_radius_km=max(0.0, float(urban_outline_skyscraper_radius_km)),
        urban_outline_min_height_m=max(0.0, float(urban_outline_min_height_m)),
        urban_outline_feature_type=str(urban_outline_feature_type),
        urban_outline_skyscraper_only=bool(urban_outline_skyscraper_only),
        cloud_stripe_style=cloud_stripe_style,
        cloud_stripe_mode=("alpha" if str(cloud_stripe_mode) == "alpha" else "width"),
        cloud_missing_tint_opacity=min(1.0, max(0.0, cloud_missing_tint_opacity)),
        star_render_expected_width=max(1, int(star_render_expected_width)),
        content_fov_deg=max(90.0, min(127.0, float(content_fov_deg))),
        window_geometry_arg=window_geometry_arg,
        window_frame_mode=("window" if str(window_frame_mode).strip().lower() == "window" else "frameless"),
    )
