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
from ..search.models import SearchJumpTarget
from .famous_star_shortcuts import (
    NamedStarShortcut,
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

    sky_disc_alpha: float = 0.1
    sky_disc_style: str = "smooth"
    sky_disc_alt_rings: bool = False
    cloud_disc_alpha: float = 0.075
    satellite_opacity: float = 0.7
    aircraft_opacity: float = 0.4
    terrain_horizon_opacity: float = 0.003
    earth_guide_opacity: float = 0.028
    urban_outline_opacity: float = 0.2
    ground_tint_opacity: float = 0.04
    enlarge_moon: bool = False
    bright_bodies_mode: str = "outline"
    star_base_radius: float = 4.0
    vmag_limit: float = 7.0
    visual_preset: str = "night"
    star_visibility_boost: float = 1.0
    asterism_visibility_boost: float = 1.0
    earth_guide_visibility_boost: float = 1.0
    show_dso_initial: Optional[bool] = None
    show_asterisms_initial: Optional[bool] = None
    show_guidelines_initial: Optional[bool] = None
    # 'auto'|'top'|'bottom'|'off' or None. None means use default behavior (auto).
    observation_info_mode: Optional[str] = None
    sky_disc_gui_allowed: bool = True
    cloud_gui_allowed: bool = True
    satellite_gui_allowed: bool = True
    aircraft_gui_allowed: bool = True
    terrain_horizon_gui_allowed: bool = True
    earth_guide_gui_allowed: bool = True
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
    content_fov_deg: float = 110.0
    window_geometry_arg: Optional[str | tuple[int, int, int, int]] = None
    window_frame_mode: str = "frameless"


def _apply_visibility_boost(opacity: float, visibility_boost: float, tier_scale: float) -> float:
    """Increase opacity by a tiered boost while keeping the value in range."""
    base_opacity = min(1.0, max(0.0, float(opacity)))
    boost = max(1.0, float(visibility_boost))
    if base_opacity <= 0.0 or boost <= 1.0 or tier_scale <= 0.0:
        return base_opacity
    boosted = base_opacity * (1.0 + (boost - 1.0) * float(tier_scale))
    return min(1.0, boosted)


def _apply_visibility_boost_scale(scale: float, visibility_boost: float, tier_scale: float) -> float:
    """Increase a multiplicative scale using the same tiered boost model."""
    base_scale = max(1.0, float(scale))
    boost = max(1.0, float(visibility_boost))
    if boost <= 1.0 or tier_scale <= 0.0:
        return base_scale
    return base_scale * (1.0 + (boost - 1.0) * float(tier_scale))


def prepare_window_viewer_data(
    city_name: str,
    city_data: tuple[float, float, str],
    view_center: tuple[float, float],
    *,
    edge_fov_deg: float,
    content_fov_deg: float,
    observer_height_m: float,
    ground_elevation_m: float,
    location_height_label: str | None,
    location_height_m: float,
    show_observer_height: bool,
) -> ViewerData:
    """Build the viewer input consumed by SkyWindow."""
    lat, lon, tz_name = city_data
    return ViewerData(
        location=(lat, lon),
        timezone_name=tz_name,
        city_name=city_name,
        view_center=view_center,
        edge_fov_deg=float(edge_fov_deg),
        content_fov_deg=float(content_fov_deg),
        observer_height_m=float(observer_height_m),
        ground_elevation_m=max(0.0, float(ground_elevation_m)),
        location_height_label=location_height_label,
        location_height_m=max(0.0, float(location_height_m)),
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
        named_stars_by_band=build_named_star_shortcuts(star_catalog, max_vmag=2.0, include_satellites=False),
        # Search dialog resolves ISS and JPL bodies through their dedicated callbacks.
        named_stars_search_all=build_search_jump_targets(star_catalog, include_satellites=False),
    )


def prepare_window_user_options(
    *,
    sky_disc_alpha: float,
    sky_disc_style: str,
    sky_disc_alt_rings: bool,
    cloud_disc_alpha: float,
    satellite_opacity: float,
    aircraft_opacity: float,
    terrain_horizon_opacity: float,
    earth_guide_opacity: float,
    urban_outline_opacity: float,
    ground_tint_opacity: float,
    enlarge_moon: bool,
    bright_bodies_mode: str,
    star_base_radius: float,
    vmag_limit: float,
    visual_preset: str,
    star_visibility_boost: float,
    visibility_boost: float,
    show_dso_initial: Optional[bool],
    show_asterisms_initial: Optional[bool],
    show_guidelines_initial: Optional[bool],
    observation_info_mode: Optional[str],
    sky_disc_gui_allowed: bool,
    cloud_gui_allowed: bool,
    satellite_gui_allowed: bool,
    aircraft_gui_allowed: bool,
    terrain_horizon_gui_allowed: bool,
    earth_guide_gui_allowed: bool,
    urban_outline_gui_allowed: bool,
) -> SkyWindowUserOptions:
    """Normalize user-facing options before constructing SkyWindow."""
    visibility_boost = max(1.0, float(visibility_boost))
    return SkyWindowUserOptions(
        sky_disc_alpha=_apply_visibility_boost(sky_disc_alpha, visibility_boost, 1.0),
        sky_disc_style=str(sky_disc_style).strip().lower(),
        sky_disc_alt_rings=bool(sky_disc_alt_rings),
        cloud_disc_alpha=_apply_visibility_boost(cloud_disc_alpha, visibility_boost, 1.0),
        satellite_opacity=_apply_visibility_boost(satellite_opacity, visibility_boost, 1.0),
        aircraft_opacity=_apply_visibility_boost(aircraft_opacity, visibility_boost, 1.0),
        terrain_horizon_opacity=_apply_visibility_boost(terrain_horizon_opacity, visibility_boost, 1.0),
        earth_guide_opacity=_apply_visibility_boost(earth_guide_opacity, visibility_boost, 1.0),
        urban_outline_opacity=_apply_visibility_boost(urban_outline_opacity, visibility_boost, 1.0),
        ground_tint_opacity=_apply_visibility_boost(ground_tint_opacity, visibility_boost, 1.0),
        enlarge_moon=bool(enlarge_moon),
        bright_bodies_mode=str(bright_bodies_mode).strip().lower(),
        star_base_radius=max(2.0, star_base_radius),
        vmag_limit=vmag_limit,
        visual_preset=visual_preset,
        star_visibility_boost=star_visibility_boost,
        asterism_visibility_boost=_apply_visibility_boost_scale(1.0, visibility_boost, 1.0),
        earth_guide_visibility_boost=visibility_boost,
        show_dso_initial=show_dso_initial,
        show_asterisms_initial=show_asterisms_initial,
        show_guidelines_initial=show_guidelines_initial,
        observation_info_mode=observation_info_mode,
        sky_disc_gui_allowed=bool(sky_disc_gui_allowed),
        cloud_gui_allowed=bool(cloud_gui_allowed),
        satellite_gui_allowed=bool(satellite_gui_allowed),
        aircraft_gui_allowed=bool(aircraft_gui_allowed),
        terrain_horizon_gui_allowed=bool(terrain_horizon_gui_allowed),
        earth_guide_gui_allowed=bool(earth_guide_gui_allowed),
        urban_outline_gui_allowed=bool(urban_outline_gui_allowed),
    )


def prepare_window_runtime_options(
    *,
    delta_t: timedelta,
    sky_update_interval: int,
    urban_outline_radius_km: float,
    urban_outline_skyscraper_radius_km: float,
    urban_outline_min_height_m: float,
    urban_outline_feature_type: str,
    urban_outline_skyscraper_only: bool,
    cloud_stripe_style: tuple[int, float],
    cloud_stripe_mode: str,
    cloud_missing_tint_opacity: float,
    visibility_boost: float,
    star_render_expected_width: int,
    content_fov_deg: float,
    window_geometry_arg: Optional[str | tuple[int, int, int, int]],
    window_frame_mode: str,
) -> SkyWindowRuntimeOptions:
    """Normalize runtime options before constructing SkyWindow."""
    visibility_boost = max(1.0, float(visibility_boost))
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
        cloud_missing_tint_opacity=_apply_visibility_boost(cloud_missing_tint_opacity, visibility_boost, 1.0),
        star_render_expected_width=max(1, int(star_render_expected_width)),
        content_fov_deg=max(90.0, min(127.0, float(content_fov_deg))),
        window_geometry_arg=window_geometry_arg,
        window_frame_mode=("window" if str(window_frame_mode).strip().lower() == "window" else "frameless"),
    )
