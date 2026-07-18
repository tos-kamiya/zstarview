# -*- coding: utf-8 -*-
"""Prepared input data for SkyWindow construction."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import timedelta
from typing import Callable, Optional

import numpy as np
import polars as pl

from ..astro import (
    DeepSkyCatalogArrays,
    StarCatalogArrays,
    prepare_deep_sky_catalog_arrays,
    prepare_star_catalog_arrays,
    prepare_star_catalog_meta,
)
from ..data.skyscraper_tiles import SKYSCRAPER_OUTER_RADIUS_KM
from ..data.import_overture_buildings import DEFAULT_DOWNLOAD_TIMEOUT_SECONDS
from ..data.urban_outline_from_buildings import MAX_URBAN_OUTLINE_CANDIDATES
from ..paths import (
    CLOUD_DEFAULT_OPACITY,
    CLOUD_MISSING_TINT_RGBA,
    NIGHT_LIGHT_DEFAULT_OPACITY,
    OVERLAY_FONT_SIZE_DEFAULT,
    OVERLAY_FONT_SIZE_MAX,
    OVERLAY_FONT_SIZE_MIN,
    RIDGE_GLOW_DEFAULT_OPACITY,
    TROPICAL_CYCLONE_DEFAULT_OPACITY,
)
from ..search.models import SearchJumpTarget
from ..types import StarCatalogMeta, ViewerData
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

    presentation_id: str = "scenic"
    star_data_policy: str = "scenic_view_scoped"
    sky_disc_alpha: float = 0.15
    sky_disc_style: str = "smooth"
    sky_disc_altaz_rings: str = "dimalt"
    sky_disc_altaz_rings_hover: str = "altaz"
    night_light_opacity: float = NIGHT_LIGHT_DEFAULT_OPACITY
    ridge_glow_opacity: float = RIDGE_GLOW_DEFAULT_OPACITY
    cloud_disc_alpha: float = CLOUD_DEFAULT_OPACITY
    geo_satellite: bool = False
    satellite_opacity: float = 0.7
    aircraft_opacity: float = 0.4
    tropical_cyclone_opacity: float = TROPICAL_CYCLONE_DEFAULT_OPACITY
    terrain_horizon_opacity: float = 0.003
    earth_guide_opacity: float = 0.028
    urban_outline_opacity: float = 0.2
    water_overlay_opacity: float = 0.4
    ground_tint_opacity: float = 0.025
    overlay_font_size: float = float(OVERLAY_FONT_SIZE_DEFAULT)
    enlarge_moon: bool = False
    bright_bodies_mode: str = "outline"
    star_base_radius: float = 4.0
    vmag_limit: float = 7.0
    visual_preset: str = "night"
    star_visibility_boost: float = 1.0
    light_background_star_outline: bool = False
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
    tropical_cyclone_gui_allowed: bool = True
    terrain_horizon_gui_allowed: bool = True
    earth_guide_gui_allowed: bool = True
    night_light_gui_allowed: bool = True
    urban_outline_gui_allowed: bool = True


@dataclass(frozen=True)
class SkyWindowRuntimeOptions:
    """Runtime and window-hosting options for SkyWindow."""

    delta_t: timedelta = timedelta(0)
    sky_update_interval: int = 60
    urban_outline_radius_km: float = 2.5
    urban_outline_skyscraper_radius_km: float = SKYSCRAPER_OUTER_RADIUS_KM
    urban_outline_min_height_m: float = 0.0
    urban_outline_max_candidates: int = MAX_URBAN_OUTLINE_CANDIDATES
    urban_outline_feature_type: str = "both"
    urban_outline_skyscraper_only: bool = False
    urban_outline_download_timeout_seconds: float = DEFAULT_DOWNLOAD_TIMEOUT_SECONDS
    cloud_stripe_style: tuple[int, float] = (30, 1.7)
    cloud_stripe_mode: str = "halftone"
    cloud_missing_tint_opacity: float = float(CLOUD_MISSING_TINT_RGBA[3]) / 255.0
    star_render_expected_width: int = 600
    content_fov_deg: float = 110.0
    window_geometry_arg: Optional[str | tuple[int, int, int, int]] = None
    window_frame_mode: str = "frameless"
    load_last_window_geometry: (
        Callable[[], Optional[tuple[int, int, int, int]]] | None
    ) = None
    save_last_window_geometry: Callable[[int, int, int, int], None] | None = None



def _normalize_cloud_stripe_mode(mode: str) -> str:
    """Normalize cloud stripe mode to one of 'width', 'alpha', or 'halftone'."""
    mode = mode.strip().lower()
    if mode in {"alpha", "halftone"}:
        return mode
    return "width"

def _apply_visibility_boost(
    opacity: float, visibility_boost: float, tier_scale: float
) -> float:
    """Increase opacity by a tiered boost while keeping the value in range."""
    base_opacity = min(1.0, max(0.0, float(opacity)))
    boost = max(1.0, float(visibility_boost))
    if base_opacity <= 0.0 or boost <= 1.0 or tier_scale <= 0.0:
        return base_opacity
    boosted = base_opacity * (1.0 + (boost - 1.0) * float(tier_scale))
    return min(1.0, boosted)


def _apply_visibility_boost_scale(
    scale: float, visibility_boost: float, tier_scale: float
) -> float:
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
    ground_elevation_m: float,
    location_height_label: str | None,
    location_height_m: float,
    height_add_m: float,
) -> ViewerData:
    """Build the viewer input consumed by SkyWindow."""
    lat, lon, tz_name = city_data
    location_height_m = max(0.0, float(location_height_m))
    height_add_m = max(0.0, float(height_add_m))
    return ViewerData(
        location=(lat, lon),
        timezone_name=tz_name,
        city_name=city_name,
        view_center=view_center,
        edge_fov_deg=float(edge_fov_deg),
        content_fov_deg=float(content_fov_deg),
        observer_height_m=location_height_m + height_add_m,
        height_add_m=height_add_m,
        ground_elevation_m=max(0.0, float(ground_elevation_m)),
        location_height_label=location_height_label,
        location_height_m=location_height_m,
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
        star_catalog_lod6_indices=(star_catalog_np["vmag"] <= 6.0)
        .nonzero()[0]
        .astype("int32", copy=False),
        star_catalog_meta=prepare_star_catalog_meta(star_catalog),
        dso_catalog_np=None
        if dso_catalog is None
        else prepare_deep_sky_catalog_arrays(dso_catalog),
        named_stars_by_band=build_named_star_shortcuts(
            star_catalog, max_vmag=2.0, include_satellites=False
        ),
        # Search dialog resolves ISS and JPL bodies through their dedicated callbacks.
        named_stars_search_all=build_search_jump_targets(
            star_catalog, include_satellites=False
        ),
    )


def prepare_window_user_options(
    *,
    presentation_id: str = "scenic",
    star_data_policy: str = "scenic_view_scoped",
    sky_disc_alpha: float,
    sky_disc_style: str,
    sky_disc_altaz_rings: str,
    sky_disc_altaz_rings_hover: str,
    night_light_opacity: float = NIGHT_LIGHT_DEFAULT_OPACITY,
    cloud_disc_alpha: float,
    geo_satellite: bool = False,
    satellite_opacity: float,
    aircraft_opacity: float,
    tropical_cyclone_opacity: float = TROPICAL_CYCLONE_DEFAULT_OPACITY,
    terrain_horizon_opacity: float,
    earth_guide_opacity: float,
    urban_outline_opacity: float,
    water_overlay_opacity: float = 0.4,
    ground_tint_opacity: float,
    overlay_font_size: float = float(OVERLAY_FONT_SIZE_DEFAULT),
    enlarge_moon: bool,
    bright_bodies_mode: str,
    star_base_radius: float,
    vmag_limit: float,
    visual_preset: str,
    star_visibility_boost: float,
    visibility_boost: float,
    light_background_star_outline: bool = False,
    show_dso_initial: Optional[bool],
    show_asterisms_initial: Optional[bool],
    show_guidelines_initial: Optional[bool],
    observation_info_mode: Optional[str],
    sky_disc_gui_allowed: bool,
    cloud_gui_allowed: bool,
    satellite_gui_allowed: bool,
    aircraft_gui_allowed: bool,
    tropical_cyclone_gui_allowed: bool,
    terrain_horizon_gui_allowed: bool,
    earth_guide_gui_allowed: bool,
    ridge_glow_opacity: float = RIDGE_GLOW_DEFAULT_OPACITY,
    night_light_gui_allowed: bool = True,
    urban_outline_gui_allowed: bool = True,
) -> SkyWindowUserOptions:
    """Normalize user-facing options before constructing SkyWindow."""
    visibility_boost = max(1.0, float(visibility_boost))
    overlay_font_size = min(
        float(OVERLAY_FONT_SIZE_MAX),
        max(float(OVERLAY_FONT_SIZE_MIN), float(overlay_font_size)),
    )
    return SkyWindowUserOptions(
        presentation_id=str(presentation_id).strip().lower() or "scenic",
        star_data_policy=str(star_data_policy).strip().lower() or "scenic_view_scoped",
        sky_disc_alpha=_apply_visibility_boost(sky_disc_alpha, visibility_boost, 1.0),
        sky_disc_style=str(sky_disc_style).strip().lower(),
        sky_disc_altaz_rings=str(sky_disc_altaz_rings).strip().lower(),
        sky_disc_altaz_rings_hover=str(sky_disc_altaz_rings_hover).strip().lower(),
        night_light_opacity=_apply_visibility_boost(
            night_light_opacity, visibility_boost, 1.0
        ),
        ridge_glow_opacity=_apply_visibility_boost(
            ridge_glow_opacity, visibility_boost, 1.0
        ),
        cloud_disc_alpha=_apply_visibility_boost(
            cloud_disc_alpha, visibility_boost, 1.0
        ),
        geo_satellite=bool(geo_satellite),
        satellite_opacity=_apply_visibility_boost(
            satellite_opacity, visibility_boost, 1.0
        ),
        aircraft_opacity=_apply_visibility_boost(
            aircraft_opacity, visibility_boost, 1.0
        ),
        tropical_cyclone_opacity=_apply_visibility_boost(
            tropical_cyclone_opacity,
            visibility_boost,
            1.0,
        ),
        terrain_horizon_opacity=_apply_visibility_boost(
            terrain_horizon_opacity, visibility_boost, 1.0
        ),
        earth_guide_opacity=_apply_visibility_boost(
            earth_guide_opacity, visibility_boost, 1.0
        ),
        urban_outline_opacity=_apply_visibility_boost(
            urban_outline_opacity, visibility_boost, 1.0
        ),
        water_overlay_opacity=_apply_visibility_boost(
            water_overlay_opacity, visibility_boost, 1.0
        ),
        ground_tint_opacity=_apply_visibility_boost(
            ground_tint_opacity, visibility_boost, 1.0
        ),
        overlay_font_size=overlay_font_size,
        enlarge_moon=bool(enlarge_moon),
        bright_bodies_mode=str(bright_bodies_mode).strip().lower(),
        star_base_radius=max(2.0, star_base_radius),
        vmag_limit=vmag_limit,
        visual_preset=visual_preset,
        star_visibility_boost=star_visibility_boost,
        light_background_star_outline=bool(light_background_star_outline),
        asterism_visibility_boost=_apply_visibility_boost_scale(
            1.0, visibility_boost, 1.0
        ),
        earth_guide_visibility_boost=visibility_boost,
        show_dso_initial=show_dso_initial,
        show_asterisms_initial=show_asterisms_initial,
        show_guidelines_initial=show_guidelines_initial,
        observation_info_mode=observation_info_mode,
        sky_disc_gui_allowed=bool(sky_disc_gui_allowed),
        cloud_gui_allowed=bool(cloud_gui_allowed),
        satellite_gui_allowed=bool(satellite_gui_allowed),
        aircraft_gui_allowed=bool(aircraft_gui_allowed),
        tropical_cyclone_gui_allowed=bool(tropical_cyclone_gui_allowed),
        terrain_horizon_gui_allowed=bool(terrain_horizon_gui_allowed),
        earth_guide_gui_allowed=bool(earth_guide_gui_allowed),
        night_light_gui_allowed=bool(night_light_gui_allowed),
        urban_outline_gui_allowed=bool(urban_outline_gui_allowed),
    )


def prepare_window_runtime_options(
    *,
    delta_t: timedelta,
    sky_update_interval: int,
    urban_outline_radius_km: float,
    urban_outline_skyscraper_radius_km: float,
    urban_outline_min_height_m: float,
    urban_outline_max_candidates: int,
    urban_outline_feature_type: str,
    urban_outline_skyscraper_only: bool,
    urban_outline_download_timeout_seconds: float = DEFAULT_DOWNLOAD_TIMEOUT_SECONDS,
    cloud_stripe_style: tuple[int, float],
    cloud_stripe_mode: str,
    cloud_missing_tint_opacity: float,
    visibility_boost: float,
    star_render_expected_width: int,
    content_fov_deg: float,
    window_geometry_arg: Optional[str | tuple[int, int, int, int]],
    window_frame_mode: str,
    load_last_window_geometry: Callable[[], Optional[tuple[int, int, int, int]]]
    | None = None,
    save_last_window_geometry: Callable[[int, int, int, int], None] | None = None,
) -> SkyWindowRuntimeOptions:
    """Normalize runtime options before constructing SkyWindow."""
    visibility_boost = max(1.0, float(visibility_boost))
    return SkyWindowRuntimeOptions(
        delta_t=delta_t,
        sky_update_interval=max(1, int(sky_update_interval)),
        urban_outline_radius_km=max(0.0, float(urban_outline_radius_km)),
        urban_outline_skyscraper_radius_km=max(
            0.0, float(urban_outline_skyscraper_radius_km)
        ),
        urban_outline_min_height_m=max(0.0, float(urban_outline_min_height_m)),
        urban_outline_max_candidates=max(0, int(urban_outline_max_candidates)),
        urban_outline_feature_type=str(urban_outline_feature_type),
        urban_outline_skyscraper_only=bool(urban_outline_skyscraper_only),
        urban_outline_download_timeout_seconds=max(0.0, float(urban_outline_download_timeout_seconds)),
        cloud_stripe_style=cloud_stripe_style,
        cloud_stripe_mode=_normalize_cloud_stripe_mode(str(cloud_stripe_mode)),
        cloud_missing_tint_opacity=_apply_visibility_boost(
            cloud_missing_tint_opacity, visibility_boost, 1.0
        ),
        star_render_expected_width=max(1, int(star_render_expected_width)),
        content_fov_deg=max(90.0, min(127.0, float(content_fov_deg))),
        window_geometry_arg=window_geometry_arg,
        window_frame_mode=(
            "window"
            if str(window_frame_mode).strip().lower() == "window"
            else "frameless"
        ),
        load_last_window_geometry=load_last_window_geometry,
        save_last_window_geometry=save_last_window_geometry,
    )
