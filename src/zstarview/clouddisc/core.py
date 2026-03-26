# -*- coding: utf-8 -*-
"""
Core module for CloudDisc rendering.

This module contains the main `CloudDisc` class, which orchestrates the entire
process of fetching satellite data, processing it, and rendering a cloud image
from a specific observer's perspective.
"""

import datetime as dt
import logging
from typing import Optional, Tuple
from dataclasses import replace

import numpy as np
from PIL import Image

from .config import CloudDiscConfig
from .projectors.az import az_project_lonlat_grid
from .providers.goes import GoesProvider
from .providers.hima import HimaProvider
from .providers.select import GOES_SATELLITES, SUPPORTED_SATELLITES, pick_satellite, visible_satellites
from .render.grayscale import convert_bt_to_la_image
from .sampling.bt_sampler import build_bt_sampler
from .sampling.estimate_bt_warm_cold import estimate_bt_warm_from_equator_band, estimate_bt_cold_hybrid
from .types import CloudMeta, CloudSourceData, RenderKey, SourceKey, VisibilityError, round_down_utc_to_slot


logger = logging.getLogger(__name__)


class CloudDisc:
    """
    A class to render cloud images from satellite data.

    This class brings together all the components of the clouddisc library:
    - Data providers (GOES, Himawari)
    - Projection logic
    - Sampling and interpolation
    - Image rendering
    """

    def __init__(self, cfg: CloudDiscConfig, **kwargs):
        """
        Initializes the CloudDisc renderer.

        Args:
            cfg: An instance of CloudDiscConfig containing the configuration.
            **kwargs: Additional configuration options to override `cfg`.
        """
        if kwargs:
            cfg = replace(cfg, **kwargs)

        self.cfg: CloudDiscConfig = cfg
        self.goes: GoesProvider = GoesProvider(cfg)
        self.hima: HimaProvider = HimaProvider(cfg)

    def _now_rounded(self) -> dt.datetime:
        """
        Gets the current UTC time rounded down to the nearest 10-minute interval.
        Satellite data is typically published in 10-minute cycles.

        Returns:
            The rounded datetime object.
        """
        t = dt.datetime.now(dt.timezone.utc).replace(second=0, microsecond=0)
        return t.replace(minute=(t.minute // 10) * 10)

    def _select_satellite(self, lat: float, lon: float) -> str:
        supported_visible = tuple(visible_satellites(lat, lon, SUPPORTED_SATELLITES))
        if not supported_visible:
            raise VisibilityError("No supported satellite for this region")
        sat = pick_satellite(
            lat,
            lon,
            priority=self.cfg.sat_priority,
        )
        logger.debug("Selected satellite=%s for observer at (lat=%.2f, lon=%.2f)", sat, lat, lon)
        return sat

    def make_source_key(
        self,
        *,
        lat: float,
        lon: float,
        when_utc: Optional[dt.datetime] = None,
    ) -> SourceKey:
        """Build a source key for cloud data fetch/cache lookup."""
        when = self._now_rounded() if when_utc is None else round_down_utc_to_slot(when_utc)
        sat = self._select_satellite(lat, lon)
        provider = "GOES" if sat in GOES_SATELLITES else "HIMAWARI"
        return SourceKey(
            satellite=sat,
            provider=provider,
            timeslot_utc=when,
            sat_priority=self.cfg.sat_priority,
        )

    def fetch_source(
        self,
        *,
        lat: float,
        lon: float,
        when_utc: Optional[dt.datetime] = None,
    ) -> CloudSourceData:
        """Fetch cloud source data independently from camera-dependent rendering."""
        source_key = self.make_source_key(lat=lat, lon=lon, when_utc=when_utc)
        sat = source_key.satellite
        when = source_key.timeslot_utc
        sat_used = sat
        if sat in GOES_SATELLITES:
            goes_visible = tuple(visible_satellites(lat, lon, GOES_SATELLITES))
            res, sat_used = self.goes.fetch_bt_c13_with_failover(
                sat=sat,
                when_utc=when,
                allowed_sats=goes_visible,
            )
            da, used_time, src_paths = res
            product = "CMIPF-C13"
        elif sat == "HIMAWARI":
            da, used_time, src_paths = self.hima.fetch_bt_c13(
                when_utc=when,
                observer_lat=lat,
                observer_lon=lon,
                cloud_shell_km=6371.0 + 5.0,
            )
            product = "ISatSS-B13"
        else:
            raise VisibilityError(f"No suitable satellite provider found for '{sat}'")
        logger.info("Using %s (%s) data from time=%s", sat_used, product, used_time.isoformat())
        return CloudSourceData(
            source_key=SourceKey(
                satellite=sat_used,
                provider=("GOES" if sat_used in GOES_SATELLITES else "HIMAWARI"),
                timeslot_utc=source_key.timeslot_utc,
                sat_priority=source_key.sat_priority,
            ),
            data_array=da,
            satellite=sat_used,
            product=product,
            time_utc=used_time,
            src_paths=src_paths,
        )

    def render_from_source(
        self,
        *,
        source: CloudSourceData,
        lat: float,
        lon: float,
        alt: float,
        az: float,
        radius_px: int,
        edge_fov_deg: float = 90.0,
        mask_fov_deg: float = 90.0,
        cloud_shell_km: float = 6371.0 + 5.0,
    ) -> Tuple[Image.Image, CloudMeta]:
        """Render a cloud image from pre-fetched source data."""
        img, meta, _missing_mask, _coverage_ratio = self.render_from_source_with_coverage(
            source=source,
            lat=lat,
            lon=lon,
            alt=alt,
            az=az,
            radius_px=radius_px,
            edge_fov_deg=edge_fov_deg,
            mask_fov_deg=mask_fov_deg,
            cloud_shell_km=cloud_shell_km,
        )
        return img, meta

    def render_from_source_with_coverage(
        self,
        *,
        source: CloudSourceData,
        lat: float,
        lon: float,
        alt: float,
        az: float,
        radius_px: int,
        edge_fov_deg: float = 90.0,
        mask_fov_deg: float = 90.0,
        cloud_shell_km: float = 6371.0 + 5.0,
    ) -> Tuple[Image.Image, CloudMeta, np.ndarray, float]:
        """Render from pre-fetched source and return missing-data mask/coverage."""
        render_key = RenderKey(
            source=source.source_key,
            alt_deg=alt,
            az_deg=az,
            radius_px=radius_px,
            edge_fov_deg=edge_fov_deg,
            mask_fov_deg=mask_fov_deg,
        )
        return self._render_from_source_impl(
            source=source,
            lat=lat,
            lon=lon,
            render_key=render_key,
            cloud_shell_km=cloud_shell_km,
        )

    def _render_from_source_impl(
        self,
        *,
        source: CloudSourceData,
        lat: float,
        lon: float,
        render_key: RenderKey,
        cloud_shell_km: float,
    ) -> Tuple[Image.Image, CloudMeta, np.ndarray, float]:
        # Step 3: Create a sampler function: (lon, lat) -> Brightness Temperature [K]
        sampler = build_bt_sampler(source.data_array)

        # Step 4: Project a grid of longitude/latitude points that corresponds to the pixels
        # in the final image, as seen from the observer's perspective.
        lon_grid, lat_grid, mask_inside = az_project_lonlat_grid(
            lat0_deg=lat,
            lon0_deg=lon,
            alt0_deg=render_key.alt_deg,
            az0_deg=render_key.az_deg,
            radius_px=render_key.radius_px + 1,  # Project slightly larger grid for better interpolation at edges
            cloud_shell_km=cloud_shell_km,
            alt_min_deg=self.cfg.alt_min_deg,
            edge_fov_deg=render_key.edge_fov_deg,
            mask_fov_deg=render_key.mask_fov_deg,
        )

        # Step 5: Sample the brightness temperatures for each point in the projected grid.
        bt = sampler(lon_grid, lat_grid)
        finite_bt = np.isfinite(bt)
        inside_count = int(np.count_nonzero(mask_inside))
        if inside_count > 0:
            valid_inside = int(np.count_nonzero(mask_inside & finite_bt))
            coverage_ratio = float(valid_inside) / float(inside_count)
        else:
            coverage_ratio = 1.0
        missing_mask = mask_inside & ~finite_bt

        # Step 6: Estimate the warm and cold brightness temperature thresholds. These values
        # are used to map the temperature range to the grayscale color range, enhancing contrast.
        bt_warm, sample_arr = estimate_bt_warm_from_equator_band(
            source.data_array,
            lon_center_deg=lon,
            delta_lon=60.0,
            equator_lat=0.0,
            warm_p=97.0,
            half=5,
        )
        bt_cold = estimate_bt_cold_hybrid(bt, mask_inside, sample_arr, bt_warm, cold_local_p=5.0, cold_eq_p=3.0)
        inside_vals = bt[mask_inside & finite_bt].astype(np.float64)
        if inside_vals.size > 0:
            p05, p25, p50, p75, p95 = np.percentile(inside_vals, [5, 25, 50, 75, 95])
            logger.debug(
                (
                    "Cloud BT stats: sat=%s product=%s coverage=%.1f%% "
                    "warm=%.2f cold=%.2f p05=%.2f p25=%.2f p50=%.2f p75=%.2f p95=%.2f "
                    "inside_valid=%d eq_samples=%d"
                ),
                source.satellite,
                source.product,
                coverage_ratio * 100.0,
                bt_warm,
                bt_cold,
                p05,
                p25,
                p50,
                p75,
                p95,
                inside_vals.size,
                sample_arr.size,
            )
        else:
            logger.debug(
                "Cloud BT stats: sat=%s product=%s coverage=%.1f%% warm=%.2f cold=%.2f inside_valid=0 eq_samples=%d",
                source.satellite,
                source.product,
                coverage_ratio * 100.0,
                bt_warm,
                bt_cold,
                sample_arr.size,
            )

        # Step 7: Render the final grayscale (Luminance-Alpha) image from the BT data.
        img = convert_bt_to_la_image(bt, mask_inside, bt_warm, bt_cold)

        meta = CloudMeta(
            satellite=source.satellite,
            product=source.product,
            time_utc=source.time_utc,
            src_paths=source.src_paths,
        )
        return img, meta, missing_mask.astype(np.uint8), coverage_ratio
