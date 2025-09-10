# -*- coding: utf-8 -*-
"""
Core module for CloudDisc rendering.

This module contains the main `CloudDisc` class, which orchestrates the entire
process of fetching satellite data, processing it, and rendering a cloud image
from a specific observer's perspective.
"""

import datetime as dt
import logging
from typing import Tuple
from dataclasses import replace

from PIL import Image

from .config import CloudDiscConfig
from .projectors.az import az_project_lonlat_grid
from .providers.goes import GoesProvider
from .providers.hima import HimaProvider
from .providers.select import pick_satellite
from .render.grayscale import convert_bt_to_la_image
from .sampling.bt_sampler import build_bt_sampler
from .sampling.estimate_bt_warm_cold import estimate_bt_warm_from_equator_band, estimate_bt_cold_hybrid
from .types import CloudMeta, VisibilityError


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

    def render_now(
        self,
        lat: float,
        lon: float,
        alt: float,
        az: float,
        radius_px: int,
        edge_fov_deg: float = 90.0,
        mask_fov_deg: float = 90.0,
        cloud_shell_km: float = 6371.0 + 5.0,  # 5km above Earth's surface
    ) -> tuple[Image.Image, CloudMeta]:
        """
        Renders a cloud image for the current time from the observer's perspective.

        The method performs the following steps:
        1. Selects the most suitable satellite based on the observer's location.
        2. Fetches the brightness temperature (BT) data from the selected satellite provider.
        3. Creates a sampler to interpolate BT values for any given longitude and latitude.
        4. Projects a lon/lat grid representing the observer's view.
        5. Samples the BT values onto this grid.
        6. Estimates warm/cold temperature thresholds to map temperatures to pixel values.
        7. Renders the final grayscale (Luminance-Alpha) image.

        Args:
            lat: Observer's latitude in degrees.
            lon: Observer's longitude in degrees.
            alt: Observer's viewing altitude in degrees (90 is zenith).
            az: Observer's viewing azimuth in degrees (0 is North).
            radius_px: The radius of the output circular image in pixels.
            edge_fov_deg: The field of view angle to the edge of the rendered disc.
            mask_fov_deg: The field of view angle for the visibility mask.

        Returns:
            A tuple containing the rendered LA (Luminance-Alpha) PIL Image and a CloudMeta object.
        """
        when = self._now_rounded()

        # Clamp altitude to avoid math errors near the poles (altitudes of +/- 90).
        eps = 1e-3
        alt = max(-(90.0 - eps), min(90.0 - eps, alt))

        # Step 1: Automatically select the best satellite based on location or config priority.
        sat = pick_satellite(lat, lon, priority=self.cfg.sat_priority)
        logger.info("Selected satellite=%s for observer at (lat=%.2f, lon=%.2f)", sat, lat, lon)

        # Step 2: Fetch data from the provider. This returns a DataArray of brightness temperatures.
        sat_used = sat
        if sat in ("G16", "G18"):
            res, sat_used = self.goes.fetch_bt_c13_with_failover(sat=sat, when_utc=when)
            da, used_time, src_paths = res
            product = "CMIPF-C13"
        elif sat == "HIMAWARI":
            da, used_time, src_paths = self.hima.fetch_bt_c13(when_utc=when)
            product = "HSD-B13" if len(src_paths) > 1 else "ISatSS-B13"
        else:
            raise VisibilityError(f"No suitable satellite provider found for '{sat}'")
        logger.info("Using %s (%s) data from time=%s", sat_used, product, used_time.isoformat())

        # Step 3: Create a sampler function: (lon, lat) -> Brightness Temperature [K]
        sampler = build_bt_sampler(da)

        # Step 4: Project a grid of longitude/latitude points that corresponds to the pixels
        # in the final image, as seen from the observer's perspective.
        lon_grid, lat_grid, mask_inside = az_project_lonlat_grid(
            lat0_deg=lat,
            lon0_deg=lon,
            alt0_deg=alt,
            az0_deg=az,
            radius_px=radius_px + 1,  # Project slightly larger grid for better interpolation at edges
            cloud_shell_km=cloud_shell_km,
            alt_min_deg=self.cfg.alt_min_deg,
            edge_fov_deg=edge_fov_deg,
            mask_fov_deg=mask_fov_deg,
        )

        # Step 5: Sample the brightness temperatures for each point in the projected grid.
        bt = sampler(lon_grid, lat_grid)

        # Step 6: Estimate the warm and cold brightness temperature thresholds. These values
        # are used to map the temperature range to the grayscale color range, enhancing contrast.
        bt_warm, sample_arr = estimate_bt_warm_from_equator_band(da, lon_center_deg=lon, delta_lon=60.0, equator_lat=0.0, warm_p=97.0, half=5)
        bt_cold = estimate_bt_cold_hybrid(bt, mask_inside, sample_arr, bt_warm, cold_local_p=5.0, cold_eq_p=3.0)

        # Step 7: Render the final grayscale (Luminance-Alpha) image from the BT data.
        img = convert_bt_to_la_image(bt, mask_inside, bt_warm, bt_cold)

        meta = CloudMeta(satellite=sat_used, product=product, time_utc=used_time, src_paths=src_paths)
        return img, meta
