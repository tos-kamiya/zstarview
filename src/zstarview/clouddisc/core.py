"""
Core module for CloudDisc rendering.
"""

import datetime as dt
import logging
from typing import Tuple

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
    """

    def __init__(self, cfg: CloudDiscConfig):
        """
        Initializes the CloudDisc renderer.

        Args:
            cfg: An instance of CloudDiscConfig containing the configuration.
        """
        self.cfg: CloudDiscConfig = cfg
        self.goes: GoesProvider = GoesProvider(cfg)
        self.hima: HimaProvider = HimaProvider(cfg)

    def _now_rounded(self) -> dt.datetime:
        """
        Gets the current UTC time rounded down to the nearest 10 minutes.

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
    ) -> tuple[Image.Image, CloudMeta]:
        """
        Renders a cloud image for the current time (rounded to 10 minutes).

        The method performs the following steps:
        1. Selects the most suitable satellite based on the observer's location.
        2. Fetches the brightness temperature data from the selected satellite provider.
        3. Creates a sampler to interpolate brightness temperature values for any given longitude and latitude.
        4. Projects a lon/lat grid based on the observer's viewing angle.
        5. Samples the brightness temperature on the grid.
        6. Renders the final grayscale (Luminance-Alpha) image.

        Args:
            lat: Observer's latitude.
            lon: Observer's longitude.
            alt: Observer's altitude in degrees.
            az: Observer's azimuth in degrees.
            radius_px: The radius of the output image in pixels. The final image will be (2*R, 2*R).
            edge_fov_deg: Field of view in degrees.
            mask_fov_deg: Field-of-view angle for the visibility mask.

        Returns:
            A tuple containing the rendered LA image and metadata.
        """
        when = self._now_rounded()

        # Clamp altitude to avoid math errors at the poles
        eps = 1e-3
        alt = max(-(90.0 - eps), min(90.0 - eps, alt))

        # 1) Automatically select the best satellite based on location, or follow priority list.
        sat = pick_satellite(lat, lon, priority=self.cfg.sat_priority)
        logger.info("Selected satellite=%s", sat, extra={"lat": lat, "lon": lon})

        # 2) Fetch data from the provider (gets BT[K] DataArray and projection info)
        sat_used = sat
        if sat in ("G16", "G18"):
            res, sat_used = self.goes.fetch_bt_c13_with_failover(sat=sat, when_utc=when)
            da, used_time, src_paths = res
            product = "CMIPF-C13"
        elif sat == "HIMAWARI":
            da, used_time, src_paths = self.hima.fetch_bt_c13(when_utc=when)
            product = "HSD-B13" if len(src_paths) > 1 else "ISatSS-B13"
        else:
            raise VisibilityError("No suitable satellite")
        logger.info("Using %s (%s) time=%s", sat_used if sat in ("G16", "G18") else "HIMAWARI", product, used_time.isoformat())

        # 3) Create a sampler for (lon,lat) -> BT[K]
        sampler = build_bt_sampler(da)

        # 4) Project a lon/lat grid for the given viewing perspective
        lon_grid, lat_grid, mask_inside = az_project_lonlat_grid(
            lat0_deg=lat,
            lon0_deg=lon,
            alt0_deg=alt,
            az0_deg=az,
            radius_px=radius_px + 1,
            cloud_shell_km=6371.0 + 5.0,  # Earth radius + cloud height
            alt_min_deg=self.cfg.alt_min_deg,
            edge_fov_deg=edge_fov_deg,
            mask_fov_deg=mask_fov_deg,
        )

        # 5) Sample the brightness temperature on the grid
        bt = sampler(lon_grid, lat_grid)

        bt_warm, sample_arr = estimate_bt_warm_from_equator_band(da, lon_center_deg=lon, delta_lon=60.0, equator_lat=0.0, warm_p=97.0, half=5)

        bt_cold = estimate_bt_cold_hybrid(bt, mask_inside, sample_arr, bt_warm, cold_local_p=5.0, cold_eq_p=3.0)

        # 6) Render the grayscale (LA) image from the brightness temperature data
        img = convert_bt_to_la_image(bt, mask_inside, bt_warm, bt_cold)

        meta = CloudMeta(satellite=sat_used, product=product, time_utc=used_time, src_paths=src_paths)
        return img, meta
