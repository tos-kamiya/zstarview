# -*- coding: utf-8 -*-
"""
Core module for CloudDisc rendering.

This module contains the main `CloudDisc` class, which orchestrates the entire
process of fetching satellite data, processing it, and rendering a cloud image
from a specific observer's perspective.
"""

import datetime as dt
import logging
import threading
from dataclasses import replace
from typing import Optional, Sequence, Tuple

import numpy as np

from .config import CloudDiscConfig
from .projectors.az import az_project_lonlat_grid
from .providers.goes import GoesProvider
from .providers.hima import HimaProvider
from .providers.select import (
    GOES_SATELLITES,
    SUPPORTED_SATELLITES,
    pick_satellite,
    visible_satellites,
)
from .render.grayscale import (
    _bt_to_weight,
    _suppress_low_cloud_weight,
    convert_bt_to_rgba_image,
)
from .sampling.bt_sampler import build_bt_sampler
from .sampling.estimate_bt_warm_cold import (
    estimate_bt_cold_hybrid,
    estimate_bt_warm_from_equator_band,
)
from .types import (
    CloudMeta,
    CloudSourceData,
    RenderKey,
    SourceKey,
    VisibilityError,
    round_down_utc_to_slot,
)

logger = logging.getLogger(__name__)

DEFAULT_CLOUD_SHELLS_KM: tuple[float, ...] = (6371.0 + 3.0, 6371.0 + 5.0, 6371.0 + 7.0)
DEFAULT_CLOUD_SHELL_WEIGHTS: tuple[float, ...] = (0.20, 0.60, 0.20)
DEFAULT_CLOUD_SHELL_LOW_WEIGHTS: tuple[float, ...] = (0.0, 1.0, 0.0)
DEFAULT_CLOUD_SHELL_BLEND_LOW_CLOUD_AMOUNT: float = 0.25
DEFAULT_CLOUD_SHELL_BLEND_HIGH_CLOUD_AMOUNT: float = 0.65


def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    """Return a smooth 0..1 ramp between the given thresholds."""
    denom = max(1e-6, float(edge1) - float(edge0))
    t = float(np.clip((float(x) - float(edge0)) / denom, 0.0, 1.0))
    return t * t * (3.0 - 2.0 * t)


def _estimate_scene_cloud_amount(
    bt: np.ndarray,
    mask_inside: np.ndarray,
    bt_warm: float,
    bt_cold: float,
) -> float:
    """Estimate scene-wide cloudiness from a representative shell."""
    inside = mask_inside & np.isfinite(bt)
    if not np.any(inside):
        return 0.0
    weight = _bt_to_weight(bt, bt_warm, bt_cold)
    weight = _suppress_low_cloud_weight(weight)
    return float(np.mean(weight[inside]))


def _blend_cloud_shell_weights(
    cloud_amount: float,
    shell_radii_km: Sequence[float],
) -> tuple[float, ...]:
    """Blend shell weights from middle-only to the default 3-layer mix."""
    if len(shell_radii_km) != 3:
        return tuple(1.0 / float(len(shell_radii_km)) for _ in shell_radii_km)
    t = _smoothstep(
        DEFAULT_CLOUD_SHELL_BLEND_LOW_CLOUD_AMOUNT,
        DEFAULT_CLOUD_SHELL_BLEND_HIGH_CLOUD_AMOUNT,
        cloud_amount,
    )
    low = np.asarray(DEFAULT_CLOUD_SHELL_LOW_WEIGHTS, dtype=np.float64)
    high = np.asarray(DEFAULT_CLOUD_SHELL_WEIGHTS, dtype=np.float64)
    blended = low * (1.0 - t) + high * t
    return tuple(float(v) for v in blended)


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
        self._last_bt_warm: float = float(cfg.bt_warm_k)

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
        cloud_shells_km: Sequence[float] = DEFAULT_CLOUD_SHELLS_KM,
        abort_event: threading.Event | None = None,
    ) -> CloudSourceData:
        """Fetch cloud source data independently from camera-dependent rendering."""
        source_key = self.make_source_key(lat=lat, lon=lon, when_utc=when_utc)
        sat = source_key.satellite
        when = source_key.timeslot_utc
        sat_used = sat
        shell_max_km = max(float(v) for v in cloud_shells_km) if cloud_shells_km else (6371.0 + 5.0)
        if sat in GOES_SATELLITES:
            goes_visible = tuple(visible_satellites(lat, lon, GOES_SATELLITES))
            res, sat_used = self.goes.fetch_bt_c13_with_failover(
                sat=sat,
                when_utc=when,
                allowed_sats=goes_visible,
                abort_event=abort_event,
            )
            da, used_time, src_paths = res
            product = "CMIPF-C13"
        elif sat == "HIMAWARI":
            da, used_time, src_paths = self.hima.fetch_bt_c13(
                when_utc=when,
                observer_lat=lat,
                observer_lon=lon,
                cloud_shell_km=shell_max_km,
                abort_event=abort_event,
            )
            product = "ISatSS-B13"
        else:
            raise VisibilityError(f"No suitable satellite provider found for '{sat}'")
        logger.info("Using %s (%s) data from time=%s", sat_used, product, used_time.isoformat())
        source_expected_count = getattr(da, "attrs", {}).get("source_expected_count")
        source_available_count = getattr(da, "attrs", {}).get("source_available_count")
        source_completeness_ratio = getattr(da, "attrs", {}).get("source_completeness_ratio")
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
            source_expected_count=int(source_expected_count) if source_expected_count is not None else None,
            source_available_count=int(source_available_count) if source_available_count is not None else None,
            source_completeness_ratio=(
                float(source_completeness_ratio) if source_completeness_ratio is not None else None
            ),
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
        cloud_shells_km: Sequence[float] = DEFAULT_CLOUD_SHELLS_KM,
    ) -> Tuple[np.ndarray, CloudMeta]:
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
            cloud_shells_km=cloud_shells_km,
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
        cloud_shells_km: Sequence[float] = DEFAULT_CLOUD_SHELLS_KM,
    ) -> Tuple[np.ndarray, CloudMeta, np.ndarray, float]:
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
            cloud_shells_km=cloud_shells_km,
        )

    def _render_from_source_impl(
        self,
        *,
        source: CloudSourceData,
        lat: float,
        lon: float,
        render_key: RenderKey,
        cloud_shells_km: Sequence[float],
    ) -> Tuple[np.ndarray, CloudMeta, np.ndarray, float]:
        shell_radii_km = tuple(float(v) for v in cloud_shells_km) if cloud_shells_km else DEFAULT_CLOUD_SHELLS_KM
        # Step 3: Create a sampler function: (lon, lat) -> Brightness Temperature [K]
        sampler = source.sampler
        if sampler is None:
            sampler = build_bt_sampler(source.data_array)
            source.sampler = sampler

        projected_layers: list[tuple[np.ndarray, np.ndarray]] = []
        combined_inside_mask: np.ndarray | None = None
        combined_valid_mask: np.ndarray | None = None
        bt_for_threshold: np.ndarray | None = None
        threshold_mask_inside: np.ndarray | None = None

        for cloud_shell_km in shell_radii_km:
            lon_grid, lat_grid, mask_inside = az_project_lonlat_grid(
                lat0_deg=lat,
                lon0_deg=lon,
                alt0_deg=render_key.alt_deg,
                az0_deg=render_key.az_deg,
                radius_px=render_key.radius_px + 1,
                cloud_shell_km=cloud_shell_km,
                alt_min_deg=self.cfg.alt_min_deg,
                edge_fov_deg=render_key.edge_fov_deg,
                mask_fov_deg=render_key.mask_fov_deg,
            )
            bt = sampler(lon_grid, lat_grid)
            finite_bt = np.isfinite(bt)
            projected_layers.append((bt, mask_inside))
            combined_inside_mask = mask_inside if combined_inside_mask is None else (combined_inside_mask | mask_inside)
            valid_inside = mask_inside & finite_bt
            combined_valid_mask = valid_inside if combined_valid_mask is None else (combined_valid_mask | valid_inside)
            if bt_for_threshold is None:
                bt_for_threshold = bt
                threshold_mask_inside = mask_inside

        assert bt_for_threshold is not None
        assert threshold_mask_inside is not None
        assert combined_inside_mask is not None
        assert combined_valid_mask is not None

        inside_count = int(np.count_nonzero(combined_inside_mask))
        if inside_count > 0:
            valid_inside = int(np.count_nonzero(combined_valid_mask))
            coverage_ratio = float(valid_inside) / float(inside_count)
        else:
            coverage_ratio = 1.0
        missing_mask = combined_inside_mask & ~combined_valid_mask

        equator_band_missing = bool(getattr(source.data_array, "attrs", {}).get("equator_band_missing", False))
        if equator_band_missing:
            bt_warm = float(self._last_bt_warm)
            sample_arr = np.array([], dtype=np.float32)
            logger.info("Reusing previous Himawari bt_warm=%.2f because equator-band tiles are missing", bt_warm)
        else:
            bt_warm, sample_arr = estimate_bt_warm_from_equator_band(
                source.data_array,
                lon_center_deg=lon,
                delta_lon=60.0,
                equator_lat=0.0,
                warm_p=97.0,
                half=5,
            )
            if not np.isfinite(bt_warm):
                bt_warm = float(self._last_bt_warm)
                logger.info("Reusing previous Himawari bt_warm=%.2f because equator-band estimate was unavailable", bt_warm)
        bt_cold = estimate_bt_cold_hybrid(
            bt_for_threshold,
            threshold_mask_inside,
            sample_arr,
            bt_warm,
            cold_local_p=5.0,
            cold_eq_p=3.0,
        )
        cloud_amount = _estimate_scene_cloud_amount(bt_for_threshold, threshold_mask_inside, bt_warm, bt_cold)
        blend_weights = _blend_cloud_shell_weights(cloud_amount, shell_radii_km)
        inside_vals = bt_for_threshold[threshold_mask_inside & np.isfinite(bt_for_threshold)].astype(np.float64)
        if inside_vals.size > 0:
            p05, p25, p50, p75, p95 = np.percentile(inside_vals, [5, 25, 50, 75, 95])
            logger.debug(
                (
                    "Cloud BT stats: sat=%s product=%s coverage=%.1f%% cloud_amount=%.3f "
                    "warm=%.2f cold=%.2f p05=%.2f p25=%.2f p50=%.2f p75=%.2f p95=%.2f "
                    "inside_valid=%d eq_samples=%d"
                ),
                source.satellite,
                source.product,
                coverage_ratio * 100.0,
                cloud_amount,
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
                "Cloud BT stats: sat=%s product=%s coverage=%.1f%% cloud_amount=%.3f warm=%.2f cold=%.2f inside_valid=0 eq_samples=%d",
                source.satellite,
                source.product,
                coverage_ratio * 100.0,
                cloud_amount,
                bt_warm,
                bt_cold,
                sample_arr.size,
            )

        # Step 7: Render each shell separately, then blend them with fixed weights.
        img_acc: np.ndarray | None = None
        for (bt, mask_inside), blend_weight in zip(projected_layers, blend_weights, strict=True):
            layer_img = convert_bt_to_rgba_image(bt, mask_inside, bt_warm, bt_cold).astype(np.float32)
            if img_acc is None:
                img_acc = layer_img * blend_weight
            else:
                img_acc += layer_img * blend_weight
        assert img_acc is not None
        img = np.clip(np.round(img_acc), 0, 255).astype(np.uint8)
        self._last_bt_warm = float(bt_warm)

        meta = CloudMeta(
            satellite=source.satellite,
            product=source.product,
            time_utc=source.time_utc,
            src_paths=source.src_paths,
        )
        return img, meta, missing_mask.astype(np.uint8), coverage_ratio
