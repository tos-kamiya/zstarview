"""Fetch collocated B13/B16 data and write calibration diagnostics.

This exploratory script does not alter or feed the normal cloud renderer.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from zstarview.clouddisc.config import CloudDiscConfig
from zstarview.clouddisc.core import CloudDisc
from zstarview.clouddisc.providers.select import GOES_SATELLITES
from zstarview.clouddisc.render.grayscale import (
    _bt_to_weight,
    _suppress_low_cloud_weight,
)
from zstarview.clouddisc.sampling.b13_b16 import (
    collocate_b13_b16,
    diagnostic_summary,
)
from zstarview.clouddisc.sampling.bt_sampler import build_bt_sampler
from zstarview.clouddisc.sampling.estimate_bt_warm_cold import (
    estimate_bt_cold_hybrid,
    estimate_bt_warm_from_equator_band,
    estimate_bt_warm_hybrid,
)


def _parse_utc(value: str) -> dt.datetime:
    parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed.astimezone(dt.timezone.utc)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fetch same-scan B13/B16 data and write diagnostic artifacts."
    )
    parser.add_argument("--lat", required=True, type=float)
    parser.add_argument("--lon", required=True, type=float)
    parser.add_argument("--when-utc", type=_parse_utc)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--extent-deg", type=float, default=12.0)
    parser.add_argument("--step-deg", type=float, default=0.1)
    parser.add_argument("--search-back-minutes", type=int, default=120)
    parser.add_argument(
        "--sat-priority",
        nargs="+",
        default=("AUTO",),
        help="AUTO or an ordered list containing HIMAWARI, G19, and G18.",
    )
    return parser


def _cloud_amount(
    source_data,
    bt13: np.ndarray,
    lat_grid: np.ndarray,
    lon_grid: np.ndarray,
    observer_lat: float,
    observer_lon: float,
    extent_deg: float,
) -> tuple[np.ndarray, float, float]:
    _, eq_samples = estimate_bt_warm_from_equator_band(
        source_data,
        lon_center_deg=observer_lon,
        delta_lon=60.0,
        equator_lat=0.0,
        warm_p=97.0,
        half=5,
        equator_lat_half_band_deg=5.0,
    )
    central = (
        (np.abs(lat_grid - observer_lat) <= extent_deg * 0.25)
        & (np.abs(lon_grid - observer_lon) <= extent_deg * 0.25)
    )
    central_bt = bt13[central]
    central_mask = np.isfinite(central_bt)
    bt_warm = estimate_bt_warm_hybrid(
        central_bt,
        central_mask,
        eq_samples,
        fallback_bt_warm=310.0,
    )
    bt_cold = estimate_bt_cold_hybrid(
        central_bt,
        central_mask,
        eq_samples,
        bt_warm,
    )
    amount = _suppress_low_cloud_weight(_bt_to_weight(bt13, bt_warm, bt_cold))
    return amount.astype(np.float32), float(bt_warm), float(bt_cold)


def _write_figure(
    path: Path,
    bt13: np.ndarray,
    bt16: np.ndarray,
    delta: np.ndarray,
    amount13: np.ndarray,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    panels = (
        (bt13, "B13 brightness temperature", "turbo", None, None),
        (bt16, "B16 brightness temperature", "turbo", None, None),
        (delta, "BT13 - BT16", "coolwarm", -5.0, 30.0),
        (amount13, "Existing B13 cloud amount", "gray", 0.0, 1.0),
    )
    for axis, (values, title, cmap, vmin, vmax) in zip(
        axes.ravel(), panels, strict=True
    ):
        image = axis.imshow(values, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        axis.set_title(title)
        axis.set_xlabel("longitude sample")
        axis.set_ylabel("latitude sample")
        fig.colorbar(image, ax=axis, shrink=0.8)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> int:
    args = _parser().parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = CloudDiscConfig(
        cache_dir=args.cache_dir,
        sat_priority=tuple(args.sat_priority),
        search_back_minutes=args.search_back_minutes,
    )
    disc = CloudDisc(config)
    source = disc.fetch_source(lat=args.lat, lon=args.lon, when_utc=args.when_utc)

    if source.satellite in GOES_SATELLITES:
        bt16_data, bt16_time, bt16_paths = disc.goes.fetch_bt_c16_for_c13(
            source.satellite,
            source.time_utc,
            source.src_paths[0],
        )
    else:
        bt16_data, bt16_time, bt16_paths = disc.hima.fetch_bt_b16_for_b13(
            source.time_utc,
            source.src_paths,
        )
    if bt16_time != source.time_utc:
        raise RuntimeError("B13/B16 observation times do not match")

    extent = max(0.1, float(args.extent_deg))
    step = max(0.01, float(args.step_deg))
    lats = np.arange(args.lat - extent, args.lat + extent + step * 0.5, step)
    lons = np.arange(args.lon - extent, args.lon + extent + step * 0.5, step)
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    bt13 = build_bt_sampler(source.data_array)(lon_grid, lat_grid)
    bt16 = build_bt_sampler(bt16_data)(lon_grid, lat_grid)
    diagnostic = collocate_b13_b16(bt13, bt16)
    amount13, bt_warm, bt_cold = _cloud_amount(
        source.data_array,
        bt13,
        lat_grid,
        lon_grid,
        args.lat,
        args.lon,
        extent,
    )

    summary = diagnostic_summary(diagnostic)
    summary.update(
        {
            "satellite": source.satellite,
            "b13_product": source.product,
            "b16_product": "CMIPF-C16"
            if source.satellite in GOES_SATELLITES
            else "ISatSS-B16",
            "observation_time_utc": source.time_utc.isoformat(),
            "observer_lat": float(args.lat),
            "observer_lon": float(args.lon),
            "extent_deg": extent,
            "step_deg": step,
            "bt_warm_k": bt_warm,
            "bt_cold_k": bt_cold,
            "b13_paths": [str(path) for path in source.src_paths],
            "b16_paths": [str(path) for path in bt16_paths],
            "note": "Diagnostic only; no high-score calibration is applied.",
        }
    )
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    np.savez_compressed(
        output_dir / "samples.npz",
        latitude=lat_grid.astype(np.float32),
        longitude=lon_grid.astype(np.float32),
        bt13_k=diagnostic.bt13_k,
        bt16_k=diagnostic.bt16_k,
        delta_bt_k=diagnostic.delta_bt_k,
        valid_mask=diagnostic.valid_mask,
        amount13=amount13,
    )
    _write_figure(
        output_dir / "overview.png",
        diagnostic.bt13_k,
        diagnostic.bt16_k,
        diagnostic.delta_bt_k,
        amount13,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
