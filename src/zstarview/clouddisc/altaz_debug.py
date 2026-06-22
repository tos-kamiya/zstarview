# -*- coding: utf-8 -*-
"""Debug helpers for the alt/az cloud grid.

These utilities are intended for development and diagnostics only. They save
matplotlib visualisations of the `CloudAltAzGrid` so that the ingestion
result can be inspected without running the full GUI.
"""

from __future__ import annotations

import datetime as dt
import logging
from pathlib import Path
from typing import Tuple



try:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
except ImportError:  # pragma: no cover - debug-only dependency
    plt = None  # type: ignore[assignment]

from .altaz_grid import CloudAltAzGrid
from .altaz_constants import ALT_AZ_DEBUG_SUBDIR

logger = logging.getLogger(__name__)


def save_altaz_grid_debug_image(
    grid: CloudAltAzGrid,
    cache_root: Path,
    *,
    view_center: Tuple[float, float] | None = None,
    figure_size_inches: Tuple[float, float] = (12.0, 6.0),
) -> Path | None:
    """Save a debug visualisation of a `CloudAltAzGrid` to disk.

    The image contains two panels:

    1. A polar-ish heatmap of the full grid (azimuth on x, altitude on y).
    2. A small text block with source metadata and coverage statistics.

    Args:
        grid: the grid to visualise.
        cache_root: top-level cache directory; the image is written under
            ``cache_root / ALT_AZ_DEBUG_SUBDIR``.
        view_center: optional (alt_deg, az_deg) to mark on the heatmap.
        figure_size_inches: matplotlib figure size.

    Returns:
        Path to the written PNG, or ``None`` if matplotlib is unavailable.
    """
    if plt is None:
        logger.warning("matplotlib is not installed; skipping alt/az debug image")
        return None

    out_dir = cache_root / ALT_AZ_DEBUG_SUBDIR
    out_dir.mkdir(parents=True, exist_ok=True)

    ts = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%d_%H%M%S")
    filename = f"altaz_grid_{grid.satellite}_{ts}.png"
    out_path = out_dir / filename

    fig, axes = plt.subplots(1, 2, figsize=figure_size_inches, width_ratios=[3, 1])
    ax_img = axes[0]
    ax_text = axes[1]

    extent = [
        grid.az_min_deg,
        grid.az_max_deg,
        grid.alt_min_deg,
        grid.alt_max_deg,
    ]
    im = ax_img.imshow(
        grid.amount,
        origin="lower",
        aspect="auto",
        extent=extent,
        cmap="Greys",
        vmin=0.0,
        vmax=1.0,
    )
    ax_img.set_xlabel("Azimuth (deg)")
    ax_img.set_ylabel("Altitude (deg)")
    ax_img.set_title(f"CloudAltAzGrid — {grid.satellite} @ {grid.time_utc:%H:%M} UTC")
    fig.colorbar(im, ax=ax_img, label="Cloud amount")

    if view_center is not None:
        alt_c, az_c = view_center
        ax_img.axvline(az_c, color="red", linestyle="--", alpha=0.5)
        ax_img.axhline(alt_c, color="red", linestyle="--", alpha=0.5)

    valid = grid.amount[grid.missing_mask == 0]
    mean_amount = float(valid.mean()) if valid.size else 0.0
    coverage = float(grid.coverage_ratio)
    completeness = grid.source_completeness_ratio
    meta_lines = [
        f"Satellite: {grid.satellite}",
        f"Product: {grid.product}",
        f"Time: {grid.time_utc:%Y-%m-%d %H:%M:%S} UTC",
        f"Observer: ({grid.observer_lat:.4f}, {grid.observer_lon:.4f})",
        f"Shells: {grid.shells_km}",
        f"Grid: {grid.amount.shape[1]} x {grid.amount.shape[0]} "
        f"({grid.grid_resolution_deg:g} deg)",
        f"Coverage: {coverage * 100.0:.1f}%",
        f"Source completeness: {completeness * 100.0:.1f}%"
        if completeness is not None
        else "Source completeness: N/A",
        f"Mean amount: {mean_amount:.3f}",
    ]
    if view_center is not None:
        meta_lines.append(f"View center: ({view_center[0]:.2f}, {view_center[1]:.2f})")

    ax_text.axis("off")
    ax_text.text(
        0.05,
        0.95,
        "\n".join(meta_lines),
        transform=ax_text.transAxes,
        fontsize=9,
        verticalalignment="top",
        fontfamily="monospace",
        wrap=True,
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    logger.info("Saved alt/az grid debug image to %s", out_path)
    return out_path
