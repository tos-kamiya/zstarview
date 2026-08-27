#!/usr/bin/env python3
"""Export a B13-only and B13+B16 cloud-layer comparison image."""

from __future__ import annotations

import argparse
import datetime as dt
from dataclasses import replace
from pathlib import Path

from PIL import Image

from zstarview.clouddisc import CloudDisc, CloudDiscConfig
from zstarview.clouddisc.types import CloudBandData, CloudSourceData, SourceKey
from zstarview.gui.cloud_render import _render_halftone_cloud_rgba_from_altaz_grid
from zstarview.paths import CACHE_PATH, CLOUD_HATCH_DEFAULT, CLOUD_SHELLS_KM
from zstarview.render.geometry import get_screen_geometry
from zstarview.types import ViewProjection


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lat", type=float, default=34.6551)
    parser.add_argument("--lon", type=float, default=133.9195)
    parser.add_argument("--alt", type=float, default=35.0)
    parser.add_argument("--az", type=float, default=0.0)
    parser.add_argument("--when-utc", type=str, default=None)
    parser.add_argument(
        "--local-slot",
        type=str,
        default=None,
        help="Use an existing cache slot such as 2026/08/27/1450 without network access.",
    )
    parser.add_argument("--size", type=int, default=600)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    disc = CloudDisc(
        CloudDiscConfig(
            cache_dir=CACHE_PATH,
            sat_priority=("AUTO",),
            bt_warm_k=310.0,
            bt_cold_k=190.0,
            alt_min_deg=1.0,
            search_back_minutes=120,
        )
    )
    when_utc = None
    if args.when_utc is not None:
        when_utc = dt.datetime.fromisoformat(args.when_utc.replace("Z", "+00:00"))
    if args.local_slot is not None:
        slot = Path(CACHE_PATH) / "hima_isatss" / "noaa-himawari9" / "AHI-L2-FLDK-ISatSS" / args.local_slot
        b13_paths = sorted(slot.glob("*M1C13*.nc"))
        b16_paths = sorted(slot.glob("*M1C16*.nc"))
        if not b13_paths or not b16_paths:
            raise RuntimeError(f"No cached B13/B16 tiles found in {slot}")
        b13_data = disc.hima._stitch_local_paths(
            b13_paths, source_label=str(slot), observer_lat=args.lat, observer_lon=args.lon
        )
        b16_data = disc.hima._stitch_local_paths(
            b16_paths, source_label=str(slot), observer_lat=args.lat, observer_lon=args.lon
        )
        slot_time = dt.datetime.strptime(args.local_slot, "%Y/%m/%d/%H%M").replace(tzinfo=dt.timezone.utc)
        source = CloudSourceData(
            source_key=SourceKey("HIMAWARI", slot_time, "HIMAWARI"),
            data_array=b13_data,
            satellite="HIMAWARI",
            product="ISatSS-B13",
            time_utc=slot_time,
            src_paths=b13_paths,
            auxiliary_bands={
                "B16": CloudBandData(16, b16_data, "ISatSS-B16", slot_time, b16_paths)
            },
        )
    else:
        source = disc.fetch_source(lat=args.lat, lon=args.lon, when_utc=when_utc)
    source_b13_only = replace(source, auxiliary_bands={})

    common = {
        "lat": args.lat,
        "lon": args.lon,
        "cloud_shells_km": CLOUD_SHELLS_KM,
    }
    grid_b13 = disc.build_altaz_grid_from_source(source=source_b13_only, **common)
    grid_b13_b16 = disc.build_altaz_grid_from_source(source=source, **common)

    geometry = get_screen_geometry(
        args.size, args.size, args.alt, edge_fov_deg=90.0, content_fov_deg=115.0
    )
    projection = ViewProjection(
        view_center=(args.alt, args.az), edge_fov_deg=90.0, content_fov_deg=115.0
    )
    shell_phases = ((0.000, 0.000), (0.500, 0.500), (-0.183, 0.683))

    def render_grid(grid: object) -> object:
        shell_images = []
        for shell_index, shell_amount in enumerate(grid.shell_amounts or ()):
            shell_grid = replace(grid, amount=shell_amount, shell_amounts=None)
            shell_images.append(
                _render_halftone_cloud_rgba_from_altaz_grid(
                    shell_grid,
                    args.size,
                    args.size,
                    CLOUD_HATCH_DEFAULT,
                    geometry=geometry,
                    projection=projection,
                    target_stripes=30,
                    width_factor=1.7,
                    density_reference_size=(args.size, args.size),
                    grid_phase=shell_phases[shell_index % len(shell_phases)],
                )
            )
        # Match the cloud compositing at a diagnostic level: shell opacity is
        # represented by the current default cloud alpha (0.11).
        import numpy as np

        out = np.zeros((args.size, args.size, 4), dtype=np.uint8)
        for image in reversed(shell_images):
            alpha = (image[..., 3].astype(np.float32) * 0.11).astype(np.uint8)
            inv = 255 - alpha
            out[..., :3] = (
                image[..., :3].astype(np.uint16) * alpha[..., None]
                + out[..., :3].astype(np.uint16) * inv[..., None]
            ) // 255
            out[..., 3] = alpha + (out[..., 3].astype(np.uint16) * inv) // 255
        return out

    rendered: dict[str, object] = {}
    for name, grid in (
        ("okayama_north_b13_only.png", grid_b13),
        ("okayama_north_b13_b16.png", grid_b13_b16),
    ):
        rgba = render_grid(grid)
        Image.fromarray(rgba, mode="RGBA").save(output_dir / name)
        rendered[name] = rgba

    # Also provide opaque alpha visualizations, since a transparent white layer
    # looks nearly binary when viewed over a black background.
    for name, rgba in rendered.items():
        alpha = rgba[..., 3]
        Image.fromarray(alpha, mode="L").save(
            output_dir / name.replace(".png", "_alpha.png")
        )
    b13_alpha = rendered["okayama_north_b13_only.png"][..., 3].astype("int16")
    b16_alpha = rendered["okayama_north_b13_b16.png"][..., 3].astype("int16")
    difference = ((b16_alpha - b13_alpha) * 4 + 128).clip(0, 255).astype("uint8")
    Image.fromarray(difference, mode="L").save(output_dir / "okayama_north_b13_b16_difference_x4.png")

    print(f"source_time_utc={source.time_utc.isoformat()}")
    print(f"b13_algorithm={grid_b13.algorithm_version}")
    print(f"b13_b16_algorithm={grid_b13_b16.algorithm_version}")
    print(f"output_dir={output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
