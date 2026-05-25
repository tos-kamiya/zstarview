#!/usr/bin/env python3
"""Fit a geostationary pixel-to-lat/lon map from RGB control points."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image

from zstarview.utils.geostationary_image_mapping import (
    fit_pixel_lonlat_tps_from_image,
)


def _parse_query(value: str) -> tuple[float, float]:
    text = (value or "").strip()
    parts = [part.strip() for part in text.split(",")]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("Expected a comma-separated pixel coordinate pair, e.g. '123,456'.")
    try:
        return (float(parts[0]), float(parts[1]))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Pixel coordinates must be numeric.") from exc


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate a geostationary projection from a Meteosat or geosatellite API image "
            "and control colors."
        ),
    )
    parser.add_argument(
        "image",
        type=Path,
        help="Input Meteosat or geosatellite API image path.",
    )
    parser.add_argument(
        "control_map",
        type=Path,
        help="Control-point map such as latlonmap.txt, with lat/lon and RGB markers.",
    )
    parser.add_argument(
        "--dump-grid",
        type=Path,
        default=None,
        metavar="PATH",
        help="Write full-resolution lon/lat grids for the estimated Meteosat/geosatellite projection to a compressed .npz file.",
    )
    parser.add_argument(
        "--query",
        action="append",
        type=_parse_query,
        default=[],
        metavar="X,Y",
        help="Query an individual pixel coordinate in the input image. May be repeated.",
    )
    parser.add_argument(
        "--lon0",
        type=float,
        default=0.0,
        help="Reserved compatibility metadata (default: 0.0).",
    )
    parser.add_argument(
        "--height-m",
        type=float,
        default=35785831.0,
        help="Reserved compatibility metadata (default: 35785831.0).",
    )
    parser.add_argument(
        "--sweep-axis",
        choices=("x", "y"),
        default="x",
        help="Reserved compatibility metadata (default: x).",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit machine-readable summary and query results as JSON.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    image_path = args.image.expanduser()
    control_map_path = args.control_map.expanduser()
    if not image_path.exists():
        parser.error(f"Image file does not exist: {image_path}")
    if not control_map_path.exists():
        parser.error(f"Control map file does not exist: {control_map_path}")

    with Image.open(image_path) as image:
        pixel_map, located = fit_pixel_lonlat_tps_from_image(image, control_map_path)
        summary = {
            "model": "thin_plate_spline",
            "control_points": len(located),
            "x_center": float(pixel_map.x_center),
            "y_center": float(pixel_map.y_center),
            "x_scale": float(pixel_map.x_scale),
            "y_scale": float(pixel_map.y_scale),
        }
        residuals = []
        for point in located:
            lon_deg, lat_deg = pixel_map.pixel_to_lonlat(point.x_px, point.y_px)
            residuals.append(
                {
                    "rgb": "#{:02x}{:02x}{:02x}".format(*point.rgb),
                    "pixel_x": float(point.x_px),
                    "pixel_y": float(point.y_px),
                    "lat_deg": float(point.lat_deg),
                    "lon_deg": float(point.lon_deg),
                    "fit_lat_deg": float(np.asarray(lat_deg).item()),
                    "fit_lon_deg": float(np.asarray(lon_deg).item()),
                }
            )

        query_results = []
        for x_px, y_px in args.query:
            lon_q, lat_q = pixel_map.pixel_to_lonlat(x_px, y_px)
            query_results.append(
                {
                    "x_px": float(x_px),
                    "y_px": float(y_px),
                    "lon_deg": float(np.asarray(lon_q).item()),
                    "lat_deg": float(np.asarray(lat_q).item()),
                }
            )

        if args.dump_grid is not None:
            dump_path = args.dump_grid.expanduser()
            dump_path.parent.mkdir(parents=True, exist_ok=True)
            lon_deg_grid, lat_deg_grid = pixel_map.pixel_grid_to_lonlat(image.width, image.height)
            np.savez_compressed(
                dump_path,
                lon_deg=lon_deg_grid,
                lat_deg=lat_deg_grid,
                model=np.asarray(summary["model"], dtype=np.str_),
                x_center=np.asarray(summary["x_center"], dtype=np.float64),
                y_center=np.asarray(summary["y_center"], dtype=np.float64),
                x_scale=np.asarray(summary["x_scale"], dtype=np.float64),
                y_scale=np.asarray(summary["y_scale"], dtype=np.float64),
                x_samples=np.asarray(pixel_map.x_samples, dtype=np.float64),
                y_samples=np.asarray(pixel_map.y_samples, dtype=np.float64),
                lon_weights_radial=np.asarray(pixel_map.lon_weights["radial"], dtype=np.float64),
                lon_weights_affine=np.asarray(pixel_map.lon_weights["affine"], dtype=np.float64),
                lat_weights_radial=np.asarray(pixel_map.lat_weights["radial"], dtype=np.float64),
                lat_weights_affine=np.asarray(pixel_map.lat_weights["affine"], dtype=np.float64),
            )

    if args.json:
        print(
            json.dumps(
                {
                    "summary": summary,
                    "control_residuals": residuals,
                    "queries": query_results,
                    "dump_grid": str(args.dump_grid.expanduser()) if args.dump_grid is not None else None,
                },
                ensure_ascii=True,
                indent=2,
            )
        )
    else:
        print("Model summary:")
        print(f"  model: {summary['model']}")
        print(f"  control points: {summary['control_points']}")
        print(f"  x center/scale: {summary['x_center']}, {summary['x_scale']}")
        print(f"  y center/scale: {summary['y_center']}, {summary['y_scale']}")
        if residuals:
            xs = [abs(item["fit_lon_deg"] - item["lon_deg"]) for item in residuals]
            ys = [abs(item["fit_lat_deg"] - item["lat_deg"]) for item in residuals]
            print(f"  max control-point error: lon={max(xs):.6g} deg, lat={max(ys):.6g} deg")
        if query_results:
            print("Queries:")
            for item in query_results:
                print(
                    "  ({x_px:.3f}, {y_px:.3f}) -> lat={lat_deg:.6f}, lon={lon_deg:.6f}".format(
                        **item,
                    )
                )
        if args.dump_grid is not None:
            print(f"Wrote: {args.dump_grid.expanduser()}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
