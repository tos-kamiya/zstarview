#!/usr/bin/env python3
"""Fit an Equidistant Conic pixel-to-lat/lon map from RGB control points."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image

from zstarview.utils.geostationary_image_mapping import (
    fit_equidistant_conic_projection_from_image,
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
        description="Estimate an Equidistant Conic projection from a Meteosat or geosatellite API image and control colors.",
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
        help="Write full-resolution lon/lat grids for the best Equidistant Conic fit to a compressed .npz file.",
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
        fit, located = fit_equidistant_conic_projection_from_image(image, control_map_path)
        pixel_map = fit.pixel_map
        summary = {
            "model": "equidistant_conic",
            "control_points": len(located),
            "longitude_of_projection_origin": float(fit.longitude_of_projection_origin),
            "latitude_of_projection_origin": float(fit.latitude_of_projection_origin),
            "standard_parallel_1": float(fit.standard_parallel_1),
            "standard_parallel_2": float(fit.standard_parallel_2),
            "rms_pixel_error": float(fit.rms_pixel_error),
            "max_pixel_error": float(fit.max_pixel_error),
            "pixel_to_proj_x_coeffs": tuple(float(value) for value in np.asarray(pixel_map.pixel_to_proj_x_coeffs, dtype=np.float64)),
            "pixel_to_proj_y_coeffs": tuple(float(value) for value in np.asarray(pixel_map.pixel_to_proj_y_coeffs, dtype=np.float64)),
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
                longitude_of_projection_origin=np.asarray(summary["longitude_of_projection_origin"], dtype=np.float64),
                latitude_of_projection_origin=np.asarray(summary["latitude_of_projection_origin"], dtype=np.float64),
                standard_parallel_1=np.asarray(summary["standard_parallel_1"], dtype=np.float64),
                standard_parallel_2=np.asarray(summary["standard_parallel_2"], dtype=np.float64),
                rms_pixel_error=np.asarray(summary["rms_pixel_error"], dtype=np.float64),
                max_pixel_error=np.asarray(summary["max_pixel_error"], dtype=np.float64),
                pixel_to_proj_x_coeffs=np.asarray(pixel_map.pixel_to_proj_x_coeffs, dtype=np.float64),
                pixel_to_proj_y_coeffs=np.asarray(pixel_map.pixel_to_proj_y_coeffs, dtype=np.float64),
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
        print(f"  lon0/lat0: {summary['longitude_of_projection_origin']}, {summary['latitude_of_projection_origin']}")
        print(f"  standard parallels: {summary['standard_parallel_1']}, {summary['standard_parallel_2']}")
        print(f"  rms pixel error: {summary['rms_pixel_error']:.6g}")
        print(f"  max pixel error: {summary['max_pixel_error']:.6g}")
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
