from __future__ import annotations

import argparse
import os
from pathlib import Path

from osgeo import gdal, ogr

gdal.UseExceptions()

DEFAULT_INPUT_SHP = "water_polygons.shp"
DEFAULT_BASE_TILE_ROOT = "water_tiles"
RESOLUTION_TO_OUTPUT_DIR = {
    125: "water_tiles_125m",
    250: "water_tiles_250m",
    500: "water_tiles_500m",
}


def _cleanup_existing_variants(base_path: Path) -> None:
    for suffix in (".tif", ".0", ".1"):
        candidate = base_path.with_suffix(suffix)
        if candidate.exists():
            candidate.unlink()


def _write_marker(marker_path: Path) -> None:
    marker_path.parent.mkdir(parents=True, exist_ok=True)
    marker_path.write_bytes(b"")


def _resolution_config(resolution_m: int) -> tuple[int, int, float, str]:
    if resolution_m == 125:
        return 32, 16, 4.0, RESOLUTION_TO_OUTPUT_DIR[125]
    if resolution_m == 250:
        return 16, 8, 2.0, RESOLUTION_TO_OUTPUT_DIR[250]
    if resolution_m == 500:
        return 8, 4, 1.0, RESOLUTION_TO_OUTPUT_DIR[500]
    raise ValueError("resolution_m must be 125, 250, or 500")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rasterize OSM water polygons into tile masks.")
    parser.add_argument(
        "--input-shp",
        default=DEFAULT_INPUT_SHP,
        help="Input water_polygons shapefile.",
    )
    parser.add_argument(
        "--resolution-m",
        type=int,
        choices=(125, 250, 500),
        default=125,
        help="Target raster resolution in meters.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Defaults to water_tiles_<resolution>m.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    cols, rows, px_scale, default_output_dir = _resolution_config(int(args.resolution_m))
    output_dir = str(args.output_dir or default_output_dir)
    os.makedirs(output_dir, exist_ok=True)

    ds_shp = ogr.Open(str(args.input_shp))
    layer = ds_shp.GetLayer()
    extent = layer.GetExtent()  # (min_x, max_x, min_y, max_y)
    _s_min_x, _s_max_x, s_min_y, s_max_y = extent
    spatial_ref = layer.GetSpatialRef().ExportToWkt()

    px_per_deg = (86400 / 360) * px_scale
    deg_x = 360 / cols
    deg_y_step = 180 / rows

    for r in range(rows):
        for c in range(cols):
            min_x = -180 + (c * deg_x)
            max_x = min_x + deg_x
            base_path = Path(output_dir) / f"tile_y{r}_x{c}"
            output_path = base_path.with_suffix(".tif")

            tile_max_y = 90 - (r * deg_y_step)
            tile_min_y = tile_max_y - deg_y_step
            current_min_y = max(tile_min_y, s_min_y)
            current_max_y = min(tile_max_y, s_max_y)
            if current_min_y >= current_max_y:
                generated_path = base_path.with_suffix(".0")
                _write_marker(generated_path)
                print(
                    f"Generating {generated_path.name}: 0x0 (Lat: {tile_min_y:.2f} to {tile_max_y:.2f})"
                )
                continue

            tile_w = int((max_x - min_x) * px_per_deg)
            tile_h = int((current_max_y - current_min_y) * px_per_deg)
            if tile_w <= 0 or tile_h <= 0:
                generated_path = base_path.with_suffix(".0")
                _write_marker(generated_path)
                print(
                    f"Generating {generated_path.name}: {tile_w}x{tile_h} (Lat: {current_min_y:.2f} to {current_max_y:.2f})"
                )
                continue

            _cleanup_existing_variants(base_path)

            layer.ResetReading()
            layer.SetSpatialFilterRect(min_x, current_min_y, max_x, current_max_y)

            feature_count = layer.GetFeatureCount()
            if feature_count == 0:
                generated_path = base_path.with_suffix(".0")
                _write_marker(generated_path)
                print(
                    f"Generating {generated_path.name}: {tile_w}x{tile_h} (Lat: {current_min_y:.2f} to {current_max_y:.2f})"
                )
                layer.SetSpatialFilter(None)
                continue

            mem_driver = gdal.GetDriverByName("MEM")
            raster = mem_driver.Create("", tile_w, tile_h, 1, gdal.GDT_Byte)
            raster.SetGeoTransform(
                (
                    min_x,
                    (max_x - min_x) / tile_w,
                    0.0,
                    current_max_y,
                    0.0,
                    -(current_max_y - current_min_y) / tile_h,
                )
            )
            raster.SetProjection(spatial_ref)
            band = raster.GetRasterBand(1)
            band.Fill(0)
            gdal.RasterizeLayer(raster, [1], layer, burn_values=[1])

            array = band.ReadAsArray()
            generated_path = output_path
            if array is not None and array.size > 0:
                if (array == 0).all():
                    generated_path = base_path.with_suffix(".0")
                    _write_marker(generated_path)
                elif (array != 0).all():
                    generated_path = base_path.with_suffix(".1")
                    _write_marker(generated_path)
                else:
                    gdal.Translate(
                        str(output_path),
                        raster,
                        format="GTiff",
                        creationOptions=["COMPRESS=CCITTFAX4", "NBITS=1", "PHOTOMETRIC=MINISBLACK"],
                    )
            else:
                gdal.Translate(
                    str(output_path),
                    raster,
                    format="GTiff",
                    creationOptions=["COMPRESS=CCITTFAX4", "NBITS=1", "PHOTOMETRIC=MINISBLACK"],
                )
            print(
                f"Generating {generated_path.name}: {tile_w}x{tile_h} (Lat: {current_min_y:.2f} to {current_max_y:.2f})"
            )
            layer.SetSpatialFilter(None)

    print("Done.")


if __name__ == "__main__":
    main()
