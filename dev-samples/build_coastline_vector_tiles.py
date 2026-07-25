"""Build clipped coastline GeoJSON tiles on the existing 32x16 global grid."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

from osgeo import ogr


GRID_COLS = 32
GRID_ROWS = 16


def _tile_geometry(min_lon: float, min_lat: float, max_lon: float, max_lat: float) -> ogr.Geometry:
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(min_lon, min_lat)
    ring.AddPoint(max_lon, min_lat)
    ring.AddPoint(max_lon, max_lat)
    ring.AddPoint(min_lon, max_lat)
    ring.AddPoint(min_lon, min_lat)
    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)
    return polygon


def _cleanup(base: Path) -> None:
    for suffix in (".geojson", ".geojson.gz", ".0"):
        path = base.with_suffix(suffix)
        if path.exists():
            path.unlink()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-shp", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--gzip", action="store_true")
    args = parser.parse_args()

    source = ogr.Open(args.input_shp)
    if source is None:
        raise RuntimeError(f"cannot open {args.input_shp}")
    source_layer = source.GetLayer()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    geojson_driver = ogr.GetDriverByName("GeoJSON")

    tile_width = 360.0 / GRID_COLS
    tile_height = 180.0 / GRID_ROWS
    for row in range(GRID_ROWS):
        max_lat = 90.0 - row * tile_height
        min_lat = max_lat - tile_height
        for col in range(GRID_COLS):
            min_lon = -180.0 + col * tile_width
            max_lon = min_lon + tile_width
            base = output_dir / f"tile_y{row}_x{col}"
            _cleanup(base)
            source_layer.SetSpatialFilterRect(min_lon, min_lat, max_lon, max_lat)
            tile_geom = _tile_geometry(min_lon, min_lat, max_lon, max_lat)
            output_path = base.with_suffix(".geojson")
            datasource = None
            output_layer = None
            feature_count = 0
            for source_feature in source_layer:
                geometry = source_feature.GetGeometryRef()
                if geometry is None:
                    continue
                boundary = geometry.Boundary()
                clipped = boundary.Intersection(tile_geom)
                if clipped is None or clipped.IsEmpty():
                    continue
                if datasource is None:
                    datasource = geojson_driver.CreateDataSource(str(output_path))
                    if datasource is None:
                        raise RuntimeError(f"cannot create {output_path}")
                    output_layer = datasource.CreateLayer("coastline", source_layer.GetSpatialRef(), ogr.wkbUnknown)
                output_feature = ogr.Feature(output_layer.GetLayerDefn())
                output_feature.SetGeometry(clipped)
                output_layer.CreateFeature(output_feature)
                output_feature = None
                feature_count += 1
            output_layer = None
            datasource = None
            source_layer.SetSpatialFilter(None)
            if feature_count == 0:
                if output_path.exists():
                    output_path.unlink()
                marker = base.with_suffix(".0")
                marker.write_bytes(b"")
                print(f"empty tile: {marker.name}")
                continue
            if args.gzip:
                gzip_path = base.with_suffix(".geojson.gz")
                with output_path.open("rb") as source_file, gzip.open(gzip_path, "wb", compresslevel=9) as target_file:
                    shutil.copyfileobj(source_file, target_file)
                output_path.unlink()
                print(f"created {gzip_path}: features={feature_count}")
            else:
                print(f"created {output_path}: features={feature_count}")


if __name__ == "__main__":
    main()
