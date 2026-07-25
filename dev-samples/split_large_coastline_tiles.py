"""Split the largest coastline GeoJSON tiles into 4x4 child tiles."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from osgeo import ogr


GRID_COLS = 32
GRID_ROWS = 16
PARENT_PATTERN = re.compile(r"tile_y(\d+)_x(\d+)\.geojson$")


def _box(min_lon: float, min_lat: float, max_lon: float, max_lat: float) -> ogr.Geometry:
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(min_lon, min_lat)
    ring.AddPoint(max_lon, min_lat)
    ring.AddPoint(max_lon, max_lat)
    ring.AddPoint(min_lon, max_lat)
    ring.AddPoint(min_lon, min_lat)
    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)
    return polygon


def _parent_bounds(row: int, col: int) -> tuple[float, float, float, float]:
    width = 360.0 / GRID_COLS
    height = 180.0 / GRID_ROWS
    min_lon = -180.0 + col * width
    max_lat = 90.0 - row * height
    return min_lon, max_lat - height, min_lon + width, max_lat


def _cleanup(base: Path) -> None:
    for suffix in (".geojson", ".0"):
        path = base.with_suffix(suffix)
        if path.exists():
            path.unlink()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--count", type=int, default=8)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    parents = sorted(input_dir.glob("tile_y*_x*.geojson"), key=lambda path: path.stat().st_size, reverse=True)[: args.count]
    driver = ogr.GetDriverByName("GeoJSON")

    for parent_path in parents:
        match = PARENT_PATTERN.fullmatch(parent_path.name)
        if match is None:
            continue
        row, col = int(match.group(1)), int(match.group(2))
        min_lon, min_lat, max_lon, max_lat = _parent_bounds(row, col)
        parent_source = ogr.Open(str(parent_path))
        if parent_source is None:
            raise RuntimeError(f"cannot open {parent_path}")
        parent_layer = parent_source.GetLayer()
        parent_features = [feature.Clone() for feature in parent_layer]
        for child_row in range(4):
            child_min_lat = min_lat + (max_lat - min_lat) * child_row / 4.0
            child_max_lat = min_lat + (max_lat - min_lat) * (child_row + 1) / 4.0
            for child_col in range(4):
                child_min_lon = min_lon + (max_lon - min_lon) * child_col / 4.0
                child_max_lon = min_lon + (max_lon - min_lon) * (child_col + 1) / 4.0
                base = output_dir / f"tile_y{row}_x{col}_q{child_row}{child_col}"
                _cleanup(base)
                child_geom = _box(child_min_lon, child_min_lat, child_max_lon, child_max_lat)
                output_path = base.with_suffix(".geojson")
                datasource = None
                output_layer = None
                feature_count = 0
                for parent_feature in parent_features:
                    geometry = parent_feature.GetGeometryRef()
                    if geometry is None:
                        continue
                    clipped = geometry.Intersection(child_geom)
                    if clipped is None or clipped.IsEmpty():
                        continue
                    if datasource is None:
                        datasource = driver.CreateDataSource(str(output_path))
                        output_layer = datasource.CreateLayer("coastline", None, ogr.wkbUnknown)
                    output_feature = ogr.Feature(output_layer.GetLayerDefn())
                    output_feature.SetGeometry(clipped)
                    output_layer.CreateFeature(output_feature)
                    output_feature = None
                    feature_count += 1
                output_layer = None
                datasource = None
                if feature_count == 0:
                    if output_path.exists():
                        output_path.unlink()
                    base.with_suffix(".0").write_bytes(b"")
                print(f"{parent_path.name} -> {base.name}: features={feature_count}")


if __name__ == "__main__":
    main()
