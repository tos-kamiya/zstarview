# Earth guide and coastline data update

This document records the current data split and the procedure for updating
the optional coastline overlay. The Earth guide and the coastline overlay use
different products from the same general map-data workflow.

This is a data-maintenance reference only. It is not used by the application
at runtime, and its commands are not run automatically. The application reads
the already bundled `earth_guide_land_110m.json`; coastline data, when enabled,
is loaded from the separately downloaded and versioned coastline cache. The
source files and build directories described below are not runtime inputs.

The directory roles are:

- `raw-data/`: externally downloaded source data;
- `dev-samples/`: conversion, inspection, and preview scripts;
- `src/zstarview/data/`: generated data bundled with the application.

Source data must not be stored under `dev-samples/`.

## Current data flow

### Earth guide

The bundled Earth guide data remains derived from the Natural Earth 1:110m land
source. The runtime data is planned as four LOD-specific JSON artifacts, with
LOD 0 being the smallest and LOD 3 the most detailed. The current single
`src/zstarview/data/earth_guide_land_110m.json` is the baseline artifact until
the LOD split is implemented.

The baseline is derived from the Natural Earth 1:110m land source:

- source: `raw-data/natural-earth/ne_110m_land.json`;
- source data date: 2017-11;
- simplification: `--simplify-deg 2.0 --min-ring-area-deg2 8.0`;
- zstarview runtime generation: 2026-04-12.

When the LOD split is implemented, all four artifacts must retain the same
coordinate convention (`lon_lat_deg`), ring metadata, seam handling, and
source attribution. Only simplification and small-ring retention may differ.
The intended runtime selection is based on `ScreenGeometry.radius`:

- LOD 0: `radius <= 350 px`;
- LOD 1: `350 px < radius <= 700 px`;
- LOD 2: `700 px < radius <= 1050 px`;
- LOD 3: `radius > 1050 px`.

`fast_mode` always selects LOD 0, regardless of window size. Generation and
promotion of the LOD artifacts must compare file size, ring/vertex counts,
small-island retention, polar and seam behavior, and representative rendered
output before replacing the baseline runtime asset.

The Earth guide is not regenerated from OSM Water Polygons. OSM provides a
separate land-polygon product, but replacing the current guide would require
additional handling for large rings, holes, antimeridian wrapping, and the
guide's coarse display resolution. Natural Earth remains the stable input.

### Optional coastline overlay

The coastline overlay uses OSM Water Polygons, which are ocean and sea
polygons derived from OpenStreetMap coastline ways. The latest source obtained
on 2026-07-31 is:

```text
raw-data/water-polygons-split-4326-20260730.zip
```

The package README reports:

- source date: 2026-07-30T00:00:00Z;
- 53,330 Polygon features;
- WGS84 / EPSG:4326;
- ODbL licensing.

The source archive SHA-256 is:

```text
2f705bae3f7a7af7a7ef4740cfd6621ce8004df0aaf7fe01d914ff1498afdf17
```

The source download URL is:

```text
https://osmdata.openstreetmap.de/download/water-polygons-split-4326.zip
```

## Updating the coastline source

Use a separate output directory and keep the downloaded source under
`raw-data/`:

```bash
mkdir -p raw-data/water-polygons-split-4326-YYYYMMDD
curl -L --fail --remote-time \
  --output raw-data/water-polygons-split-4326-YYYYMMDD.zip \
  https://osmdata.openstreetmap.de/download/water-polygons-split-4326.zip
unzip -q raw-data/water-polygons-split-4326-YYYYMMDD.zip \
  -d raw-data/water-polygons-split-4326-YYYYMMDD
```

Record the date, feature count, source URL, and SHA-256 from the archive and
its `README.txt` before generating derived data.

## Regenerating coastline vector tiles

The existing experimental pipeline reads the source Shapefile, clips its
polygon boundaries to the 32x16 world grid, and packages the resulting vector
tiles. Run it in a separate build directory:

```bash
uv run -p .venv/bin/python dev-samples/build_coastline_vector_tiles.py \
  --input-shp raw-data/water-polygons-split-4326-YYYYMMDD/water-polygons-split-4326/water_polygons.shp \
  --output-dir build/coastline-YYYYMMDD/parent \
  --gzip
```

If large parent tiles need additional subdivision:

```bash
uv run -p .venv/bin/python dev-samples/split_large_coastline_tiles.py \
  --input-dir build/coastline-YYYYMMDD/parent \
  --output-dir build/coastline-YYYYMMDD/children
```

Package the tiles and create the manifest and checksums:

```bash
uv run -p .venv/bin/python dev-samples/package_coastline_release.py \
  --parent-dir build/coastline-YYYYMMDD/parent \
  --child-dir build/coastline-YYYYMMDD/children \
  --output-dir build/coastline-YYYYMMDD/release \
  --version YYYYMMDD
```

Before publishing or changing the downloader's release tag, inspect the
manifest, verify all 32 longitude-column assets, compare representative
tiles with the previous release, and check the SHA-256 file.

## Regenerating sea-mask tiles

The same Shapefile can be rasterized for the sea-mask path:

```bash
uv run -p .venv/bin/python src/zstarview/data/osm_water_polygons_to_tiff.py \
  --input-shp raw-data/water-polygons-split-4326-YYYYMMDD/water-polygons-split-4326/water_polygons.shp \
  --resolution-m 25 \
  --output-dir build/water-YYYYMMDD/tiles-25m
```

The 25m global output is large. Generate and inspect it outside
`src/zstarview/data/`, then promote it only after checking the tile counts,
empty/full markers, boundary tiles, and representative coastal locations.

## Validation and promotion

For each source update:

- compare the source date and feature count with the previous package;
- inspect the antimeridian, polar areas, islands, and tile boundaries;
- verify that clipping produces intentionally split coastline segments and
  that no segment is silently lost at a tile boundary;
- compare SVG or rendered previews at Japan, Europe, and island-rich areas;
- record the source date, release version, generation date, and checksums.

Do not replace the bundled Earth guide JSON as part of a coastline update.
Only replace coastline release assets and sea-mask assets after their
corresponding validation is complete.

## References

- OSM Water Polygons:
  `https://osmdata.openstreetmap.de/data/water-polygons.html`
- OSM Water Polygons download:
  `https://osmdata.openstreetmap.de/download/water-polygons-split-4326.zip`
- Natural Earth 1:110m land source:
  `https://www.naturalearthdata.com/downloads/110m-physical-vectors/`
