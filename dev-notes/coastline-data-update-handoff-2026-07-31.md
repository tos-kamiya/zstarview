# Coastline data update handoff

This note is a handoff for the next session. The goal is to replace the
optional coastline and sea-mask data with the latest OSM Water Polygons source.

## Current status

- The latest source archive has already been downloaded and extracted under
  `raw-data/`:
  `raw-data/water-polygons-split-4326-20260730/`.
- The source archive is
  `raw-data/water-polygons-split-4326-20260730.zip`.
- Source date: `2026-07-30T00:00:00Z`.
- Source contents: 53,330 Polygon features, WGS84 / EPSG:4326.
- Archive SHA-256:
  `2f705bae3f7a7af7a7ef4740cfd6621ce8004df0aaf7fe01d914ff1498afdf17`
- The bundled coastline release assets and 25m sea-mask assets have not yet
  been replaced.
- `src/zstarview/data/earth_guide_land_110m.json` must remain unchanged. It
  is based on Natural Earth 1:110m land data, not OSM Water Polygons.

## Important files

- Procedure:
  `src/zstarview/data/earth_guide_land_regeneration.md`
- Existing coastline tile builder:
  `dev-samples/build_coastline_vector_tiles.py`
- Existing large-tile splitter:
  `dev-samples/split_large_coastline_tiles.py`
- Existing release packager:
  `dev-samples/package_coastline_release.py`
- Existing sea-mask generator:
  `src/zstarview/data/osm_water_polygons_to_tiff.py`
- Current decision log:
  `dev-notes/session-2026-07-31.md`

## Planned script organization

The user agreed that the data-maintenance scripts may initially be copied,
rather than moved, from `dev-samples/` into a dedicated directory such as
`src/zstarview/data/tools/`. This keeps existing commands working while the
new path becomes canonical in the procedure document.

The likely files to copy are:

- `build_coastline_vector_tiles.py`
- `split_large_coastline_tiles.py`
- `package_coastline_release.py`
- `osm_water_polygons_to_tiff.py`

Avoid maintaining two divergent implementations. After the new copies are
validated, either make the old files thin compatibility wrappers or explicitly
declare one location as the source of truth.

## Environment constraint

The existing coastline and raster scripts import `osgeo.ogr` and
`osgeo.gdal`. The current environment did not provide a usable system GDAL
installation, so the regeneration pipeline has not been run successfully.
An exploratory `uv run --with geopandas --with pyogrio ...` process was left
running in a Codex background terminal at the time of handoff; inspect it
with `/ps` and stop it with `/stop` if it is no longer needed.

## Suggested next steps

1. Check `/ps` and stop any stale background data-reading process.
2. Decide whether to install/provide GDAL for the current environment, or
   implement and validate a replacement using the available geospatial stack.
3. Copy the maintenance scripts into `src/zstarview/data/tools/` and update
   `earth_guide_land_regeneration.md` to use that path.
4. Generate coastline parent tiles from the 2026-07-30 Shapefile.
5. Inspect parent tiles and split the largest tiles if necessary.
6. Package and validate the coastline release assets before changing release
   tags or downloader defaults.
7. Generate 25m sea-mask tiles into `build/`, validate them, then promote the
   results to the bundled data directory only after review.
8. Record source date, generation date, release version, tile counts, and
   checksums in the session log and relevant provenance documentation.

## Commands

The source Shapefile is expected at:

```text
raw-data/water-polygons-split-4326-20260730/water-polygons-split-4326/water_polygons.shp
```

The procedure document contains the current generation commands. Replace
`dev-samples/` with `src/zstarview/data/tools/` after the scripts are copied.

## Working tree caution

At handoff, the intentional repository changes are:

- modified `dev-notes/session-2026-07-31.md`;
- new `src/zstarview/data/earth_guide_land_regeneration.md`;
- this handoff note.

Preserve unrelated user changes in the working tree. Run `git diff --check`
before committing.
