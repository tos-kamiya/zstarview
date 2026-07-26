# Handoff: 2026-07-26

## Current state

- Repository: `/home/toshihiro/playground/zstarview`
- Branch: `dev`
- HEAD: `b0588b34 feat: expand tower viewpoint sources`
- `main` and `origin/main`: `f3575c4d`
- `origin/dev`: `94b10475`
- No Git commit or push was performed for the current water-boundary implementation.
- No `build_grid_water_tile.py` or `build_regional_water_tiles.py` process remains running.

## User-visible context

The user is improving water visibility in zstarview screenshots. Inland water polygons are sampled into dots, and a boundary polyline is now derived from simplified polygon rings and rendered with the same water color. Taipei 101 was restored to the local `dev` branch by cherry-picking `2a4ee5b6`, resulting in commit `b0588b34`; `resolve_tower_viewpoint("taipei 101")` resolves correctly.

The user supplied the refreshed source archive:

- `raw-data/water-polygons-split-4326-20260725T0338.zip` (about 861 MB)
- Contains WGS84/EPSG:4326 OSM Water Polygons, 53,328 polygon features, data date 2026-07-25.
- Expanded source is at `raw-data/water-polygons-split-4326/` and includes a 1.2 GB `water_polygons.shp`.
- Source page: https://osmdata.openstreetmap.de/data/water-polygons.html

## Implemented but uncommitted water-boundary changes

The following files contain the water boundary polyline implementation and must be preserved:

- `src/zstarview/water_overlay.py`
- `src/zstarview/gui/water_overlay_controller.py`
- `src/zstarview/gui/water_overlay_state.py`
- `src/zstarview/gui/window_render.py`
- `src/zstarview/gui/window_updates.py`
- `src/zstarview/render/render_types.py`
- `src/zstarview/render/terrain.py`
- `src/zstarview/render/atlas_pipeline.py`
- `src/zstarview/render/zstarview_pipeline.py`
- `src/zstarview/cli/export_image.py`
- `tests/test_water_overlay.py`

Behavior:

- `WaterOverlayPolyline` stores projected points from simplified polygon rings.
- Rings are clipped to the observer distance circle before projection.
- GUI and export paths pass and render polylines after water dots.
- Normal rendering uses the existing water theme color; boundary opacity is multiplied by 0.85.
- Persistent cache still stores simplified footprints; polylines are derived in memory.

Validation already completed before this handoff:

- Targeted water/controller/export/render tests: 76 passed.
- Full test suite: 1306 passed.
- Ruff and `git diff --check` passed for the touched implementation.
- Full mypy still reports pre-existing repository-wide errors, including existing water-overlay typing errors.

## 5m tile experiments

The user clarified the desired global-grid concept:

- Pixel resolution: 5m.
- Refine the existing 32x16 global tile grid by 5 in both axes.
- Total tiles: `32 * 5 * 16 * 5 = 12,800`.
- Each tile covers 2.25 degrees in longitude and latitude.

The earlier regional prototype used 0.02-degree tiles and is not the desired global-grid layout:

- `raw-data/water_tiles_5m_taipei/`
- `raw-data/water-tiles-5m-taipei-20260726.zip`
- 625 tiles, 577 empty markers, 48 partial GeoTIFFs, about 216 KB unpacked and 173 KB zipped.

The current global-grid Taipei prototype is the relevant one:

- `raw-data/water_tiles_5m_global_grid_taipei/tile_y28_x134.tif`
- Grid bounds: longitude 121.5-123.75, latitude 24.75-27.0.
- This tile contains Taipei 101.
- Raster dimensions: 45072 x 49759 pixels.
- Pixel sizes: about 5m in both directions in EPSG:4326.
- Features intersecting tile: 13.
- File size: about 1.5 MB due to 1-bit DEFLATE compression.
- It was generated without loading the full raster array into Python; GDAL writes directly to a tiled BigTIFF.

Generation helper:

- `dev-samples/build_grid_water_tile.py`
- Example:

  ```text
  /usr/bin/python3 dev-samples/build_grid_water_tile.py \
    --input-shp raw-data/water-polygons-split-4326/water_polygons.shp \
    --output-dir raw-data/water_tiles_5m_global_grid_taipei \
    --x 134 --y 28 --resolution-m 5
  ```

The helper uses the 5x-refined grid, EPSG:4326, 1-bit DEFLATE GeoTIFF, tiled blocks, and BIGTIFF=YES. Empty tiles are emitted as `.0` markers. It currently does not classify all-water tiles as `.1`; that can be added if needed.

## Important caveats for next work

- The generated 5m global-grid tile is an experiment only; the application does not load it yet.
- The app's current runtime sea-mask code expects the existing 125m/250m/500m roots and tile naming assumptions. Integrating the 5m root requires a new resolution/source selection path and likely a regional or manifest-based download policy.
- The 5m tile's geographic extent is large, but compression is effective for this sample. More coastal-complex tiles must be measured before extrapolating worldwide size.
- Do not delete or reset the existing uncommitted water-boundary changes. There are also unrelated user artifacts: `build/night-lights/2025/`, `shot1.png`, and `shot2.png`.
- `raw-data/` is ignored by Git, so generated data will not appear in normal Git status. The helper scripts may also be ignored by repository rules; check explicitly before committing.

## Likely next steps

1. Generate adjacent Taipei tiles, especially `x133/y28` and `x134/y27`, and compare compressed sizes.
2. Inspect the 5m tile visually or sample it through a temporary mask reader.
3. Decide whether all-water markers, a manifest, and release archive layout are needed.
4. If the 5m path is adopted, update the specification/design before wiring runtime downloads.
5. Commit the water-boundary implementation separately from any 5m data-generation tooling.
