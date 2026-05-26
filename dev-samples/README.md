# dev-samples

Exploratory scripts and one-off investigation helpers live here.

These files are not part of the main application surface. They are useful for
debugging, data inspection, and reproducible experiments.

For the Geo-satellite-related helpers, see:

- [`geo-satellite-index.md`](geo-satellite-index.md)

## Water tile intersection check

- `find_water_tile_intersections.py`
- Finds coarse-grid intersection points where the four surrounding water tiles
  are all `.tif` across the 125m, 250m, and 500m tile sets.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/find_water_tile_intersections.py
```

Add `--json` to emit machine-readable output.

## Water tile overview

- `water_tile_overview.py`
- Prints the tile grid, suffix counts, local probe samples, and per-band tile
  open statistics for the sea mask around an observer location.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/water_tile_overview.py --lat 0 --lon 135
```

Use `--observer-ground-m 8651` to inspect a high mountain site such as
Everest. Add `--ray-scan` to show the full ray-scan summary that mirrors the
water mask sampling path.

## Geo-satellite cloud proxy

- `geo_satellite_cloud_proxy.py`
- Converts a single color geosatellite-style image into a grayscale cloud
  likelihood proxy.
- This is the main exploratory path for white-ish cloud extraction from one
  frame.
- It emphasizes bright, low-saturation cloud regions while suppressing blue
  sea, green land, and the top-left logo area by default.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/geo_satellite_cloud_proxy.py input.png -o output.png
```

## Geo-satellite gray-common mask

- `geo_satellite_gray_common_mask.py`
- Extracts pixels that stay gray-ish and not too white across multiple images.
- This is the current line-art extraction path: it removes strong sea/land
  color first, then intersects the remaining gray-like pixels to keep the
  boundary-like structure that stays common across the provided frames.
- The output is a grayscale mask that can be reused as a boundary mask or as a
  guide for later fill/inpaint-style experiments.
- The default thresholds keep somewhat darker gray tones while rejecting more
  of the brightest whites.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/geo_satellite_gray_common_mask.py \
  Europe-IR-20260525130000.png \
  Europe-IR-20260525174500.png \
  Europe-IR-20260525204500.png \
  -o Europe-IR-gray-common-mask.png
```

## Equidistant Conic fit probe

- `fit_equidistant_conic_image_mapping.py`
- Fits an Equidistant Conic projection against the same RGB control-point
  markers and reports the pixel residuals. Use this when you want to test
  whether the image geometry looks more like a conic map than a geostationary
  view.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/fit_equidistant_conic_image_mapping.py \
  raw-data/Europe-IR.png \
  raw-data/latlonmap.txt \
  --dump-grid raw-data/eqdc_lonlat.npz
```

## Equidistant Conic grid overlay

- `draw_equidistant_conic_latlon_grid.py`
- Draws the same 10-degree lat/lon grid on top of an image using a precomputed
  Equidistant Conic lon/lat grid from `fit_equidistant_conic_image_mapping.py`.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/draw_equidistant_conic_latlon_grid.py \
  raw-data/Europe-IR.png \
  --grid-npz raw-data/eqdc_lonlat.npz \
  -o raw-data/eqdc_grid.png
```
