# dev-samples

Exploratory scripts and one-off investigation helpers live here.

These files are not part of the main application surface. They are useful for
debugging, data inspection, and reproducible experiments.

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
- Converts a color geosatellite-style image into a grayscale proxy that
  emphasizes bright, low-saturation cloud regions while suppressing blue sea,
  green land, and the top-left logo area by default.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/geo_satellite_cloud_proxy.py input.png -o output.png
```

## Geo-satellite cloud proxy from 3 frames

- `geo_satellite_cloud_proxy_temporal.py`
- Combines 3 images, suppresses temporally static overlays, and then builds a
  grayscale cloud proxy. This is useful for removing thin map borders and
  other persistent line art.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/geo_satellite_cloud_proxy_temporal.py \
  Europe-IR-20260525130000.png \
  Europe-IR-20260525143000.png \
  Europe-IR-20260525174500.png \
  -o Europe-IR-proxy.png
```

## Geo-satellite static RGB mask

- `geo_satellite_static_rgb_mask.py`
- Extracts bright line-like pixels that remain in the same place across 3
  images by subtracting a local background first. This is useful for isolating
  map overlays or other static line art before further processing.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/geo_satellite_static_rgb_mask.py \
  Europe-IR-20260525130000.png \
  Europe-IR-20260525143000.png \
  Europe-IR-20260525174500.png \
  -o Europe-IR-static-mask.png
```

## Meteosat / geosatellite API pixel map fit

- `fit_geostationary_image_mapping.py`
- Fits a geostationary pixel-to-lat/lon transform from a Meteosat or other
  geosatellite API image plus a `latlonmap.txt`-style control-point file. It
  locates the exact RGB marker pixels, solves the affine transform from
  projection space to image pixels, and can optionally dump full-resolution
  lon/lat grids to `.npz`.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/fit_geostationary_image_mapping.py \
  meteosat.png \
  latlonmap.txt \
  --dump-grid meteosat_lonlat.npz \
  --query 640,360
```

## Geostationary lat/lon grid overlay

- `draw_geostationary_latlon_grid.py`
- Draws a 10-degree lat/lon grid on top of a geostationary image. Lines at
  30-degree multiples are red; the rest are black. It can read a precomputed
  `.npz` lon/lat grid or fit the image from a companion `latlonmap.txt`.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/draw_geostationary_latlon_grid.py \
  meteosat.png \
  -o meteosat_grid.png
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
