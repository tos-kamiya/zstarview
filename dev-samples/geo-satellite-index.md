# Geo-satellite Sample Index

This page lists the exploratory Geo-satellite-related scripts in `dev-samples/`
and the inputs and outputs they expect.

## Cloud proxy from one frame

- Script: `geo_satellite_cloud_proxy.py`
- Purpose: build a grayscale cloud-likelihood proxy from a single color
  Geo-satellite-style image.
- Input:
  - one PNG image, such as `Europe-IR-*.png`
- Output:
  - one grayscale PNG cloud proxy
- Notes:
  - Designed for white-ish cloud extraction from one frame.
  - Suppresses blue sea, green land, and the top-left logo area by default.

Example:

```bash
uv run -p .venv/bin/python dev-samples/geo_satellite_cloud_proxy.py input.png -o output.png
```

## Gray-common boundary mask from multiple frames

- Script: `geo_satellite_gray_common_mask.py`
- Purpose: keep pixels that are gray-ish and not too white across several
  images, so boundary-like line art remains while sea/land colors are removed.
- Input:
  - 2 or more PNG images, such as several `Europe-IR-*.png` frames
- Output:
  - one grayscale PNG mask
- Notes:
  - This is the current line-art extraction path.
  - It filters out strong sea/land color first, then intersects the remaining
    gray-like pixels across all frames.
  - The result can be reused as a boundary mask or as a guide for later
    fill/inpaint experiments.

Example:

```bash
uv run -p .venv/bin/python dev-samples/geo_satellite_gray_common_mask.py \
  Europe-IR-20260525130000.png \
  Europe-IR-20260525174500.png \
  Europe-IR-20260525204500.png \
  -o Europe-IR-gray-common-mask.png
```

## Geostationary mapping fit

- Script: `fit_equidistant_conic_image_mapping.py`
- Purpose: estimate an Equidistant Conic style pixel-to-lat/lon mapping from a
  geostationary-style image and a control-point list.
- Input:
  - one source image, such as `raw-data/Europe-IR.png`
  - one control-point file, such as `raw-data/latlonmap.txt`
- Output:
  - a fitted lon/lat grid file when `--dump-grid` is used
- Notes:
  - Used to infer the geometry of the source image.

Example:

```bash
uv run -p .venv/bin/python dev-samples/fit_equidistant_conic_image_mapping.py \
  raw-data/Europe-IR.png \
  raw-data/latlonmap.txt \
  --dump-grid raw-data/eqdc_lonlat.npz
```

## Geostationary grid overlay

- Script: `draw_equidistant_conic_latlon_grid.py`
- Purpose: draw a latitude/longitude grid on top of a geostationary-style
  image using a precomputed lon/lat grid.
- Input:
  - one source image, such as `raw-data/Europe-IR.png`
  - one `--grid-npz` file produced by the fit script
- Output:
  - one PNG image with the grid overlay
- Notes:
  - Useful for checking whether the inferred projection matches the source
    image geometry.

Example:

```bash
uv run -p .venv/bin/python dev-samples/draw_equidistant_conic_latlon_grid.py \
  raw-data/Europe-IR.png \
  --grid-npz raw-data/eqdc_lonlat.npz \
  -o raw-data/eqdc_grid.png
```
