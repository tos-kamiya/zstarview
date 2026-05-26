# Geo-satellite Sample Index

This page lists the exploratory Geo-satellite-related scripts in `dev-samples/`
and the inputs and outputs they expect.

## Recommended pipeline

For the current Europe workflow, the default artifacts are:

- `raw-data/geosatellite/Europe-IR-gray-common-mask.png`
- `raw-data/geosatellite/eqdc_lonlat.npz`

The recommended sequence is:

1. Start from one or more `Europe-IR-*.png` frames.
2. Build a cloud proxy from a single frame with `geo_satellite_cloud_proxy.py`.
3. Build or refresh the shared gray-common mask with
   `geo_satellite_gray_common_mask.py`.
4. Fill the mask into the cloud proxy with `geo_satellite_cloud_inpaint.py`.
5. Project the inpainted cloud image into the observer-centric disc with
   `geo_satellite_cloudimage.py`.
6. Optionally overlay the lon/lat grid with
   `draw_equidistant_conic_latlon_grid.py`.

The inpaint step defaults to `raw-data/geosatellite/Europe-IR-gray-common-mask.png`, and the
cloudimage step defaults to `raw-data/geosatellite/eqdc_lonlat.npz`.

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
  raw-data/geosatellite/Europe-IR-20260525130000.png \
  raw-data/geosatellite/Europe-IR-20260525174500.png \
  raw-data/geosatellite/Europe-IR-20260525204500.png \
  -o raw-data/geosatellite/Europe-IR-gray-common-mask.png
```

## Cloud inpaint from a cloud image and a mask

- Script: `geo_satellite_cloud_inpaint.py`
- Purpose: fill masked regions in a grayscale cloud image by propagating
  nearby pixel values inward.
- Input:
  - one grayscale cloud image, such as the output of `geo_satellite_cloud_proxy.py`
  - one grayscale mask image, such as `raw-data/geosatellite/Europe-IR-gray-common-mask.png`
- Output:
  - one grayscale PNG with the masked regions filled
- Notes:
  - This is the next step after cloud proxy + mask extraction.
  - White mask pixels are treated as holes that should be filled from the
    surrounding cloud image.
  - The default mask path is `raw-data/geosatellite/Europe-IR-gray-common-mask.png`.

Example:

```bash
uv run -p .venv/bin/python dev-samples/geo_satellite_cloud_inpaint.py \
  Europe-IR-cloud-proxy.png \
  raw-data/geosatellite/Europe-IR-gray-common-mask.png \
  -o Europe-IR-cloud-inpainted.png
```

## Cloudimage projection from a cloud image

- Script: `geo_satellite_cloudimage.py`
- Purpose: map a grayscale cloud image back into an observer-centric cloud
  disc using the fitted lon/lat grid from the Europe image workflow.
- Input:
  - one grayscale cloud image, such as the output of `geo_satellite_cloud_inpaint.py`
  - one observer spec in the form `@lat,lon`
  - one fitted grid `.npz`, typically `raw-data/geosatellite/eqdc_lonlat.npz`
- Output:
  - one grayscale PNG cloud disc
- Notes:
  - Uses the observer's `lat`, `lon`, `alt`, `az`, and `fov` in the same
    style as `zstarview.cloudimage`.
  - The fitted grid `.npz` supplies the Europe image geometry needed to map
    lon/lat back into the source cloud image.
  - `--cloud-height-km` sets the cloud-top height above the Earth's surface,
    not the final shell radius.
  - The observer location must stay inside the Europe workflow band:
    latitude `32N` to `73N`, longitude `15W` to `35E`.
  - The default grid path is `raw-data/geosatellite/eqdc_lonlat.npz`.
  - This is a `dev-samples` validation helper, not a core app test target.

Example:

```bash
uv run -p .venv/bin/python dev-samples/geo_satellite_cloudimage.py \
  Europe-IR-cloud-inpainted.png \
  '@51.5074,-0.1278' \
  --alt 25 \
  --az 180 \
  --fov 60 \
  --cloud-height-km 5 \
  -o Europe-IR-cloud-disc.png
```

## Geostationary mapping fit

- Script: `fit_equidistant_conic_image_mapping.py`
- Purpose: estimate an Equidistant Conic style pixel-to-lat/lon mapping from a
  geostationary-style image and a control-point list.
- Input:
  - one source image, such as `raw-data/geosatellite/Europe-IR.png`
  - one control-point file, such as `raw-data/geosatellite/latlonmap.txt`
- Output:
  - a fitted lon/lat grid file when `--dump-grid` is used
- Notes:
  - Used to infer the geometry of the source image.

Example:

```bash
uv run -p .venv/bin/python dev-samples/fit_equidistant_conic_image_mapping.py \
  raw-data/geosatellite/Europe-IR.png \
  raw-data/geosatellite/latlonmap.txt \
  --dump-grid raw-data/geosatellite/eqdc_lonlat.npz
```

## Geostationary grid overlay

- Script: `draw_equidistant_conic_latlon_grid.py`
- Purpose: draw a latitude/longitude grid on top of a geostationary-style
  image using a precomputed lon/lat grid.
- Input:
  - one source image, such as `raw-data/geosatellite/Europe-IR.png`
  - one `--grid-npz` file produced by the fit script
- Output:
  - one PNG image with the grid overlay
- Notes:
  - Useful for checking whether the inferred projection matches the source
    image geometry.

Example:

```bash
uv run -p .venv/bin/python dev-samples/draw_equidistant_conic_latlon_grid.py \
  raw-data/geosatellite/Europe-IR.png \
  --grid-npz raw-data/geosatellite/eqdc_lonlat.npz \
  -o raw-data/geosatellite/eqdc_grid.png
```
