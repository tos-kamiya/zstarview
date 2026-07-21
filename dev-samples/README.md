# dev-samples

Exploratory scripts and one-off investigation helpers live here.

These files are not part of the main application surface. They are useful for
debugging, data inspection, and reproducible experiments.

## Overture Transportation road probe

- `overture_transportation_probe.py`
- Downloads a small Overture `segment` GeoJSON extract around a point and
  reports road-segment counts, `width_rules` coverage, width values, and road
  classes with width data.
- The default bbox radius is 50m. The downloaded file is temporary unless
  `--output` is supplied.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/overture_transportation_probe.py \
  --lat 35.681236 --lon 139.767125
```

Keep the raw probe result for inspection with `--output path/to/roads.geojson`.

## EOG VNL night-light tile builder

- `build_vnl_night_lights.py`
- Converts an EOG VIIRS Nighttime Lights VNL v2.2 annual GeoTIFF into the
  eight 90-degree EPSG:4326 GeoTIFF tiles used by zstarview.
- Accepts either a locally downloaded `.tif`/`.tif.gz` file or a URL. The
  source is not copied into the repository. The output also includes
  `manifest.json` with tile dimensions, bounds, SHA-256 checksums, source
  metadata, and the EOG CC BY 4.0 attribution. It also writes
  `NOTICE-night-lights.txt` for inclusion with the release assets.
- The default source band is band 1. Select the masked radiance band from the
  EOG file before running the conversion if the downloaded file contains
  multiple bands.

Run it with a local source file:

```bash
uv run -p .venv/bin/python dev-samples/build_vnl_night_lights.py \
  --source VNL_npp_2025_global_vcmslcfg_v2_c202604011200.average_masked.dat.tif.gz \
  --year 2025 \
  --output-dir build/night-lights/2025
```

An authenticated URL can be used without putting the token in the command
line history:

```bash
export EOG_BEARER_TOKEN='...'
uv run -p .venv/bin/python dev-samples/build_vnl_night_lights.py \
  --url 'https://eogdata.mines.edu/...' \
  --bearer-token-env EOG_BEARER_TOKEN \
  --output-dir build/night-lights/2020
```

Do not commit the raw download, bearer token, or generated raster directory.

## Night-light boundary scan

`find_night_light_boundaries.py` scans VNL GeoTIFF tiles for a roughly 128 km
square containing a strong brightness difference across a straight boundary.
The boundary direction is tested at 0, 45, 90, and 135 degrees. The result
reports latitude, longitude, angle, and the two brightness sums, in the same
coordinate order used by zstarview.

```bash
uv run -p .venv/bin/python dev-samples/find_night_light_boundaries.py \
  build/night-lights/2025/C1.tif
```

The default output contains ten spatially separated candidates; candidates
must be at least 128 km apart. Use `--top 20` to print more candidates or
`--min-separation-km 250` to spread them farther apart. The default scan uses
a 4 km working resolution and moves the window center by 16 km; use
`--sample-km` and `--step-km` to trade speed for detail.

Each result also includes a recommended observer location 60 km toward the
dark side, the bright/dark side azimuths, the absolute brightness difference,
and the bright-to-dark brightness ratio. Pass the recommended coordinates to
zstarview as `@recommended_latitude,recommended_longitude` and look toward the
reported bright-side azimuth at a low view altitude.

Limit the scan to a latitude/longitude rectangle by specifying both bounds:

```bash
uv run -p .venv/bin/python dev-samples/find_night_light_boundaries.py \
  build/night-lights/2025/A1.tif build/night-lights/2025/B1.tif \
  --bound-lat "24,50" --bound-lon=-125,-66
```

For the Geo-satellite-related helpers, see:

- [`geo-satellite-index.md`](geo-satellite-index.md)

Current default artifacts for the Europe workflow:

- `raw-data/geosatellite/Europe-IR-gray-common-mask.png`
- `raw-data/geosatellite/eqdc_lonlat.npz`

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

## Public ArcGIS hurricane feed probe

- `arcgis_active_hurricanes_probe.py`
- Fetches the public `Active_Hurricanes_v1` ArcGIS FeatureServer without an
  API key or token and prints the current position, forecast position, and
  wind polygon layers for the active storm.
- The `Observed Position` layer is queried in latest-DTG order so the script
  can show the current position explicitly.
- This is useful for checking whether a public tropical cyclone feed can be
  consumed directly before wiring up any authenticated ArcGIS workflow.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/arcgis_active_hurricanes_probe.py
```

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
  raw-data/geosatellite/Europe-IR-20260525130000.png \
  raw-data/geosatellite/Europe-IR-20260525174500.png \
  raw-data/geosatellite/Europe-IR-20260525204500.png \
  -o raw-data/geosatellite/Europe-IR-gray-common-mask.png
```

## Geo-satellite cloud inpaint

- `geo_satellite_cloud_inpaint.py`
- Fills masked regions in a grayscale cloud image by propagating nearby pixel
  values inward.
- This is the next step after building a cloud proxy and a gray-common mask:
  use the cloud image as the base, and use the mask image to mark line-art or
  boundary regions that should be filled from surrounding pixels.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/geo_satellite_cloud_inpaint.py \
  Europe-IR-cloud-proxy.png \
  raw-data/geosatellite/Europe-IR-gray-common-mask.png \
  -o Europe-IR-cloud-inpainted.png
```

## Geo-satellite cloudimage projection

- `geo_satellite_cloudimage.py`
- Projects a grayscale cloud image back into an observer-centric cloud disc
  using the fitted lon/lat grid from the Europe image workflow.
- Input:
  - one grayscale cloud image, such as the output of `geo_satellite_cloud_inpaint.py`
  - one observer spec in the form `@lat,lon`
  - the fitted grid `.npz` produced by the geometry-fit helper
- Output:
  - one grayscale PNG cloud disc with the source image mapped into the
    observer's `alt/az` view
- Constraints:
  - `@lat,lon` must fall inside the Europe workflow band: latitude `32N` to
    `73N`, longitude `15W` to `35E`.
- This script is intended for exploratory use under `dev-samples/`, not as a
  core app test target.
- The grid file defaults to `raw-data/geosatellite/eqdc_lonlat.npz`.

Run it with:

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

## Equidistant Conic fit probe

- `fit_equidistant_conic_image_mapping.py`
- Fits an Equidistant Conic projection against the same RGB control-point
  markers and reports the pixel residuals. Use this when you want to test
  whether the image geometry looks more like a conic map than a geostationary
  view.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/fit_equidistant_conic_image_mapping.py \
  raw-data/geosatellite/Europe-IR.png \
  raw-data/geosatellite/latlonmap.txt \
  --dump-grid raw-data/geosatellite/eqdc_lonlat.npz
```

## Equidistant Conic grid overlay

- `draw_equidistant_conic_latlon_grid.py`
- Draws the same 10-degree lat/lon grid on top of an image using a precomputed
  Equidistant Conic lon/lat grid from `fit_equidistant_conic_image_mapping.py`.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/draw_equidistant_conic_latlon_grid.py \
  raw-data/geosatellite/Europe-IR.png \
  --grid-npz raw-data/geosatellite/eqdc_lonlat.npz \
  -o raw-data/geosatellite/eqdc_grid.png
```
