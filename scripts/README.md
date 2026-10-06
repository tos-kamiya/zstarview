# Rendering scripts

## Experimental cloud voxel rendering

Render a sky image with the experimental cloud voxels:

```sh
uv run -p .venv/bin/python scripts/prototype_cloud_cell_shading.py \
  @35.68,139.69 \
  --view-alt 45 \
  --view-az 180 \
  --opacity 0.15 \
  --sky-opacity 0.4 \
  --size 1200x800 \
  --output /tmp/cloud-prototype.png
```

The location may be a city name or a coordinate accepted by zstarview. The
script uses the current UTC time by default; pass
`--datetime 2026-10-04T03:00:00Z` to reproduce a render. It saves a JSON
sidecar next to the PNG with the render time, selected satellite observation
time, source, view, Sun position, and model coefficients. `--fov` sets the
edge field of view; `--content-fov` controls the projection overscan.

The base image includes the sky-color disc and stars. The Sun, Moon, planets,
and unrelated overlay layers are hidden. Cloud data uses the existing
CloudDisc cache and satellite retrieval path. The satellite pixel brightness
temperatures are converted to estimated cloud amounts; those amounts and the
height allocation are model approximations, not measured water content or
optical depth.

`--opacity` adjusts cloud transmission from 0 to 1. Select
`--cloud-model altaz` for the previous angular cell experiment. Its default is
an 18-by-72 alt/az grid (5-degree cells), with nine shell fields from 3 to
11 km shaded separately. `--show-grid` overlays its angular boundaries,
including clear cells.

The default `--cloud-model voxel` uses each native satellite pixel as
one horizontal column, split into nine 1-km slabs (centers 3 through 11 km).
Raw pixel values are used without bilinear interpolation or blur. Rays traverse
voxel faces, and transmission depends on distance inside each voxel. Each voxel
gets a color from direct sunlight reaching it, replenishment by environmental
light, and a sunlight-to-whiteness curve. Cloud amount and path length determine
sunlight and view-ray extinction. Sky color and stars are attenuated by the
accumulated view-ray transmission.

Use `--cloud-only` to render the voxel clouds over a black background without
compositing sky color, stars, planets, or other base-image layers.
Use `--sky-color-only` to composite the atmospheric sky-color disc under the
clouds while omitting stars, planets, and other base-image layers. This mode
keeps the sky color visible without restoring stars.
For a reversible cloud-amount experiment, `--cloud-amount-subtract 0.1`
subtracts 0.1 from every estimated cloud amount and clamps at zero. Omitting
the option keeps the usual low-cloud suppression.

This experiment uses a local tangent volume with a 200-km horizontal extent.
By default, it intersects the satellite pixel-center and neighbor rays with
each altitude-offset Earth ellipsoid, then fits one affine horizontal grid per
1-km layer. Each layer remains a stack of rectangular prisms; the footprint
changes at layer boundaries instead of tapering within a voxel. This captures
approximate pixel shift and size changes while avoiding individual frustum
intersection tests. It still approximates each layer locally and omits B16
redistribution. Use `--flat-height-grid` to compare with the previous fixed
grid. JSON records each layer's offset and pixel basis, along with native pixel
spacing, cropped window, thresholds and height edges.
`--show-grid` draws dark edges near cloud voxel face boundaries in this mode;
it does not overlay the angular grid on clear sky. The old renderer remains
available with `--cloud-model altaz`; `--alt-bins` and `--az-bins` apply only
there. The voxel implementation is in `src/zstarview/render/cloud_voxels.py`.

Use `--sky-opacity 0.5` to set the base sky-color opacity independently of
cloud opacity (`--opacity`). The sky setting is passed to the image exporter
and recorded in the JSON sidecar.

The clear-sky sunlight baseline is `max(0, sin(solar_altitude))`. Along the
Sun ray, cloud extinction reduces this light while environmental light adds
5 percent of the baseline. The bounded formula `baseline*T + environment*(1-T)`
keeps deeply shaded daytime clouds from becoming black or fully white. A
transfer table maps sunlight `[0, 0.01, 0.03, 0.10, 0.30, 1]` to whiteness
`[0.18, 0.60, 0.85, 0.95, 0.99, 1]`; the 0.18 floor also sets the no-sun
appearance. These are visual tuning values, not radiometric measurements.
Cloud opacity defaults to 0.85; sky opacity defaults to 1.0. Set them
independently with `--opacity` and `--sky-opacity`. The JSON sidecar records
both opacities, the transfer table, environment fraction, Sun position, and
satellite observation time.
