# Gaia EDR3 colour background

This directory contains the 2048x1024 display texture used as a possible
integrated stellar background. It is an equirectangular all-sky map in
Galactic coordinates, not a point-source catalog.

The source image is the Gaia EDR3 colour map from ESA/Gaia/DPAC, with
acknowledgement to A. Moitinho. The image is used under the CC BY-SA 3.0 IGO
or ESA Standard Licence. See `gaia-edr3-colour-manifest.json` for the source
and derived-asset metadata.

`prepare_gaia_background.py` can regenerate the packaged-size texture from
the downloaded source image:

```text
.venv/bin/python src/zstarview/data/gaia-edr3/prepare_gaia_background.py \
  The_colour_of_the_sky_from_Gaia_s_Early_Data_Release_3-copy.png
```
