# Geo-satellite Traceability

This directory collects the Europe Geo-satellite investigation assets used for
the gray-common-mask workflow and its supporting projection files.

The source Europe infrared frames were obtained from MET Norway Geo-Satellite
API:

`https://api.met.no/weatherapi/geosatellite/1.4/?area=europe&type=infrared`

## Mask generation

The refreshed gray-common mask was generated with:

```bash
uv run -p .venv/bin/python dev-samples/geo_satellite_gray_common_mask.py \
  raw-data/geosatellite/Europe-IR-20260525130000.png \
  raw-data/geosatellite/Europe-IR-20260526020000.png \
  raw-data/geosatellite/Europe-IR-20260526110000.png \
  raw-data/geosatellite/Europe-IR-20260526201500.png \
  raw-data/geosatellite/Europe-IR-20260527120000.png \
  raw-data/geosatellite/Europe-IR-20260528114500.png \
  raw-data/geosatellite/Europe-IR-20260604141500.png \
  -o raw-data/geosatellite/Europe-IR-gray-common-mask.png \
  --gray-spread 0.35 \
  --white-brightness 0.63
```

## Files in this directory

- `Europe-IR-20260525130000.png`
- `Europe-IR-20260526020000.png`
- `Europe-IR-20260526110000.png`
- `Europe-IR-20260526201500.png`
- `Europe-IR-20260527120000.png`
- `Europe-IR-20260528114500.png`
- `Europe-IR-20260604141500.png`
- `Europe-IR-gray-common-mask.png`
- `Europe-IR.png`
- `eqdc_grid.png`
- `eqdc_lonlat.npz`
- `latlonmap.txt`

## Notes

- `Europe-IR-gray-common-mask.png` is the generated output from the command
  above.
- `Europe-IR.png`, `eqdc_lonlat.npz`, `eqdc_grid.png`, and `latlonmap.txt`
  are the companion assets for the projection and grid-fitting workflow.
- The current gray-common mask uses seven source frames, including
  `Europe-IR-20260604141500.png`, to better match the shared cloudy-area
  footprint across the selected captures.
- The source frames in this directory were downloaded from the MET Norway
  Geo-Satellite API URL above.
