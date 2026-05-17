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
  read statistics for the sea mask around an observer location.

Run it with:

```bash
uv run -p .venv/bin/python dev-samples/water_tile_overview.py --lat 0 --lon 135
```

Use `--observer-ground-m 8651` to inspect a high mountain site such as
Everest. Add `--ray-scan` to show the full ray-scan summary that mirrors the
water mask sampling path.
