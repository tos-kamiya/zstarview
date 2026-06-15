# 2026-05-15 Handoff

## Scope

Water overlay preview and coastline handling.

## Current State

- The preview script and main app now fetch water data with a 1.2x query bbox.
- The SVG preview shows a dashed 1x bbox marker.
- The earlier coastline stitching experiments were abandoned.
- The current attempt to drop coastline polygons by bbox overlap did not solve the visual issue.

## What Was Tried

- Bbox stitching along the coastline ring.
- Inserting bbox corners along the perimeter arc.
- Switching to overfetching with a 1.2x query bbox.
- Drawing a 1x bbox marker in the SVG output.
- Dropping coastline polygons that do not overlap the view bbox.
- Switching the overlap basis between:
  - cache payload bbox
  - requested query bbox
  - current view bbox

## Key Observation

- The problematic coastline geometry still survives the overlap filter.
- That means the coastline polygon still intersects the current view bbox, so bbox-overlap is too weak as a removal criterion.
- The command used was correct:

```bash
uv run dev-samples/overpass_water_overlay_probe.py \
  --lat 34.6825 \
  --lon 135.1867 \
  --radius-km 5.0 \
  --input-cache ~/.cache/zstarview/water_overlay/earth_+34.6825_+135.1867_r39.05.json \
  --output-json tmp.json \
  --output-svg tmp.svg \
  --output-svg-raw tmp-raw.svg
```

## Files Most Relevant

- `dev-samples/overpass_water_overlay_probe.py`
- `src/zstarview/water_overlay.py`
- `src/zstarview/gui/water_overlay_controller.py`
- `src/zstarview/cli/export_image.py`
- `tests/test_overpass_water_overlay_probe.py`

## Suggested Next Step

- Decide whether coastline should be:
  - hidden entirely in the preview for this diagnostic path, or
  - clipped with a proper geometry operation, or
  - kept but rendered in a separate simplified way.
- If the goal is only to inspect rivers and inland water, the simplest next step is to suppress coastline rendering in the probe and keep the main app behavior unchanged.

## Validation Already Run

- `uv run -p .venv/bin/python ruff check dev-samples/overpass_water_overlay_probe.py tests/test_overpass_water_overlay_probe.py`
- `uv run -p .venv/bin/python pytest tests/test_overpass_water_overlay_probe.py -q`
- `uv run -p .venv/bin/python -m py_compile dev-samples/overpass_water_overlay_probe.py tests/test_overpass_water_overlay_probe.py`

