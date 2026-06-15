# 2026-05-26 Handoff

## Scope

Geo-satellite GUI integration.

## Current State

- Geo-satellite support is already wired into `zstarview-export-image`.
- The export path now treats `--layer-timeout-seconds` as download/network fetch time only.
- Post-download processing such as Geo-satellite projection is no longer counted toward timeout.
- Geo-satellite raw cache is latest-only, and intermediate proxy/inpaint work is done in memory.
- The Geo-satellite pipeline lives under `src/zstarview/geosatellite/`.
- Runtime assets for the experimental path live under `src/zstarview/data/geosatellite/`.

## Relevant Files

- `src/zstarview/cli/export_image.py`
- `src/zstarview/geosatellite/`
- `src/zstarview/data/geosatellite/`
- `src/zstarview/gui/cloud_controller.py`
- `src/zstarview/gui/window_inputs.py`
- `src/zstarview/gui/window_updates.py`
- `src/zstarview/gui/viewer.py`

## What Still Needs To Be Done For GUI

1. Decide where the Geo-satellite toggle should surface in the GUI, if at all.
2. Reuse the existing Geo-satellite pipeline from the GUI cloud path instead of duplicating logic.
3. Add GUI status text that matches the export path as closely as practical:
   - `Geo-sat + Downloading`
   - `Geo-sat + Projecting`
   - `Geo-sat + <frame time or error>`
4. Keep the GUI behavior consistent with the rest of the app:
   - download failures should surface as an unavailable cloud state
   - the GUI should keep running even if Geo-satellite fails

## Key Decisions Already Made

- Geo-satellite is an experimental Europe-only path.
- It does not blend with GOES/Himawari when enabled inside its responsibility region.
- `--geo-satellite` is a `true|false` option in the batch/export path.
- The export path is batch-oriented and should abort on required-layer download failures unless partial output is explicitly allowed.
- Timeouts are intended to cover download/network fetch only.

## Validation Already Run

- `uv run -p .venv/bin/python ruff check src/zstarview/cli/export_image.py tests/test_export_image_layer_gating.py docs/specification.md docs/design.md`
- `uv run -p .venv/bin/python pytest tests/test_export_image_layer_gating.py tests/test_export_image_sixel.py tests/test_geosatellite_helpers.py -q`
- `uv run -p .venv/bin/python pytest -q`

## Notes

- The latest timeout fix is committed as `40f20853` (`fix: limit export-image timeouts to downloads`).
- No GUI timeout changes were needed for this fix.
- If Geo-satellite is added to the GUI next, the main work is wiring the existing pipeline into the GUI cloud controller and aligning status text, not changing the underlying download/cache logic.
