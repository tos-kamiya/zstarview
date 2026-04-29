# Session 2026-04-30

- Topic: Keep export-image rendering when cloud data is unavailable
  - Decision: Treat cloud-layer download/render failure as non-fatal and continue with the remaining scene layers.
  - Rationale: A missing cloud source should not suppress the entire exported image; the sky, labels, and other overlays still provide useful output.
  - Notes: Preserve the existing fatal behavior for terrain, urban, aircraft, and satellite layer failures unless `--allow-partial-data` is set.
  - Validation: `uv run -p .venv/bin/python ruff check src/zstarview/cli/export_image.py tests/test_export_image_sixel.py`
  - Validation: `uv run -p .venv/bin/python pytest tests/test_export_image_sixel.py -q`

- Topic: Clarify the partial-data abort message
  - Decision: Make the fatal partial-data path explicitly mention `--allow-partial-data` and explain that it still allows image output when non-cloud data cannot be downloaded.
  - Rationale: The exit path was warning-only before the abort, which made it harder to understand that the CLI could be retried in a partial-data mode.

- Topic: Bump package version
  - Decision: Bumped `zstarview` from `1.21.5` to `1.21.6`.
  - Rationale: This records the latest CLI behavior fixes and the export-image partial-data guidance update as a patch release.
