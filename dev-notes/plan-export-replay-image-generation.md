# Plan: Implement cached replay image generation

Date: 2026-10-16
Status: Basic sequence export implemented; see implementation progress below.

## Goal and current baseline

Extend `zstarview-export-replay` from metadata-only cache range queries to a
minute-stepped PNG sequence generated from retained local data. Keep
`--list-cache-ranges` read-only and available. The exporter must not fetch
weather, aircraft, or satellite data, and must not mutate or prune its input
cache. Use the renderer already shared with `zstarview-export-image`, while
building replay scenes from history selected for each frame time.

Current code already provides the CLI entry point, location resolution,
metadata inventory, rounded cache-scope checks, indexed time selection,
layer and common ranges, per-satellite-object reporting, and text-only
relative display of common ranges. It does not load selected layer payloads,
construct a rendering plan, initialize Qt for image export, draw replay
frames, or write image jobs. `export_image.py` has the single-image pipeline,
but its main path also performs live acquisition and should not be called by
replay.

The approved time interval is half-open `[start, end)`, sampled at each UTC
minute. `--start` and `--end` accept zoned ISO 8601 or a relative negative
offset such as `--start=-08:30`. Capture one UTC reference at command entry,
floor it to a minute, and use it to resolve every relative endpoint. A
24-hour interval emits 1440 frames. The same parser and reference-time rule
apply to range-query filters and export mode.

## Implementation stages

### 1. Split query and export modes; freeze time and settings

- Preserve `--list-cache-ranges` as its existing metadata-only mode.
- Add an export mode requiring `--start`, `--end`, and `--output-dir`; reject
  contradictory modes and missing endpoints before resolving a location.
- Parse absolute endpoints and `-HH:MM` relative offsets, allowing either
  form on either endpoint. Require minute boundaries for absolute values,
  resolve all times to UTC, and validate `start < end` after resolution.
- Capture the run reference once before location lookup so search latency does
  not shift a relative interval. Carry the resolved UTC bounds forward.
- Normalize existing view/layer arguments into an immutable settings object.
  Keep existing aliases and defaults; keep unsupported positive-opacity layers
  as errors. Add the already designed `--hud` and `--observation-info`
  validation without changing normal `zstarview` or `export-image` behavior.
- Decide whether the range-query JSON `checked_at_utc` should use the command
  reference captured before scanning. Preserve the current absolute range JSON
  schema; relative common-range labels remain text-only.

Deliverable: parser and pure time/settings helpers, plus mode-specific help and
errors. No Qt initialization, cache payload reads, or output-directory writes
in this stage.

### 2. Build a deterministic replay plan and preflight

- Add a Qt-independent plan model under `zstarview.replay` containing the
  resolved location/scope, settings digest, UTC frame sequence, per-frame
  selected history references, static dependencies, and missing/unsupported
  reasons.
- Reuse `Inventory` and `HistoryIndex` rather than scanning directories for
  every minute. Use the same selectors as range coverage and record a stable
  tie-break for same-time candidates.
- Preflight every requested minute and every enabled layer. Check the actual
  ephemeris time coverage, catalog availability, observer-height resolution,
  and static source availability needed by the selected renderer. Do not treat
  metadata-only `unverified` results as render-ready.
- By default, fail before rendering if any required frame/layer is missing or
  unsupported. With `--allow-partial-data`, retain the plan and mark each
  affected frame/layer explicitly. Never substitute current live data.
- Validate that the output directory is not inside the input cache. Estimate
  disk needs where practical. `--plan-only` reports the plan and missing
  inputs without initializing Qt or creating output files.

Deliverable: deterministic plan JSON in memory and a useful preflight report.
The plan must be independent of wall-clock time after endpoints are resolved.

### 3. Add read-only payload adapters and freeze job inputs

- Implement per-layer adapters for native/Geo-satellite clouds, aircraft,
  ISS, and Horizons. Reuse existing decoders and cache-key rules, but do not
  call fetchers, live-time loaders, cleanup/purge paths, or legacy cache
  writers.
- Verify selected payload digest/schema/type/scope and decode it into the
  existing render-layer data types. Preserve `empty`, `missing`, `corrupt`,
  `incompatible`, and `unsupported` as distinct outcomes.
- Resolve and validate the static inputs actually needed by the enabled
  rendering layers, including terrain/building/water/lighting sources and
  local catalogs/ephemeris. Record source versions/digests in the plan.
- Copy selected dynamic and static dependencies into job-local `inputs/`
  using content digests and atomic rename. Deduplicate identical payloads;
  do not depend on hard links. After freezing, rendering reads only job-local
  inputs. If a selected source disappears during copy, fail or mark partial
  according to mode; do not silently reselect a different record.
- Keep all cache-root access behind the input-root abstraction so tests can
  fail on writes and the exporter cannot accidentally update acquisition
  caches.

Deliverable: a frozen, digest-checked input bundle and plan provenance. No
rendering yet.

### 4. Render one replay frame through the shared pipeline

- Extract only the pure setup/render pieces needed from
  `export_image_inputs.py`, `export_image_layers.py`, and
  `export_image_render.py`. Do not call `export_image.main` or its
  acquisition-capable initialization path.
- Initialize one Qt GUI application, fonts, catalogs, ephemeris, compositor,
  and immutable render style for the job. Render on the main thread and use
  `ApplicationServices.native_work_lock` for Astropy/Skyfield/NumPy native
  operations.
- For the first requested minute, assemble a fresh viewer and
  `RenderSceneData` from the target-time celestial calculation and the
  selected frozen layers. Recompute time-dependent sun/moon lighting and
  projections for that frame. Explicitly set absent layers to empty so no
  previous-frame pixels or state can leak into the output.
- Apply image dimensions from `--image-size`. Default HUD/status/menu/grip and
  window chrome remain absent. Render only requested HUD fields at the
  configured observation-info position, using the frame timestamp.
- Write one PNG to a temporary path and atomically finalize it. Confirm image
  dimensions and record image digest and per-layer status in a frame sidecar.

Deliverable: a one-frame end-to-end export using the real plan and adapters.
This stage is the first point that creates output files.

### 5. Add the minute loop, manifest, and failure accounting

- Iterate UTC frame times in ascending order using `start + n minutes` while
  `frame_time < end`; do not use local-time/DST iteration.
- Rebuild only time-dependent scene data each frame; reuse safe immutable
  catalogs/style/static geometry and payloads. Do not parallelize native
  astronomy/render work in the first version.
- Write `manifest.json`, frame PNGs named with UTC timestamps, and one JSON
  sidecar per frame. Record normalized settings, resolved endpoints and
  relative-time reference, application/renderer version, frozen input
  digests, selection provenance, and `ok/empty/disabled/missing/corrupt/
  incompatible/unsupported/failed` layer states.
- Finalize each artifact atomically. Keep progress and diagnostics on stderr;
  reserve stdout for machine-readable output only if an explicit JSON mode is
  later designed.
- Exit 0 when all frames/layers are complete, 2 when partial output completed
  with missing layers, and 1 on invalid input, failed preflight, or rendering
  failure. Never count a failed frame as complete.

Deliverable: a complete PNG sequence and auditable job manifest for a short
interval.

### 6. Add safe resume and operational validation

- Reject a non-empty output directory unless `--resume` is specified. Bind
  resume to the normalized settings, resolved endpoints, plan digest, input
  bundle, and renderer version.
- Skip a frame only when its sidecar, PNG digest, and plan digest all match.
  Report orphan PNGs, stale sidecars, or digest mismatches without silently
  overwriting them. Resume only with the frozen job-local input bundle; never
  resolve a new location or read newer cache state.
- Run a small real-cache export after fixture-level preflight succeeds. Then
  measure 24-hour runtime, peak memory, number of copied inputs, output size,
  and interruption/resume behavior. Do not promise completion time before
  those measurements.
- Update user-facing specification and release notes only when image export
  is implemented and its supported layer matrix is known. Update design and
  session notes as implementation decisions change.

Deliverable: interrupted jobs can resume without changing their selected
inputs or silently accepting incomplete frames.

## Verification gates

- Time/settings: ISO offsets and UTC conversion, relative negative offsets,
  minute flooring, `-00:00`, mixed absolute/relative endpoints, 24-hour/1440
  frame count, and rejection of empty/reversed intervals.
- Plan/selection: selector and common-range agreement at every minute,
  future-acquisition exclusion, rounded scope matching, same-time tie-break,
  empty aircraft/satellite snapshots, partial Horizons objects, and bounded
  final tails.
- Read-only behavior: monkeypatch all provider fetches, cache writers,
  cleanup, and deletes to fail; assert query/preflight/render performs none.
  Confirm source cache files are unchanged after input freezing and rendering.
- Payload correctness: missing/corrupt/digest mismatch/schema mismatch/type
  mismatch and a source removed while copying; each gets the documented
  status and strict/partial behavior.
- Rendering: Qt offscreen single frame, output dimensions, transparent/opaque
  background behavior, no default HUD/interactive UI, selected HUD items,
  missing-layer residue prevention, and changing celestial/lighting output
  across adjacent minutes while a cloud source is unchanged.
- Job integrity: atomic writes, manifest/sidecar digest matching, strict
  preflight with no output, partial exit code, failed-frame accounting,
  resume validation, and refusal to overwrite mismatched/orphan files.
- Run focused tests at each stage, then relevant existing export-image and
  replay regressions; run Ruff, mypy for new modules, and diff checks.

## Scope boundaries

This plan excludes video encoding, real-time capture, reconstruction of GUI
uptime or failed-request intervals, new provider acquisition, frame-parallel
rendering, window geometry/interaction replay, and changing normal GUI cache
selection or rendering behavior.

## Implementation progress (2026-10-16)

- Implemented list/export mode separation, UTC minute endpoints, relative `-HH:MM` endpoints,
  `--plan-only`, strict preflight, partial output, image size, selected HUD, progress, manifests,
  per-frame sidecars, content-digest input copies, and guarded resume.
- Implemented native cloud, aircraft, ISS, and Horizons selection/decoding through the existing
  history index and render pipeline. Sky and celestial data are recomputed at each target minute.
- Verified a real-cache one-minute plan and PNG export. A 640x480 frame included cloud, aircraft,
  ISS, and Horizons; PNG dimensions and sidecar/manifest SHA-256 agree. A 320x240 voxel-cloud frame
  also completed after fixing lock recursion.
- Current unsupported inputs: Geo-satellite cloud rendering, water, road lights, diffuse sky/AKARI,
  meteor trails, and tropical cyclones. Terrain horizon, prepared urban outlines, Night Lights,
  and ridge glow use zstarview's CLI opacity defaults and read their local caches only.
- A local-cache real frame reports terrain, urban, and Night Lights as `ok`; with no local DEM it
  uses the flat-ground horizon fallback. Large Night Lights tiles are fingerprinted in place rather
  than copied into every job.
- Resume reuses the saved resolved location and checks settings/input plan digests, but currently
  requires the source cache records to remain readable. Full job-bundle-only resume remains follow-up.
- Validation run: compileall, Ruff, targeted mypy, real-cache plan-only, and single-frame exports.
  No test suite or 24-hour performance measurement was run in this turn.
