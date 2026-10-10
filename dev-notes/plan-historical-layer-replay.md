# Implementation Plan: Historical Layer Display from 24-Hour Caches

## Goal

Allow the existing `--datetime` and negative `--hours` launch-time selection to
display cached cloud, aircraft, and artificial-satellite layers for targets up
to 24 hours in the past. Use the normal render pipeline and the current view,
brightness, and style options. Historical display must not fetch network data or
substitute a current cache. When `--place` is used, persist the raw query and
resolved observer location as a timestamped launch metadata record.

This is temporal layer display, not complete screen capture or restoration. Do
not add per-frame scene manifests, screenshot recording, replay export, or
process session IDs.

## Time-selection contract

- A target within the existing present tolerance keeps current live behavior.
- A past target no more than 24 hours old enables the three supported cached
  layers. Other time-sensitive overlays, including precipitation and tropical
  cyclones, remain unavailable in historical mode.
- Historical workers only read local cache records. They must never call a
  source fetcher, HTTP client, S3 listing, or download path.
- A future target, a target older than 24 hours, or a missing/incompatible
  record leaves the affected layer unavailable. Current data is never used as
  fallback.
- Record acquisition time and source-valid time independently. Historical
  selection requires both the correct scope and data that was already acquired
  by the target time.

## Implementation sequence

### 1. Launch location metadata

- Add a Qt-independent `zstarview.replay` module to write and prune launch
  metadata under `CACHE_PATH/replay/v1/launches/YYYY/MM/DD/HH/mm/<uuid>.json`.
- After `resolve_launch_location` succeeds in `_StartupBootstrap._run`, write a
  record only when `args.place` is non-empty. Include schema version, UUID,
  `started_at_utc`, original parsed query, `--place-countrycode`, effective
  `--place-lang`, resolved display name, latitude/longitude, timezone, source,
  ground elevation, and effective observer height.
- Serialize ASCII-safe JSON to a sibling temporary file, then `replace()` it
  atomically. A write or cleanup error is logged and must not fail startup.
- Prune launch records whose recorded UTC start is more than 24 hours old.

### 2. Common historical-target policy

- Extend the GUI overlay-availability calculation with an explicit opt-in for
  cached historical layers, leaving image export and other existing callers on
  their current policy.
- Enable only clouds, aircraft, and satellites for past targets inside the
  24-hour window. Keep present refresh behavior and all future-time behavior
  unchanged.
- Pass the resolved target UTC to each layer controller. Historical controller
  paths must be selected before any fetch path is entered.

### 3. Aircraft history lookup

- Add a reader in `aircraft/cache.py` that searches the existing immutable
  `history/YYYY/MM/DD/HH/mm` records for the requested bbox and chooses the
  latest complete snapshot whose `fetched_at_utc <= target_time_utc`.
- Reject records with source timestamps after the requested time or snapshots
  whose `last_contact_unix` is more than the existing 10-minute live stale
  window behind the requested time.
- Add a historical-only controller branch that emits the selected snapshots
  with source `cache-history`; it must not invoke its configured fetcher.
- Continue projecting snapshots through the existing target-time renderer.

### 4. Artificial-satellite history lookup

- Add a cache reader that searches the appropriate group and, for Horizons,
  observer-scope history directory. Consider history and the compatibility
  latest file, but accept only complete records fetched by the target time.
- For orbit-element groups, select the newest applicable element set with an
  epoch not after the target and within the existing group validity window.
- For Horizons records, require the same observer scope and an applicable
  recorded epoch. Missing matches return `Satellites: unavailable` without
  calling the network fetcher.
- Update `resolve_satellite_elements_for_time` to use that reader for eligible
  historical requests; preserve the present-time code path byte-for-byte in
  behavior.

### 5. Cloud source history and offline loading

- Add a normal, timestamped cloud-source cache under the existing CloudDisc
  cache root. Store the already normalized `CloudSourceData` object as a
  pickle artifact plus an ASCII JSON sidecar containing schema, fetched time,
  source time/key, observer coordinates, source paths, artifact relative path,
  and SHA-256. The existing source assets remain the underlying downloaded
  inputs; the artifact makes the normalized completed source reusable without
  S3/network activity.
- Write the artifact and sidecar via temporary sibling files; publish the
  sidecar last so incomplete pairs cannot be selected. Cache I/O errors remain
  non-fatal to live cloud rendering.
- Add lookup by exact observer location (the same 4-decimal precision as the
  cloud grid key), with both `fetched_at_utc <= target_time_utc` and
  `source_time_utc <= target_time_utc`. Select the newest eligible source no
  older than 30 minutes relative to target; verify artifact path containment,
  SHA-256, schema, type, observer scope, and timestamps before unpickling/use.
- In the cloud controller, use this cache-only lookup for past targets and
  bypass `run_cloud_source_worker_process`. Save successful live source results
  to history once, then continue through the normal camera-dependent render
  path.
- Run cloud-source-history cleanup alongside existing cache cleanup and prune
  by the sidecar's `fetched_at_utc` after 24 hours. Preserve current raw source
  cleanup and in-progress-download behavior.

### 6. Documentation

- Update `docs/design/cache-history-replay.md` from proposal to implemented
  flow/schema, including the existing CLI controls and exact cache selection
  rules.
- Update `docs/design.md` summary and `docs/specification.md` visible behavior:
  supported time range, affected overlays, offline-only history lookup,
  current display settings, and unavailable behavior.
- Append implementation outcome and validation to the existing dated session
  log without removing its existing uncommitted entries.

## Review and validation

- Review call paths to confirm a historical miss cannot reach a fetcher or
  network operation.
- Review 24-hour retention and existing aircraft/satellite live freshness
  constants remain independent.
- Run `git diff --check` and inspect the final status/diff. Do not add or run
  tests in this implementation turn unless requested.

## Known trade-offs

- Cloud source artifacts can be larger than metadata-only records. They are
  retained for at most 24 hours and are needed to render a different view
  center from the same completed source without network access.
- The feature cannot backfill the 24-hour period before cloud-source history
  writing is introduced.
- Historical accuracy is limited by aircraft polling cadence, satellite
  element epochs, and the 30-minute cloud-source selection window.

## Implementation status

Implemented the launch metadata record, cache-only historical selectors for
clouds, aircraft, and artificial satellites, GUI availability wiring, and the
specification/design updates. `git diff --check` and AST parsing of the changed
Python modules passed. Tests were not run per the repository instruction.
