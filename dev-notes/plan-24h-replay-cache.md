# Implementation Plan: 24-Hour Screen Replay History

> Superseded for replay behavior by `plan-historical-layer-replay.md`. The
> scene-manifest/full-screen sections below record an earlier approach. The
> already implemented retention subsection remains historical implementation
> context; use the newer plan for further replay work.

## Goal

Keep normal cache data for 24 hours in timestamped directories, following GOES/Himawari's existing source layout, so sampled zstarview screens can be recreated. Do not create a second copy of layer data for replay. Keep live-data freshness separate from disk retention. Replay must load the exact files referenced by the scene sample and must not fetch current data as a fallback.

Initial recording cadence: one scene sample per minute, matching the existing periodic debug-frame interval. Aircraft source observations remain at the existing OpenSky refresh cadence; replay recomputes aircraft projection for each sample time from the archived state vectors.

## Phase 1: Temporal cache primitives

- Add shared UTC path, expiry, atomic-write, and cleanup helpers under a Qt-independent cache module.
- Move aircraft and artificial-satellite records into timestamped subdirectories under their existing cache roots. Store both `fetched_at_utc` and source observation/element time in each record.
- Add collision-safe request IDs to data filenames and UUID `sample_id`s to scene manifests. Do not require process `session_id`; scene samples are uniquely identified on their own.
- Keep a versioned scene-manifest reader/writer under `src/zstarview/replay/`. Manifests contain references and view state only, never copies of aircraft, satellite, or cloud data.
- Write cache data and manifests to temporary sibling files, validate and hash them, then atomic rename. Current cache lookup selects the newest matching temporal record, applying existing freshness/fallback rules independently from the 24-hour disk TTL.
- Keep legacy flat aircraft/satellite cache files readable during migration. New successful fetches write temporal entries without destructively rewriting legacy data.
- Expire scene manifests first, then remove layer records only when older than 24 hours and unreferenced. Ignore active temporary files and safely recover old abandoned temporary files.
- Keep normal cache roots; no database or new package dependency.

## Phase 2: Time-index the normal layer caches

### Aircraft

- Extend `src/zstarview/aircraft/cache.py` to append one immutable record in the normal cache after a successful OpenSky fetch. Current cache lookup selects the latest matching bbox record. Preserve the 5-minute refresh and 10-minute stale fallback.
- Record bbox, fetch time, complete state-vector list, schema, source, and digest. Do not create a new historical record for a cache hit or stale fallback; scene manifests can reference the prior successful record.
- Keep data in `AIRCRAFT_CACHE_ROOT_DIR/YYYY/MM/DD/HH/mm/<request>-<bbox>.json` so bounding-box changes remain distinct. No replay-only copy is written.

### Artificial satellites

- Extend `src/zstarview/satellites/cache.py` to append each successfully accepted group result to the normal satellite cache by acquisition time, group, and Horizons observer scope. Current lookup selects the newest matching record. Preserve the existing 24-hour element freshness/fallback behavior.
- Record OMM/TLE records, element epoch, fetch time, provider/source, group/scope, and digest. Archive only the final accepted group result; do not persist Horizons callback fragments as independent usable revisions.
- Keep mutable failure/backoff/rate-limit metadata in a small operational metadata file, separate from immutable time-indexed satellite data. A failed fetch does not create a data revision.
- Keep data in `SATELLITE_CACHE_ROOT_DIR/YYYY/MM/DD/HH/mm/<request>-<group>-<scope>.json`; observation time/group/scope are sufficient to select the layer data for a scene.

### Clouds

- Retain existing time-keyed GOES/Himawari files as the normal cloud cache. At cloud worker completion, record selected source key/time, actual source paths and digests, observer location, and completed grid/render revision in scene/cache metadata; do not duplicate source files for replay.
- Publish a cloud revision only after the worker result manifest, source assets, and derived grid pass validation. Keep the prior revision active until then; a partial worker result never replaces it.
- Update cloud cleanup in `src/zstarview/clouddisc/cache/cleanup.py` to honor source assets referenced by unexpired scene manifests. Keep the 24-hour retention policy for unreferenced assets.

## Phase 3: Record coherent scene samples

- Add a `ReplayRecorder` that samples once per minute and receives an immutable snapshot of the exact `FrameContext` and `RenderInputs` used for that paint. Integrate the request/flush points alongside, but independently from, `ZSTARVIEW_DEBUG_SAVE_PERIODIC_FRAME` in `gui/window.py`, `gui/window_updates.py`, and `gui/window_render.py`.
- Capture only serializable values and layer revision references on the GUI thread. Submit disk serialization, digest checks, and manifest commit to a background service so painting and timers stay responsive.
- Store UUID `sample_id`, UTC sample time, observer/location, timezone, time mode, viewport size, view center/FOV, enabled overlays, relevant style/opacity settings, app/build and renderer schema versions, and each layer's ordinary-cache relative path/digest/status/source time. Scene time and view/location state select layer history; process `session_id` is not required.
- The manifest refers to the last committed revision actually present in that rendered frame. While a fetch is in flight, retain the previous committed revision. Record `unavailable` explicitly when the original frame had no usable layer data.
- The recorder is enabled by default. It records only data from enabled/acquired layers and must not turn on overlays or initiate extra network requests.

## Phase 4: Replay and export

- Add a read-only replay loader that resolves manifests and references, validates hashes/schema versions, and creates scene inputs for an exact recorded sample time.
- Reuse the existing render pipeline and `FrameContext`/`RenderInputs` construction. Run with live controllers and network acquisition disabled. Never substitute a current cache entry when a recorded reference is missing or invalid.
- Preserve current overlay semantics outside replay; do not change the existing rule that arbitrary past GUI time disables live-only overlays.
- Add a separate `zstarview-export-replay` CLI (registered in `pyproject.toml`) accepting a UTC start/end range and output directory. Initially render one PNG per recorded minute; users can encode the sequence with their existing video tool. Keep video encoding and new codec dependencies out of this change.
- Reject unsupported renderer/schema revisions clearly. Missing/corrupt layer inputs produce an explicit unavailable layer and a non-success completeness result rather than silently rendering with today's data.

## Phase 5: Documentation and release behavior

- Update `docs/specification.md` only after the behavior ships: default recording, 24-hour availability, sample cadence, replay/export command, and incomplete-history behavior.
- Update `docs/design.md` from proposal to implemented architecture, including module responsibilities and lifecycle.
- Add a user-facing `release-notes.md` entry with the shipped replay behavior.

## Validation plan (for implementation work)

- Unit coverage for UTC directory construction, unique same-minute paths, schema round-trip, digest rejection, atomic-write recovery, exact 24-hour cutoff, and reference-aware cleanup.
- Cache tests proving aircraft history is written only for network success and live latest-cache TTLs are unchanged; satellite tests proving group/scope isolation and no archival of partial callback fragments.
- Cloud tests proving only completed worker results become revisions and referenced GOES/Himawari files survive cleanup until all scene references expire.
- GUI-level tests proving a scene manifest refers to the `FrameContext`/`RenderInputs` revisions that were rendered, while an in-flight cloud update does not alter an already queued sample.
- Replay tests with synthetic inputs proving deterministic repeated rendering, no network calls, target-time satellite propagation, aircraft projection at sample time, and explicit failure for missing/corrupt references.
- CLI checks for range selection, output naming/order, empty ranges, unsupported schema, and cleanup concurrent with an export.
- Run the focused tests first, then the full test suite and static checks. Compare a generated replay frame against the source frame for identical inputs.

## Risks and boundaries

- Only samples recorded after this feature is enabled can be replayed; existing current-value caches cannot backfill the previous 24 hours.
- Disk-full or write errors must not fail live rendering. Mark that interval unrecorded and log the reason; never delete unexpired replay data to recover space.
- Initial fidelity covers clouds, aircraft, artificial satellites, and recorded view/style settings. Other changing network overlays need their own archived inputs before they can be promised in a full-screen replay.
- A one-minute manifest captures the displayed state at that cadence, not transient camera gestures or every repaint between samples.
- Exact output across application updates is not guaranteed. Store renderer/schema versions and reject incompatible replay rather than implying bit-for-bit equivalence.

## Implementation Plan: 24-Hour Cache Retention

This implementation slice changes the storage lifetime of the three layer inputs. It does not yet add scene manifests, replay UI, or replay export.

### Retention contract

- Add an explicit 24-hour retention duration for aircraft and artificial-satellite source records. Keep aircraft's 5-minute refresh and 10-minute live stale fallback unchanged; those govern whether data is usable for the current display, not how long its file remains available.
- Keep satellite element validity, stale fallback, and retry/backoff behavior unchanged. A source record's retention must not make expired orbital elements appear fresh.
- Keep the cloud source cache's existing 24-hour cleanup. Preserve its existing newest-file safeguard and in-progress-download protection; these protect operational recovery and are independent of replay selection.
- Use each payload's `fetched_at_utc` as the age source where available; do not use directory mtime or data epoch as a substitute. Expire at `age > 24h`, retaining an entry at the exact cutoff.
- Cleanup is best-effort and must not interrupt fetching or rendering. Preserve the aircraft rate-limit/lock files and satellite failure/backoff metadata when pruning data records.

### Aircraft cache migration

- Keep `AIRCRAFT_CACHE_ROOT_DIR/<bbox>.json` as the mutable latest-value file so existing callers and old installations continue to work.
- Keep the flat latest-value file as a compatibility mirror and write every successful acquisition as an immutable record under `history/YYYY/MM/DD/HH/mm/<bbox>-<unique-id>.json`, with bbox, fetch time, source, and snapshots. Before the first overwrite of a legacy flat file, archive that prior valid payload too; avoid duplicating a revision already present in history. Do not archive cache hits, rate-limit skips, or stale fallbacks.
- Make `load_aircraft_cache` compare the legacy/current latest file and matching history records, returning the greatest `fetched_at_utc`. This permits older installations' flat cache to remain readable and allows later replay lookup to use retained records.
- Change cleanup's default data retention from the 10-minute stale-fallback constant to a dedicated 24-hour retention constant. Scan current and history payloads, prune by payload fetch time, leave `rate_limit.json` and `fetch.lock` alone, and remove empty history directories.
- Continue using `stale_fallback_seconds` exclusively for fetch-failure behavior. The cleanup call after a successful fetch uses the retention constant, not this parameter.

### Artificial-satellite cache migration

- Keep `SATELLITE_CACHE_ROOT_DIR[/<scope>]/<group>.json` as the mutable latest-value/backoff file for compatibility with current GUI startup and cache APIs.
- Write every successful group refresh as an immutable record under `history/YYYY/MM/DD/HH/mm/<group>-<unique-id>.json` within the same scope directory. Before the first overwrite of a legacy flat file, archive that prior valid payload too; avoid duplicating a revision already present in history. Preserve `element_epoch_utc`, `fetched_at_utc`, source, group records, and fetch metadata as present at archival time. Archive only accepted complete group results, never Horizons callback fragments or failed requests.
- Keep the current failure behavior and metadata in the latest file; a failed request must not create a new data revision. Historical files stay unchanged after archival. Successful refreshes continue to select/cache the newest current value through the existing flat path; that file is a compatibility mirror and history is the immutable retained source.
- Add cleanup for satellite element data, both current files and history, using payload `fetched_at_utc`, a dedicated 24-hour retention constant, scope isolation, and empty-directory pruning. Keep metadata-only failure records intact; remove a data payload only after its fetched time exceeds 24 hours. Existing data-epoch freshness checks remain unchanged.
- When the app successfully writes a new group value, invoke history cleanup best-effort. Keep operational fetch/backoff behavior separate from the history retention duration.

### Cloud cache and integration boundaries

- Do not move or duplicate GOES/Himawari files. Their existing UTC source paths and 24-hour cleanup already satisfy this storage slice.
- Do not yet change cloud cleanup's reference behavior; manifests do not exist in this slice. Add reference-aware deletion only with scene-manifest implementation, before replay can rely on references.
- No process session ID is added. Existing observer location and bbox/group/scope remain the layer selection context.
- Preserve legacy flat aircraft and satellite files in place during migration. New fetches begin building history; no destructive one-time conversion or network backfill is needed.

### Completion checks

- Inspect the diff for unchanged live freshness constants and separately named 24-hour retention constants.
- Review all save/overwrite and cleanup call sites to confirm only successful source updates create history and failures/hits do not.
- Run `git diff --check`. Tests are intentionally not run in this implementation turn.
