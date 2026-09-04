# Rendering Dataflow Refactor Handoff — 2026-09-04

## Goal

Continue the original rendering-pipeline cleanup by removing:

- arguments that are unused or no longer necessary;
- the same data passed through multiple paths;
- derived values passed alongside their source data when the receiver can derive
  them cheaply and unambiguously.

The state-ownership work that needed to precede this cleanup is complete.

## Completed work

The following commits form the baseline for the next session:

- `d8dc694a refactor: remove global rendering state`
  - Made the molecular-cloud source explicit.
  - Moved the star render cache under `SkyCompositorCache` ownership.
- `6e6ad7fc refactor: add application-owned services`
  - Added `ApplicationServices` with an executor, native-work lock, and lazy
    ephemeris provider.
  - Made `SkyWindow` own and inject the services instance into the core update
    path.
- `7e3e3828 refactor: complete GUI service ownership`
  - Routed the remaining GUI controllers and dialogs through
    `ApplicationServices`.
  - Removed `gui/worker_pool.py` and `gui/native_work_lock.py`.

The user manually confirmed application behavior before the last two commits.
After the latest implementation change, the full test suite reported 1660
passed with two existing cloud-grid cast warnings. GUI Ruff checks, compileall,
and `git diff --check` also passed.

## Decisions to preserve

- GUI background resources have one application-level owner. Do not recreate
  process-global GUI pools or locks.
- `zstarview-export-image` is intentionally treated as a one-shot process. Its
  calls to the process-lifetime `astro.load_ephemeris()` cache may remain; do not
  inject GUI `ApplicationServices` into the export CLI merely to remove that
  cache.
- Pure, bounded memoization caches are not automatically considered problematic
  global state. Focus on mutable runtime state and ambiguous data ownership.
- A derived argument should only be removed when derivation is cheap,
  deterministic, and based on data already available to the receiver. Preserve
  expensive computed snapshots and values that must be consistent for the
  entire frame.

## Current rendering dataflow findings

### 1. `FrameContext` and `RenderSceneData` duplicate frame inputs

**Status: complete (2026-09-04).** The duplicated `time_obj` and `viewer` paths
have been removed from `RenderSceneData`. `FrameContext.time_obj` and
`FrameContext.viewer` are now the canonical frame values throughout the GUI,
export renderer, top-level pipelines, and shared rendering helpers. All
production helper signatures require an explicit `ViewerData`; no temporary
direct-caller compatibility fallbacks remain. Tests that need a convenient
scene viewer use a test-only fixture wrapper.

`src/zstarview/render/render_types.py` now defines:

- `FrameContext.viewer` and `FrameContext.time_obj` as canonical frame inputs;
- `RenderSceneData` as calculated scene/layer data only;
- `FrameContext.geometry`, while downstream code may reconstruct projection or
  viewport information from viewer fields and painter state.

Choose one canonical path for viewer/time data before shrinking downstream
function signatures. The first implementation stage should be limited to this
boundary and should retain behavior.

Candidate direction: make `FrameContext` the frame/view/projection context and
make `RenderSceneData` contain calculated scene/layer data only. Confirm all GUI
and export construction sites before committing to this direction, because the
export renderer also builds these objects.

### 2. `SkyCompositorCache.draw()` receives an exploded viewer/projection

**Status: complete (2026-09-04, commit `8c98c719`).** The compositor now
receives one canonical `ViewerData` value. It derives view center, observer
location/height, and FOV values from that object, and no longer reconstructs a
partial viewer from parallel scalar arguments. `ScreenGeometry` remains a
separate calculated per-frame value.

Before Sec. 2, `src/zstarview/gui/composite.py` accepted separate values
including:

- `view_center`;
- `observer_lat_deg`, `observer_lon_deg`, and `observer_height_m`;
- `edge_fov_deg` and `content_fov_deg`;
- `geometry`, while also reading the painter viewport.

It then reconstructs a partial `ViewerData` named `glow_viewer_data`. This is a
strong sign that the source object or a smaller purpose-built context should be
passed instead of parallel scalar arguments. Avoid making the compositor depend
on all application state; use the narrowest existing value object that expresses
the required invariant.

### 3. Solar altitude is passed twice

**Status: complete (2026-09-04, commit `8c98c719`).** `sun_altaz` is now the
single solar input to `SkyCompositorCache.draw()`. The compositor derives the
solar altitude locally for night-light, cloud, and cache-key consumers. The
optional-`None` behavior is preserved when no Sun is present.

Before Sec. 3, `src/zstarview/render/zstarview_pipeline.py` derived both:

- `sun_altaz` from `scene.celestial_data`;
- `night_light_sun_alt_deg` from the same `scene.celestial_data`.

Both are passed to `SkyCompositorCache.draw()`. Since the latter is
`sun_altaz[0]`, prefer a single canonical solar value and derive altitude at the
use site. Check cache-key construction and all optional-`None` behavior when
making this change.

### 4. Projection values can be grouped

**Status: in progress (2026-09-04).** The `draw_cloud_overlay()`, dim-alt ring,
terrain-clip, common inverse-projection, molecular-cloud, and cloud-render
inverse-projection helper boundaries now accept `ViewProjection` rather than
separate view center/FOV values. The remaining projection helpers still need
an audit to distinguish useful grouping from APIs that intentionally accept
only a subset of projection data.

`src/zstarview/types.py` already contains `ViewProjection` with `view_center`,
`edge_fov_deg`, and `content_fov_deg`. Audit whether downstream projection calls
can consume this object rather than receiving its fields separately. Do not
blindly replace `ScreenGeometry`: it represents calculated screen dimensions and
may correctly remain a per-frame computed value.

### 5. Continue the unused-argument audit

Once the canonical context boundary is established, use call-site searches and
tests to remove parameters that are accepted but not read. Change one logical
API boundary at a time so failures identify the affected layer.

## Recommended implementation sequence

1. Enumerate every construction and consumption site for `FrameContext` and
   `RenderSceneData`, including GUI, tests, and `cli/export_image_render.py`.
2. Add or adjust focused tests that prove GUI and export render entrypoints use
   the same viewer/time values.
3. Remove the duplicated viewer/time path in one small commit. **The
   `time_obj` portion is complete; the `viewer` portion remains.**
4. **Sec. 2 complete:** Replace the exploded projection/observer arguments to
   `SkyCompositorCache.draw()` with an appropriate context or value object.
5. **Sec. 3 complete:** Collapse `night_light_sun_alt_deg` into `sun_altaz`
   and update compositor cache keys and tests.
6. **Sec. 4 in progress:** Continue the `ViewProjection` consumption audit;
   cloud overlay, dim-alt ring, terrain-clip, inverse-projection,
   molecular-cloud, and cloud-render inverse-projection helpers are migrated.
7. **In progress:** Run a final unused-parameter search through `render/` and
   `gui/composite.py`; the unused `time_obj` in the diffuse-sky opacity helper
   and instrument context-layer helper has been removed, while the remaining
   findings require boundary-by-boundary review.

For each stage, inspect both GUI and export-image callers. Keep commits focused;
do not combine visual behavior changes with signature cleanup.

## Validation commands

```text
uv run -p .venv/bin/python ruff check src/zstarview tests
uv run -p .venv/bin/python mypy --install-types --non-interactive src/zstarview tests
QT_QPA_PLATFORM=offscreen .venv/bin/python -m pytest -q
uv run -p .venv/bin/python -m compileall src
git diff --check
```

Run focused rendering and export tests during each step, then the full suite
before asking for manual confirmation or committing.

## Worktree caution

At handoff time, the refactor commits are complete, but the worktree contains
unrelated user-owned changes and untracked generated files under `docs/images/`,
`docs/images/scheduled/`, `build/`, and an untracked
`dev-notes/session-2026-08-20.md`. Do not stage, overwrite, or delete them. The
handoff file itself is newly created and uncommitted.

## Suggested opening task for the next session

Start with a read-only call graph of `FrameContext` and `RenderSceneData`. Propose
the exact canonical owner of `viewer`, `time_obj`, geometry, and projection before
editing. Then implement only the first duplication removal with focused tests.
