# Implementation Plan — Cloud Alt/Az Grid Rendering

Status: draft  
Target: replace per-frame lat/lon → screen cloud projection with a cached observer-centric (altitude, azimuth) grid.

---

## 0. Goals and Non-Goals

### Goals
- Build `CloudAltAzGrid`, an observer-centric intermediate representation of cloud coverage.
- Cache it by (observer location, satellite, timeslot, cloud shells) so camera rotation/zoom only re-runs the cheap render step.
- Render an MVP cloud image by drawing white circles whose radius/opacity depend on cloud amount.
- Integrate the new pipeline with `CloudController` without breaking the existing `CloudDisc.render_from_source_with_coverage` path.
- Keep CLI `zstarview-export-image` working; the export path can use the same alt/az grid internally.

### Non-Goals (deferred to later phases)
- Removing the old per-frame re-projection path.
- Implementing stripe/hatch rendering from the grid.
- Geo-satellite mode migration (can stay on its existing path initially).
- Missing-data tint implementation can initially fall back to the existing mask path.

---

## 1. High-Level Phases

```
Phase 0 ──▶ Survey & scaffolding
            │
Phase 1 ──▶ Core data + ingestion        CloudAltAzGrid, build_altaz_grid, cache I/O
            │
Phase 2 ──▶ MVP render                     altaz → screen projection, circle drawing
            │
Phase 3 ──▶ GUI integration                CloudController emits grid, view changes re-render
            │
Phase 4 ──▶ CLI/export-image integration   zstarview-export-image uses the grid
            │
Phase 5 ──▶ Polish                         stripes, performance, remove legacy path
```

---

## 2. Phase 0 — Survey & Scaffolding (1 session)

### 2.1 Read existing code
- `src/zstarview/clouddisc/core.py` — `CloudDisc.render_from_source_with_coverage`, shell blending, threshold estimation.
- `src/zstarview/clouddisc/render/grayscale.py` — `_bt_to_weight`, `_suppress_low_cloud_weight`.
- `src/zstarview/clouddisc/projectors/az.py` — ECEF/ENU/altaz math; this will be reused heavily.
- `src/zstarview/clouddisc/workers/cloud_source_worker.py` — where source fetch runs; decide whether to build the grid here or in a follow-up worker.
- `src/zstarview/clouddisc/types.py` — `CloudSourceData`, `CloudMeta`, `SourceKey`, `RenderKey`.
- `src/zstarview/gui/cloud_controller.py` — orchestration signals.
- `src/zstarview/gui/cloud_state.py` — `CloudImageState`.
- `src/zstarview/gui/composite.py` — `compose_cloud_over_sky`, existing stripe renderers.
- `src/zstarview/cli/export_image.py` — `_fetch_cloud_layer` and render flow.

### 2.2 Decide engineering constants
Record the decisions in code constants (e.g. `src/zstarview/clouddisc/altaz_constants.py`):

| Constant | Proposed default | Notes |
|---|---|---|
| `ALT_AZ_GRID_AZ_BINS` | 720 | 0.5° azimuth resolution |
| `ALT_AZ_GRID_ALT_BINS` | 90 | 1.0° altitude resolution (0–90°) |
| `ALT_AZ_GRID_AZ_MIN` | 0.0 | inclusive |
| `ALT_AZ_GRID_AZ_MAX` | 360.0 | exclusive of 360.0, wraps to 0 |
| `ALT_AZ_GRID_ALT_MIN` | 0.0 | horizon |
| `ALT_AZ_GRID_ALT_MAX` | 90.0 | zenith |
| `ALT_AZ_GEO_SAMPLE_STEP_DEG` | 0.2 | lat/lon sampling step during ingestion |
| `ALT_AZ_CIRCLE_BASE_RADIUS_PX` | 2.0 | minimum circle radius |
| `ALT_AZ_CIRCLE_MAX_RADIUS_PX` | 12.0 | maximum circle radius at amount=1 |
| `ALT_AZ_CIRCLE_OPACITY_SCALE` | 0.6 | global alpha multiplier |

These can be tuned during Phase 2 validation.

### 2.3 Add module skeletons
Create empty modules with docstrings:
- `src/zstarview/clouddisc/altaz_grid.py`
- `src/zstarview/clouddisc/altaz_projection.py`
- `src/zstarview/clouddisc/altaz_render.py`
- `tests/test_clouddisc_altaz_grid.py`
- `tests/test_clouddisc_altaz_projection.py`
- `tests/test_clouddisc_altaz_render.py`

### 2.4 Validation
- `ruff check src/zstarview/clouddisc tests` passes with only expected issues.
- No functional changes yet.

---

## 3. Phase 1 — Core Data + Ingestion (2–3 sessions)

### 3.1 `CloudAltAzGrid` dataclass

In `src/zstarview/clouddisc/altaz_grid.py`:

```python
@dataclass(frozen=True)
class CloudAltAzGrid:
    amount: np.ndarray          # (alt_bins, az_bins), float32, 0..1
    missing_mask: np.ndarray    # (alt_bins, az_bins), uint8, 0/255
    alt_min_deg: float
    alt_max_deg: float
    az_min_deg: float
    az_max_deg: float
    observer_lat: float
    observer_lon: float
    satellite: str
    product: str
    time_utc: datetime
    shells_km: tuple[float, ...]
    source_key: SourceKey
    coverage_ratio: float
    source_completeness_ratio: float | None = None
    grid_resolution_deg: float = 0.5
```

### 3.2 Geographic → alt/az conversion

In `src/zstarview/clouddisc/altaz_projection.py`:

```python
def geodetic_to_altaz(
    lat_deg: float,
    lon_deg: float,
    height_km: float,
    observer_lat_deg: float,
    observer_lon_deg: float,
) -> tuple[float, float]:
    """Convert a point on a cloud shell to (altitude, azimuth) for the observer."""


def latlon_grid_to_altaz_grid(
    source: CloudSourceData,
    observer_lat: float,
    observer_lon: float,
    shells_km: Sequence[float],
    alt_bins: int,
    az_bins: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample satellite BT at geographic grid points and accumulate into alt/az bins."""
```

Reuse `enu_basis`, `geodetic_to_ecef`, and ENU→altaz math from `clouddisc/projectors/az.py`.

Key math:
- `P_ecef = geodetic_to_ecef(lat, lon, r_km=shell_km)`
- `v = P_ecef - observer_ecef`
- `enu = [dot(v, east), dot(v, north), dot(v, up)]`
- `alt = degrees(asin(enu[2] / norm(v)))`
- `az = degrees(atan2(enu[0], enu[1]))` normalized to [0, 360)

### 3.3 BT → cloud amount + shell blending

Implement in `altaz_grid.py`:

```python
def _estimate_thresholds_from_source(
    source: CloudSourceData,
    observer_lon: float,
) -> tuple[float, float]:
    """Reuse estimate_bt_warm_hybrid / estimate_bt_cold_hybrid."""


def build_altaz_grid(
    source: CloudSourceData,
    lat: float,
    lon: float,
    *,
    shells_km: Sequence[float] = CLOUD_SHELLS_KM,
    az_bins: int = ALT_AZ_GRID_AZ_BINS,
    alt_bins: int = ALT_AZ_GRID_ALT_BINS,
    geo_sample_step_deg: float = ALT_AZ_GEO_SAMPLE_STEP_DEG,
) -> CloudAltAzGrid:
    """Build a camera-independent alt/az cloud grid from satellite data."""
```

Steps inside `build_altaz_grid`:
1. Estimate `bt_warm`, `bt_cold` using existing hybrid estimator.
2. Estimate scene cloud amount for shell weight blending.
3. Compute blended weights for the 3 shells.
4. For each shell:
   - Build a sparse lat/lon sample grid covering the shell region visible from the observer (roughly above horizon within a generous FOV).
   - Sample BT using `source.sampler` or `build_bt_sampler`.
   - Convert each sample to (alt, az) and to cloud amount.
   - Accumulate into `amount` array using `max` (MVP) or mean with pre-blur.
   - Track which alt/az bins received no valid data → `missing_mask`.
5. Apply light smoothing if blockiness appears (optional).
6. Return `CloudAltAzGrid`.

### 3.4 Cache I/O

In `src/zstarview/clouddisc/altaz_grid.py`:

```python
def altaz_grid_cache_key(
    observer_lat: float,
    observer_lon: float,
    source: CloudSourceData,
    shells_km: Sequence[float],
    grid_resolution_deg: float,
) -> str:
    """Hash key for on-disk cache."""


def save_altaz_grid(
    grid: CloudAltAzGrid,
    cache_root: Path,
) -> Path:


def load_altaz_grid(
    cache_root: Path,
    key: str,
) -> CloudAltAzGrid | None:
```

Cache location: `cache_root / "altaz_grids" / <key_hash>.npz`.

Store:
- `amount`
- `missing_mask`
- JSON metadata

### 3.5 Tests
- `test_altaz_grid_builds_with_empty_source` — all zero, all missing.
- `test_altaz_grid_known_observer_and_satellite_pixel` — place a cold pixel at a known lat/lon/height and verify the expected alt/az cell has high amount.
- `test_altaz_grid_cache_roundtrip` — save/load preserves data.
- `test_altaz_grid_shell_blend_weights` — verify weights change with cloud amount.

### 3.6 Validation
- `pytest tests/test_clouddisc_altaz_grid.py` passes.
- `mypy src/zstarview/clouddisc/altaz_grid.py src/zstarview/clouddisc/altaz_projection.py` passes.

---

## 4. Phase 2 — MVP Render: Circles (2 sessions)

### 4.1 Alt/az → screen projection

In `src/zstarview/clouddisc/altaz_projection.py`:

```python
def altaz_to_screen_coords(
    alt_deg: np.ndarray,
    az_deg: np.ndarray,
    *,
    width: int,
    height: int,
    center_alt_deg: float,
    center_az_deg: float,
    edge_fov_deg: float,
    radius_px: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Project alt/az arrays to pixel coordinates using azimuthal-equidistant model.

    Matches the math in projectors/az.py build_projection_context.
    """
```

Implementation notes:
- Compute angular distance `rho` and position angle `psi` from the center direction.
- Scale `rho` to pixel radius using `radius_px / edge_fov_deg`.
- Handle azimuth wrap at 0/360°: shift input azimuth by ±360° when the angular distance to center is smaller.
- Return `np.nan` for points outside `mask_fov` or below horizon.

### 4.2 Circle renderer

In `src/zstarview/clouddisc/altaz_render.py`:

```python
def render_altaz_grid_circles(
    grid: CloudAltAzGrid,
    width: int,
    height: int,
    *,
    center_alt_deg: float,
    center_az_deg: float,
    edge_fov_deg: float,
    mask_fov_deg: float = 90.0,
    base_radius_px: float = ALT_AZ_CIRCLE_BASE_RADIUS_PX,
    max_radius_px: float = ALT_AZ_CIRCLE_MAX_RADIUS_PX,
    opacity_scale: float = ALT_AZ_CIRCLE_OPACITY_SCALE,
) -> np.ndarray:
    """Render white circles from alt/az grid. Returns (H, W, 4) uint8 RGBA."""
```

Implementation options:
1. **Pure NumPy stencil (MVP):**
   - Iterate over cells with `amount > threshold`.
   - For each cell, compute center pixel, radius, alpha.
   - Use a precomputed circle mask/stencil and add alpha to output.
2. **Gaussian splat (preferred if fast enough):**
   - For each cell, draw a radial gaussian blob: `alpha = amount * scale * exp(-d^2 / sigma^2)`.
   - Accumulate in a float buffer, then convert to uint8.
   - Handles overlap naturally.

Recommended: start with the gaussian splat because it matches the "night light" aesthetic. If too slow, fall back to the stencil approach.

### 4.3 Missing-mask render helper

```python
def render_altaz_missing_mask(
    grid: CloudAltAzGrid,
    width: int,
    height: int,
    *,
    center_alt_deg: float,
    center_az_deg: float,
    edge_fov_deg: float,
    mask_fov_deg: float = 90.0,
) -> np.ndarray:
    """Project the alt/az missing_mask to screen-space uint8 alpha."""
```

For MVP, nearest-neighbor resampling of the grid cells is acceptable.

### 4.4 Tests
- `test_altaz_to_screen_known_directions` — zenith, north, east, horizon map to expected screen quadrants.
- `test_altaz_circle_renderer_returns_shape` — output shape matches `(height, width, 4)`.
- `test_altaz_circle_renderer_preserves_amount` — higher amount cells produce brighter/larger output.
- `test_altaz_missing_mask_wraps_at_azimuth_zero` — cells near 0° are correctly projected when center azimuth is near 360°.

### 4.5 Validation
- Visual sanity check: run a tiny script that renders a synthetic grid (one bright cell at zenith, one at horizon) and inspect output values.
- `pytest tests/test_clouddisc_altaz_render.py` passes.

---

## 5. Phase 3 — GUI Integration (3–4 sessions)

### 5.1 Extend `CloudSourceData` if needed

Currently `CloudSourceData` carries `sampler`, `data_array`, etc. The alt/az grid can be built from these. We may want to add:
- `CloudSourceData.altaz_grid: CloudAltAzGrid | None = None` (optional, computed lazily).

### 5.2 Build grid in the cloud source worker

In `src/zstarview/clouddisc/workers/cloud_source_worker.py`:

Option A: build inside `fetch_cloud_source` after data is available.
Option B: return source only, and build grid in a separate worker step.

**Recommended: Option A** — the worker already runs the expensive download/decode, and the grid build is a one-time cost that should not block the GUI thread.

Update the worker to optionally build the grid when `cfg.use_altaz_grid` is true (new flag, default false during rollout).

### 5.3 Update `CloudController`

In `src/zstarview/gui/cloud_controller.py`:

- Add `CloudController._use_altaz_grid: bool` (default false initially).
- When source is ready, if grid is present in `source`, emit it in `cloud_source_ready` signal.
- Add a new signal or reuse `cloud_ready` to emit the grid.
- On view changes, if grid exists, skip `render_from_source_with_coverage` and instead call `render_altaz_grid_circles` + `render_altaz_missing_mask`.

### 5.4 Update `CloudImageState`

In `src/zstarview/gui/cloud_state.py`:

- Add `altaz_grid` field.
- `set_result()` accepts an optional `altaz_grid`.
- If `altaz_grid` is set, `image` is derived from it on demand rather than stored as the primary source of truth.

### 5.5 Update `SkyWindow` render path

In `src/zstarview/gui/window_render.py` (or wherever cloud overlay is composited):

- If `CloudImageState.altaz_grid` exists:
  - Render cloud image from grid using current `view_center_alt/az`, `edge_fov_deg`, screen size.
  - Composite with `compose_cloud_over_sky` as before.
- Else: keep using `CloudImageState.image` (legacy path).

### 5.6 Feature flag

Introduce a runtime flag to opt into the new path:
- CLI: `--cloud-altaz-grid true` (hidden or experimental initially).
- Config: `cloud_altaz_grid` boolean.
- Default: `false` until Phase 5.

### 5.7 Tests
- `test_cloud_controller_emits_grid_when_enabled` — mock source with grid, verify signal payload.
- `test_cloud_image_state_keeps_grid` — `CloudImageState.altaz_grid` round-trip.
- `test_window_render_uses_grid_when_present` — verify legacy vs. grid path selection.

### 5.8 Validation
- Run GUI with `--cloud-altaz-grid true` and a known location.
- Verify:
  - Cloud appears after download.
  - Rotating/zooming the view re-renders quickly without re-downloading.
  - Missing-data tint appears in correct screen regions.
  - Disabling the flag falls back to the old path.

---

## 6. Phase 4 — CLI / Export-Image Integration (2 sessions)

### 6.1 Update `zstarview-export-image` cloud fetch

In `src/zstarview/cli/export_image.py`:

- In `_fetch_cloud_layer`, when alt/az grid mode is enabled, build or load the grid after source fetch.
- Render the final cloud image from the grid using the export's view center and FOV.
- Keep the existing `render_gray_image_to_cloud_rgba` geo-satellite path untouched if possible.

### 6.2 Cache reuse in export

- If a cached `CloudAltAzGrid` exists for the same observer + source, load it instead of rebuilding.
- This makes repeated exports of the same location much faster.

### 6.3 Tests
- `test_export_image_cloud_altaz_grid` — run export with the flag and verify output image has expected cloud content.
- `test_export_image_cloud_altaz_grid_cache` — second export with same args reuses grid cache.

### 6.4 Validation
- `python -m zstarview.zstarview_export_image <location> --cloud-altaz-grid true` produces a PNG with clouds.
- Timing: second export is faster due to grid cache.

---

## 7. Phase 5 — Polish and Legacy Removal (3–4 sessions)

### 7.1 Stripe rendering from grid

- Port existing `render_variable_width_cloud_stripes_rgba` / `_render_alpha_scaled_cloud_stripes_rgba` in `gui/composite.py` to accept `CloudAltAzGrid` instead of `CloudAmountField`.
- The screen-space stripe phase and sampling logic remain similar; only the input field changes from normalized `(u, v)` to `(alt, az)`.

### 7.2 Performance tuning
- Profile grid build time vs. old `render_from_source_with_coverage`.
- If ingestion is too slow, downsample the geographic sampling step or parallelize shell loops.
- If circle render is too slow, switch from gaussian splat to stencil circles or vectorize with NumPy.

### 7.3 Remove legacy path

Once the alt/az grid path is validated for both GUI and export:
- Remove `CloudDisc.render_from_source_with_coverage` (or move to `legacy.py`).
- Remove `CloudAmountField` based stripe path if fully replaced.
- Update `CloudImageState` to always carry `altaz_grid` and drop the legacy `image` field.
- Remove the `--cloud-altaz-grid` feature flag and make the grid path the default.

### 7.4 Documentation updates
- Update `docs/specification.md` user-facing cloud description.
- Update `docs/cli-overlays.md` if new CLI options are added or defaults change.
- Update `docs/design/rendering-pipeline.md` if details changed during implementation.
- Write a `release-notes.md` entry.

---

## 8. Testing Strategy

### Unit tests (focus)
- Coordinate transforms: lat/lon/height → alt/az, alt/az → screen.
- `CloudAltAzGrid` build from synthetic BT fields.
- Circle renderer output shape, alpha scaling, overlap.
- Cache key stability and roundtrip.

### Integration tests
- Mock `CloudSourceData` through the new pipeline end-to-end.
- GUI signal flow with a fake `CloudController`.
- Export-image CLI with mocked satellite data.

### Visual regression
- Capture reference renders from the old path and the new path for a few locations/times.
- Compare perceptually (SSIM or manual inspection) until they are acceptably close.

### Manual checks
- Rotate view continuously; verify no stutter from cloud re-projection.
- Zoom in/out; verify cloud density/scale adapts.
- Offline mode with cached grid; verify clouds still render.

---

## 9. Backward Compatibility and Rollback

- Keep `CloudDisc.render_from_source_with_coverage` untouched until Phase 5.
- Add `CloudAltAzGrid` side-by-side; do not modify existing `CloudSourceData` required fields.
- Use a feature flag so users can immediately disable the new path if issues arise.
- Do not delete `cloud_amount_field` / stripe rendering until the grid-based stripe renderer is fully equivalent.
- If the new path fails, `CloudController` should fall back to the legacy path and emit a warning banner.

---

## 10. Risks and Mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| Alt/az grid ingestion slower than old render | High | Build in worker process; cache aggressively; tune `geo_sample_step_deg`. |
| Circle renderer looks too sparse or blocky | Medium | Gaussian splat; tune radius/opacity; add smoothing. |
| Azimuth wrap bugs near 0°/360° | Medium | Write dedicated tests; use smallest-angular-distance logic. |
| Missing-data tint misaligned | Low | Use same alt/az → screen transform; test against old path. |
| Geo-satellite path divergence | Low | Leave it on existing path initially. |
| Feature flag complexity | Low | Make flag default false, then flip to true after validation. |

---

## 11. Definition of Done

- [ ] `CloudAltAzGrid` builds correctly from GOES and Himawari source data.
- [ ] Grid is cached on disk and reused across camera changes and repeated exports.
- [ ] Circle renderer produces a visually acceptable cloud overlay.
- [ ] GUI view rotation/zoom is smooth with no re-download/re-ingestion.
- [ ] `--cloud-altaz-grid true` works in both GUI and export-image.
- [ ] Legacy path still works with the flag disabled.
- [ ] Stripe rendering from the grid is implemented or a follow-up issue is filed.
- [ ] All new code has tests; CI passes.
- [ ] Design/spec docs are updated to match final behavior.
- [ ] Release notes entry is added.

---

## 12. Suggested Next Action

Start **Phase 0** by reading `clouddisc/projectors/az.py`, `clouddisc/core.py`, and `gui/cloud_controller.py`, then commit the module skeletons and constants.
