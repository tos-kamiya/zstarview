# Session 2026-07-11

Scope: Audit whether currently single-file modules should be split after the Atlas separation work.

- Topic: File-boundary audit
  - Findings: The strongest split candidate is `src/zstarview/night_lights.py` (about 69 KB): dataset manifest/download/cache handling, terrain sampling/key construction, and glow-field/profile computation are distinct subsystems spanning roughly lines 151-294, 297-1116, and 1117-1630.
  - Findings: `src/zstarview/cli/export_image.py` (about 81 KB) combines export orchestration, asynchronous layer acquisition, render-input construction, image encoding, terminal/sixel output, and search reporting. Its fetch/render/output boundaries are clear.
  - Findings: `src/zstarview/render/pipeline.py` (about 52 KB) contains the scenic pipeline and Atlas instrument presentation in one module. The Atlas-specific presentation class and helpers are already contiguous enough to extract, but shared private drawing helpers make this a medium-risk refactor.
  - Findings: `src/zstarview/gui/window.py` (about 128 KB) is a large stateful coordinator, but it already delegates rendering, updates, controllers, widgets, and state to separate modules. More mixins would likely obscure shared state; no immediate split is recommended without a behavior-focused boundary.
  - Findings: `src/zstarview/gui/composite.py` (about 100 KB) contains tightly coupled cloud, glow, and cache logic. It is a future candidate, but the central `SkyCompositorCache` makes an immediate split riskier than the three candidates above.
  - Test maintenance: `tests/test_window_render_sync.py` is about 212 KB and mixes base-scene, overlay, HUD, and interaction tests; splitting it by behavior would improve navigation without changing production architecture.

- Topic: Split `night_lights.py` by responsibility
  - Decision: Extract NASA tile/cache/GeoTIFF sampling into `night_light_source.py`, terrain-derived edge input into `edge_glow.py`, and shared configuration constants into `night_lights_constants.py`. Keep `night_lights.py` as the profile-calculation facade so existing imports and the combined profile contract remain stable.
  - Rationale: This reduces the main module from 1,609 to 1,319 lines without changing the renderer-facing `NightLightGlowProfile` or the existing calculation API. A complete model split would affect compositor and state plumbing, so it remains outside this refactor.
  - Validation: `ruff check` passed; targeted glow/night-light tests passed (30 tests); full suite passed (1,234 tests); direct `.venv/bin/mypy` passed for all four changed source modules.

- Topic: Move ridge glow distance correction
  - Decision: Move `_night_light_distance_boost()` and `_ridge_glow_distance_gain()` into `edge_glow.py` alongside `_terrain_sample_edge_strength_rows()`. Re-export them from `night_lights.py` for compatibility with existing private-function tests.
  - Rationale: These helpers only shape terrain-derived ridge glow and do not belong to NASA night-light sampling. Shared alpha-field construction remains in `night_lights.py` because both layers use it.
  - Validation: Ruff, mypy, and `tests/test_night_lights.py` (17 tests) passed.

- Topic: Remove unused terrain band mask helpers
  - Decision: Delete `_terrain_band_target_mask()` and `_terrain_band_target_altaz_mask()` together with their test-only coverage.
  - Rationale: Neither helper is referenced by production code; the active field builder performs the terrain fade directly from the threshold grid.
  - Validation: No remaining references; Ruff and mypy passed; `tests/test_night_lights.py` passed (15 tests).
