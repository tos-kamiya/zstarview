# Session 2026-07-11

Scope: Evaluate data classes for reducing the argument count of the night-light glow calculation.

- Topic: Night-light calculation arguments
  - Decision: Introduce separate immutable configuration and request/input data classes rather than extending `NightLightTerrainContext` with settings.
  - Rationale: `compute_night_light_glow_profile()` currently mixes observer/time inputs, terrain observations, download/cache policy, and numerical tuning parameters. Terrain context is observation data and should remain separate from reusable settings.
  - Proposed shape: `NightLightSettings` for refraction, tile inclusion, cache/download timeouts, maximum distance, and distance step; `NightLightRequest` for observer coordinates, sun altitude, and `NightLightTerrainContext`.
  - Compatibility: Migrate the two production callers (`gui/sky_worker.py` and `cli/export_image.py`) and keep the public function as a small wrapper if compatibility is desired.

- Topic: Internal grouping implementation
  - Decision: Keep the existing public keyword-only function unchanged, construct `_NightLightRequest` and `_NightLightSettings` inside it, and dispatch through `_compute_night_light_glow_profile()` using a typed common-kwargs boundary.
  - Rationale: This reduces the main calculation's conceptual argument surface without changing either production caller or the existing cache-keyed calculation functions.
  - Verification: mypy passed for `night_lights.py`, Ruff passed for the changed module and its tests, and `tests/test_night_lights.py` passed (`18 passed`).

- Topic: Remove unused night-light helpers
  - Decision: Remove `_terrain_sample_distances_key`, `_sample_ray_brightness`, `_sample_ray_brightness_curve`, and the test that exercised only the unused brightness-curve path.
  - Rationale: None had production callers; terrain sample keys already use `_float_sequence_key`, while active sampling uses `_sample_ray_night_light_samples` directly.
