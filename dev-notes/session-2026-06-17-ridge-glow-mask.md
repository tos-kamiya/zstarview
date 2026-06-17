# Session 2026-06-17

Scope: ridge glow masking update.

- Topic: Ridge glow should be clipped to the main ridge silhouette.
  - Decision: Change the glow compositor so the ridge-glow contribution is masked against the main terrain profile, not the cumulative secondary-ridge envelope.
  - Rationale: Ridge glow is meant to emphasize the primary horizon line. Allowing the secondary-ridge envelope to drive the cutoff lets the glow disappear in places where the main ridge is still present, which weakens the intended silhouette emphasis.
  - Notes: Keep the night-light contribution separate so it can still use the existing terrain-aware horizon handling, while the ridge-glow contribution gets its own main-ridge mask and opacity scaling.

## Summary

- Updated the glow compositor so ridge glow is clipped against the main terrain profile.
- Split ridge glow opacity from the night-light contribution in the shared alpha field.
- Validation: `uv run -p .venv/bin/python pytest tests/test_glow_mask.py -q` and `uv run -p .venv/bin/python ruff check src/zstarview/gui/composite.py tests/test_glow_mask.py`

- Topic: Link ridge glow opacity to night-light opacity
  - Decision: Make the effective ridge glow opacity follow the night-light opacity throughout the render pipeline and the window/CLI style builders.
  - Rationale: Ridge glow is a helper for the night-light silhouette, so it should not drift independently when the night-light control changes. Linking the values keeps the visual balance predictable and removes a separate tuning axis that is no longer useful.
  - Validation: `uv run -p .venv/bin/python pytest tests/test_window_render_sync.py -k "links_ridge_glow_to_night_light or skips_night_lights_while_press_pending" -q`, `uv run -p .venv/bin/python pytest tests/test_export_image_layer_gating.py -q`, and `uv run -p .venv/bin/python ruff check src/zstarview/render/pipeline.py src/zstarview/cli/export_image.py src/zstarview/gui/window_render.py src/zstarview/gui/window_inputs.py tests/test_window_render_sync.py tests/test_glow_mask.py`

- Topic: Reduce glow height scales
  - Decision: Lower the night-light height scale from `36.0` to `20.0` and the ridge-glow height scale from `7.5` to `4.0`.
  - Rationale: Both effects were reading too tall relative to the intended horizon-adjacent look. Tightening the height scales makes the glow sit closer to the ridge silhouette and reduces the sense of a vertical column.

- Topic: Restore glow height scales
  - Decision: Return the night-light height scale to `30.0` and the ridge-glow height scale to `6.0`.
  - Rationale: The lower values made the glow too hard to distinguish from low clouds. A slightly taller profile keeps the ridge emphasis while avoiding overlap with cloud structure.
