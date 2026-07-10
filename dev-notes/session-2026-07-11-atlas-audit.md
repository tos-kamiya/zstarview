# Session 2026-07-11

Scope: Audit for unused data structures, constants, and functions left behind while separating Atlas from zstarview.

- Topic: Audit result
  - Decision: No Atlas-specific orphan module, old object-viewer symbol, or unused Atlas style/data structure was confirmed. Keep the shared instrument presentation and theme structures.
  - Findings: `apply_atlas_profile()` repeats the `light_background_star_outline=True` default already present in `_ATLAS_DEFAULTS` (`src/zstarview/gui/atlas.py:29,44-45`). `docs/design/atlas.md` still says the Atlas cloud color connection and instrument cloud connection are unfinished (`:299,302`), although the implementation is connected and targeted Atlas tests pass.
  - Separate cleanup candidate: `_resolve_frame_context()` in `src/zstarview/render/pipeline.py:224` has no caller. It is not Atlas-specific and should be removed or deliberately retained in a separate cleanup change.
  - Verification: Ruff passed; targeted Atlas/theme tests passed (`17 passed`). Vulture only reported the unrelated `_resolve_frame_context()` plus public/module attributes and intentionally unused compatibility parameters. The broader Qt render test collection aborted in the headless environment before tests ran.

- Topic: Apply audit cleanup
  - Decision: Remove the duplicated Atlas star-outline assignment, remove the unused `_resolve_frame_context()` helper, and mark the implemented Atlas cloud steps as complete in the design document.
  - Verification: Ruff passed; targeted Atlas/theme tests passed (`17 passed`); `git diff --check` passed.
