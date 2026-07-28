# Session 2026-07-28

Scope: Review optimization opportunities in AKARI IR bands image generation; no runtime implementation.

- Topic: AKARI image-generation optimization review
  - Decision: Assess the current projection, sampling, and compositing path before proposing changes. Keep this session limited to design-level recommendations.
  - Rationale: The user requested consideration only, and the largest risks are likely to be visual correctness and cache invalidation rather than raw file loading alone.
  - Finding: The display asset is loaded once per process, but the runtime overlay is regenerated at the painter viewport size and performs coordinate conversion, nearest-neighbor texture lookup, gamma, value-knee, and color composition on each base render.
  - Finding: The current implementation has only a focused value-knee test; the projection, seam, cache-reuse, and image-quality behavior are not directly covered by dedicated tests.
  - Recommendation: Prioritize measurement and a projected-overlay cache keyed by observer/time/view/geometry/render-size parameters. Consider precomputed packed-color or mip-level assets only after measuring whether the projection or post-processing dominates.

- Topic: Opt-in AKARI runtime profiler
  - Decision: Add `ZSTARVIEW_PROFILE_AKARI=1` instrumentation to the runtime overlay generator. Log one profile record every ten successful overlay generations, split into asset load, screen projection, coordinate conversion, texture lookup, color processing, and total time.
  - Rationale: A 2000x2000 real GUI render includes the actual viewport and compositor behavior, which is more representative than a synthetic microbenchmark. Logging every tenth call limits noise while preserving repeated-render observations.
  - Validation: Targeted molecular-cloud test, Ruff, and `git diff --check` passed.

- Topic: AKARI runtime measurement at approximately 2000x2000
  - Finding: First measured generation took 1067.87 ms, including 140.41 ms for the asset load. A later same-size generation took 850.95 ms with a 0 ms asset load. A resized 1885x1803 generation took 798.09 ms.
  - Finding: In the cached 2000x2000 run, color processing was 315.63 ms, coordinate conversion 233.56 ms, texture lookup 160.96 ms, and screen projection 139.20 ms. Color processing is the largest individual stage; coordinate conversion is the next largest.
  - Decision: Do not prioritize NPZ format or asset loading. Evaluate packed/precomputed display RGB or lookup-table processing first, then reduced-resolution rendering or projection reuse. Preserve full-resolution visual comparison before adopting either change.

- Topic: Candidate AKARI render scale
  - Proposal: Render the AKARI overlay at 1/4 of the viewport width and height, matching the existing sky-disc internal render scale, then upscale it with smooth filtering for final compositing.
  - Rationale: This reduces the generated pixel count to 1/16. At a roughly 90-100 degree FOV, a 2000-pixel viewport corresponds to about 0.045-0.05 degrees per pixel, while the 1/4 scale is about 0.18-0.20 degrees per pixel, close to the 2048x1024 AKARI display asset's nominal angular sampling. The user does not require fine filament detail.
  - Risk: Nearest-neighbor enlargement would make the overlay blocky and may expose seams. Smooth enlargement and a visual comparison at the Galactic-plane boundary are required.

- Topic: Implement quarter-resolution AKARI rendering
  - Decision: Generate the AKARI overlay at 25% of the final width and height, scale its geometry accordingly, then restore the final RGBA dimensions with Qt smooth transformation before compositing.
  - Rationale: Keep the compositor's existing full-size array contract while reducing projection, coordinate conversion, lookup, and color work to approximately one sixteenth of the previous pixel count.
  - Validation: AKARI overlay tests, Ruff, and `git diff --check` passed. A broader export-layer test run has four pre-existing fixture failures because its `_Args` test double lacks the already-required `akari_ir_bands_opacity` attribute; the twelve remaining tests passed.

- Topic: Post-implementation AKARI measurement
  - Finding: At 2000x2000, quarter-resolution generation measured 243.60 ms including 129.12 ms initial asset loading; the non-loading stages totaled about 114 ms. This is substantially below the previous 850.95 ms cached full-resolution result.
  - Finding: An interactive 600x600 fast-render sample measured 12.61 ms, with only 150x150 generated pixels.
  - Finding: One 600x579 sample measured 550.39 ms despite only 150x145 generated pixels; its lookup and color stages were disproportionately high. Treat this as a transient system/render-thread contention outlier until repeated under idle conditions.
  - Decision: Keep the quarter-resolution implementation. No further optimization is justified from the current data; if needed, repeat the resize measurement several times after the window settles.

- Session summary
  - Decision: Fix AKARI rendering at quarter width and quarter height with smooth final upscaling. Remove the temporary runtime profiler after confirming the optimization in the GUI.
  - Validation: AKARI overlay tests and Ruff passed; `git diff --check` passed. The broader export-layer run retained four fixture failures unrelated to this change because the test double lacks `akari_ir_bands_opacity`.
  - Follow-up: None for the AKARI render scale unless future visual review finds unacceptable blur or boundary artifacts.
