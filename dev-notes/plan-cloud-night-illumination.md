# Voxel cloud nighttime illumination implementation plan

Status: Implementation and validation completed (2026-10-10).

## Scope and behavior

Replace the cloud-amount-scaled nighttime tint with the sum of weak ambient,
directional lunar, and vertical ground illumination. Preserve daytime shading,
view-ray extinction, cutout styling, shell rendering, and star ordering.
No dependencies, CLI options, package version changes, or publication steps.

## Implementation steps

1. Add shared appearance constants and a pure lunar-strength helper. Interpret
   the existing phase_angle as solar elongation: 0 degrees is new moon and
   180 degrees is full moon. Use a display approximation with a smooth horizon
   fade and zero intensity below the horizon. Include all celestial bodies,
   regardless of viewport visibility, when assembling inputs.
2. Extend native and Geo-satellite shading with optional lunar and ground inputs.
   Reuse actual-distance ray traversal for moon rays and vertical downward rays
   from cell centers. Include the target cell's half-path extinction. Cache
   illumination for visited cells. Sum nighttime RGB contributions, apply a
   shared exponential display curve, then blend with daytime color using the
   existing solar-altitude transition. Do not scale these colors by cloud amount
   again or change view-ray opacity.
3. Add a worker-only ground-light sampler using existing night-light tile
   acquisition and validation. Map each layer's physical cell centers into
   geographic coordinates with the same local ENU approximation. Sample the
   underlying raster pixel directly, without a horizontal blur or emitter
   integration. Read raster windows in batches rather than one I/O operation per
   cell. Use an absolute display normalization, independent of scene maximum.
   Missing tiles contribute zero, while other tiles and moonlight remain usable.
   Fetch each tile independently with bounded network timeouts and keep the
   full operation inside the existing subprocess deadline. Skip night-light
   acquisition entirely when the nighttime blend is zero.
4. Bump the private worker protocol and pass lunar position/elongation through
   shared styling, both GUI controllers, and image export. Add lunar state to
   native render request identity so latest-request handling stays correct.
   Ground raster data is read afresh in each one-shot worker; no persistent
   illumination cache is introduced. Use existing background controller paths
   and process isolation, avoiding extra native work on GUI threads.
5. Validate pure behavior with small synthetic volumes: lunar phase and horizon,
   lunar shadows without changed view transmission, vertical ground shadows
   across shifted grids, ground spatial differences, ambient behavior, twilight
   continuity, and native/Geo input propagation. Mock network and tile failures;
   verify raster sampling with small local GeoTIFFs. Generate a synthetic visual
   comparison artifact for dark, moonlit, and ground-lit clouds.
6. Run focused tests, Ruff and type checking, then the full test suite as
   appropriate. Update specification/design to describe implemented behavior
   and constants, and append decisions and validation results to the session log.

## Limits

This remains a display model, not calibrated radiative transfer. Lunar eclipse,
terrain shadowing, ground color inferred from lamps, Earth curvature beyond the
existing local geometry, and multiply scattered light are outside this change.
Appearance constants will be checked against synthetic comparisons; real
satellite scenes may require later visual tuning.

## Completion

- Added the shared lunar/display helpers, direct raster sampler, nighttime ray
  extinction and tone mapping, protocol version 5, and GUI/export inputs.
- Added illumination refresh scheduling and completion-key validation, including
  the first celestial state arriving after cloud initialization.
- Updated specification and design from unimplemented proposals to implemented
  behavior, with numeric constants and acquisition/deadline limits documented.
- Final headless regression run: 1713 passed, 5 skipped. All changed files pass
  Ruff; the new lighting/sampling modules and tests pass focused mypy and
  compile checks. Repository-wide mypy still has baseline errors; no new
  diagnostics are introduced relative to HEAD.
- Synthetic comparison generated at /tmp/zstarview-cloud-night-comparison.png.
  Regenerate using the comparison script; no live satellite visual validation
  was performed. No commits, tags, pushes, or version bump.
