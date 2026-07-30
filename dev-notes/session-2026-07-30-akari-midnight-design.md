# Session 2026-07-30 — AKARI midnight illumination design

Scope: Design-only review of raising AKARI IR bands visibility during the
deep-night period. No runtime implementation was made.

- Topic: Time-dependent balance between artificial light and AKARI IR bands
  - Decision: Treat the clock-time effect as a shared, smooth deep-night
    profile, but apply separate coefficients to night light, building/urban
    outlines, and AKARI. Keep ridge glow/edge glow outside this artificial-
    light profile because it represents natural atmospheric light.
  - Rationale: The source discussion distinguishes artificial illumination,
    natural atmospheric glow, and astronomical color. The current renderer
    already has these as separate paths: night light/ridge glow in the sky
    compositor, AKARI as an additive overlay, and urban outlines as a later
    terrain layer.
  - Proposed initial profile: night light relative strength 1.00 at 19–21,
    0.85 at 22:00, 0.65 at 00:00, 0.50 at 02:00, 0.45 at 04:00, and
    0.50–0.60 near 05:00, interpolated continuously. Apply a modest AKARI
    lift, approximately 1.00 outside deep night and no more than about 1.20–
    1.25 near the minimum. Apply a weaker 0.70–0.85 minimum range to urban
    outline opacity so geographic legibility is retained.
  - Constraint: The current VIIRS/night-light raster is a static spatial
    source and may already resemble a late-night observation. The proposed
    clock curve should therefore be treated as a visual activity correction,
    not as a claim that the raster contains hourly radiometry. Calibration
    against the current default appearance is required before selecting the
    final minimum.

- Topic: Time basis and composition boundaries
  - Decision: Use the observer's local civil time for the first design,
    with the existing timezone name, and keep solar-altitude gating as a
    separate boundary condition. A future solar-time option can be considered
    if locations far from their time-zone meridian show a noticeable issue.
  - Rationale: Artificial-light schedules follow local civil time, while
    sun altitude controls whether a night layer should exist at all. Combining
    these into one factor would make the behavior difficult to reason about,
    especially during polar day/night.
  - Boundary: The clock modulation is relevant only while the layer is in
    night mode; the existing sun-altitude fade toward civil twilight remains
    authoritative.

- Topic: Validation before implementation
  - Decision: Compare at least 19:00, 22:00, 00:00, 02:00, 04:00, and 05:00
    for a city, a dark rural location, and a high-latitude case. Inspect the
    Milky Way/AKARI boundary, star visibility, urban outline readability,
    and the absence of a visible step at midnight or the profile minimum.
  - Rationale: The intended effect is a change in layer balance rather than
    an absolute physical brightness reconstruction. Visual regression at
    these representative times is more informative than RGB equality alone.

- Topic: Separate building fill from building outlines
  - Decision: Apply the same deep-night reduction used for night light to the
    filled roof polygons only. Keep the roof/building outline strokes on the
    weaker urban-legibility reduction (or otherwise substantially brighter).
  - Rationale: The filled polygon reads as an illuminated surface in the
    current visual language, while the outline is a geographic reference.
    This preserves the user's distinction without treating the entire
    building layer as artificial illumination.
  - Decision detail: Multiply the final roof-fill contribution directly by
    the night-light time coefficient. The existing nonzero fill-alpha floor
    is not exempt; the fill may disappear at the deep-night minimum. The
    outline remains as the readability-preserving reference.

- Implementation summary
  - Added the local-clock `night_activity_factor` with the documented smooth
    profile and kept the existing solar-altitude twilight gating separate.
  - Applied the factor to night-light opacity and roof-polygon fill only;
    mapped the default AKARI final opacity continuously from `0.10` to
    `0.15`. Edge/ridge glow and building strokes remain independent.
  - Added tests for local-time factors, AKARI opacity lift, and disappearing
    roof fill with surviving outlines.
  - Validation: 74 targeted tests passed; F401/F821 Ruff checks, compileall,
    and `git diff --check` passed.

- Topic: Time-dependent sky ambient floor
  - Decision: Use the same local-clock deep-night progress as the artificial
    light and AKARI layers to interpolate the sky-disc ambient RGB from
    `[1, 2, 5]` in ordinary night to `[2, 4, 10]` around 02:00–04:00.
  - Rationale: The ambient floor is visually separate from AKARI and night
    light, but a small deep-night lift helps reveal low-contrast IR structure
    without increasing the AKARI overlay opacity.
  - Validation: 71 sky-disc, night-light, and startup tests passed after
    adding the time-dependent ambient coverage.

- Topic: AKARI display gamma tuning
  - Decision: Set the production AKARI display gamma to `0.35`; defer test
    expectation changes while the visual tuning remains provisional.
