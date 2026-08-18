# Session 2026-08-18

Scope: Investigate thin edge noise reported after star interpolation mesh rendering.

- Topic: Thin noise at the four viewport edges after star mesh interpolation
  - Finding: `shot1.png` shows narrow irregular bright content along the outer
    rectangular viewport boundary, distinct from the circular sky-disc edge.
  - Finding: The current path introduced by `69f3b6ba` renders each mesh cell
    into a transparent intermediate `QImage`, then scales that image to the
    viewport. The source image is sampled right up to its boundary while the
    mesh target boundary is displaced by sidereal motion.
  - Validation: `uv run -p .venv/bin/python pytest -q tests/test_star_interpolation.py`
    passes (10 tests). Existing tests cover affine motion and composition, but
    do not assert clean outer-edge pixels for a displaced mesh.

- Topic: Add a quick guard-band mitigation
  - Decision: Add a transparent guard band equivalent to approximately 16
    final viewport pixels around the low-resolution star surface during mesh
    warping. Offset both source and target mesh coordinates into that expanded
    surface, then crop the original image rectangle when compositing.
  - Rationale: This prevents displaced mesh boundaries from sampling directly
    at the cached QImage edge while preserving the original viewport size.
  - Validation: `uv run -p .venv/bin/python pytest -q tests/test_star_interpolation.py`
    passes (10 tests).

- Topic: Document the expanded mesh design
  - Decision: Specify one shared expanded coordinate system for the star
    surface and all star-following mesh users, with a final viewport clip that
    hides the guard region.
  - Rationale: A surface-only overscan mesh would require separate coordinate
    handling for bright stars, asterisms, twinkle, and labels. A shared
    expanded mesh keeps their motion aligned and prevents source-edge artifacts.
  - Documentation: Updated `docs/design/rendering-pipeline.md`,
    `docs/design/gui-screen-update-and-cache.md`, and the user-visible timing
    section of `docs/specification.md`.

- Topic: Correct final composition of the guarded transformed surface
  - Finding: Cropping the transformed image at the original source rectangle
    exposed transparent side bands after the mesh displaced the content.
  - Decision: Clip the painter to the real viewport, then draw the complete
    guarded transformed image with its guard offset and low-resolution-to-
    viewport scale applied.
  - Validation: `uv run -p .venv/bin/python pytest -q tests/test_star_interpolation.py`
    passes (10 tests).
