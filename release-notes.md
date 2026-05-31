# zstarview Release Notes

## 1.29.1 - 2026-05-31
- Refined the tropical cyclone overlay into a cone-shaped outer contour, with a `maxwind_kt`-driven base radius, a `400km` visibility cutoff, and the label anchored below the cone tip.
- Lowered the default tropical cyclone opacity to `0.25` while keeping `0.0` as the explicit disable value for cyclone fetch and drawing.
- Updated the cyclone specification, design notes, and overlay CLI docs to match the current render behavior.

## 1.29.0 - 2026-05-30
- Added a tropical cyclone opacity option that also disables cyclone API fetches when set to `0.0`, and wired the value through both GUI and export-image flows.
- Reordered the bottom status line so cyclone status now appears right after clouds, and updated the README overlay descriptions to follow the Earth guide, terrain horizon, water surface, urban outline, and night lights order.

## 1.28.12 - 2026-05-30
- Removed the temporary tropical cyclone scaffolding in the window and updates paths, switching the overlay refresh flow to direct method calls instead of defensive `getattr` lookups.
- Kept the cyclone render projection on the draw path and cleaned up the release/update wiring so the overlay refreshes without the earlier temporary hooks.

## 1.28.11 - 2026-05-30
- Fixed the tropical cyclone overlay so the projected storm position advances over time instead of staying pinned to the cached advisory point.
- Removed the temporary cyclone debug prints and kept the render cache keyed on the projected storm state so the display updates when the projection changes.

## 1.28.10 - 2026-05-28
- Kept `zstarview-gui` aligned with the CLI-only Geo-satellite behavior, so the startup dialog no longer exposes that toggle while the internal GUI path still honors the feature when enabled from the command line.
- Refreshed the bundled Geo-satellite gray-common mask assets and documented the regenerated six-frame workflow used to rebuild them.

## 1.28.9 - 2026-05-28
- Added an experimental, opt-in Geo-satellite cloud path for the Europe workflow band, with GUI wiring, status handling, and `--geo-satellite` CLI docs.
- Improved urban-outline pruning and status reporting, including the new `--urban-outline-max-candidates` knob, split base/skyscraper counts, and updated overlay docs.
- Updated the README credits/license section for MET Norway Geo-satellite attribution and clarified the opt-in behavior in the main README.

## 1.28.8 - 2026-05-28
- Fixed startup overlay handling so explicit `-a 0` and `-c 0` disable aircraft and cloud overlays for the session, even when Geo-satellite is enabled.
- Updated the CLI and design/specification docs to reflect the explicit-disable behavior.

## 1.28.7 - 2026-05-25
- Restored fixed-interval cloud refresh scheduling so sky updates no longer keep postponing the next cloud download.

The entries below were generated later from `releases.txt`, git history, and `pyproject.toml` diffs.
It treats the first commit that changed `src/zstarview/__about__.py` to a released version as that release commit.

## 1.28.6 - 2026-05-24
- Focused on GUI startup flow. Added a **city auto button**, **launcher startup profile**, and **transparent theme opacity presets** to the startup experience.
- Allowed both time-mode checkboxes to be turned off, making startup configuration more flexible.

## 1.27.13 - 2026-05-22
- Improved interactive updates with **finer view direction control**, **idle-gated scheduler refreshes**, and **draw-time rendering for satellites and aircraft**.
- Split cloud handling into separate fetch and projection stages to simplify the update flow.
- Packaging: Removed the `zstarview-debug` entry point.

## 1.26.1 - 2026-05-17
- Cleaned up observer-height handling with a **`height-add-m` alias** and a more compact location height summary.
- Fixed water tile open counting and prevented water-surface mode from being used together with terrain horizon rendering.

## 1.25.11 - 2026-05-16
- Expanded water overlay processing with **observer-height-aware scans**, **overlay cancellation**, **runtime thinning**, and **raw SVG preview output**.
- Packaging: Added the `zstarview-gui` GUI entry point.

## 1.23.8 - 2026-05-11
- Stability-focused release for cloud rendering and startup flow. Fixed cloud opacity compositing, stale sky snapshots, startup input gating, and source-first cloud rendering.

## 1.23.5 - 2026-05-10
- Added a **night light overlay**, cached the base profile, and included **cloud coverage in export summaries**.
- Also adjusted DSO hover highlight colors.

## 1.22.3 - 2026-05-08
- Improved usability with a **CLI Alt/Az keep toggle**, **fit-to-screen action**, **map-like sky disc grid**, and **never-rises opacity control**.

## 1.21.9 - 2026-05-01
- Introduced the **transparent theme** and added **terrain ridge preview**, **secondary ridge overlay**, and **tall-building viewpoints**.

## 1.21.2 - 2026-04-28
- Switched the earth guide to a hatch-line style and tuned celestial reference line sampling, label placement, and the startup black-frame issue.

## 1.20.20 - 2026-04-25
- Added an **outline mode for bright bodies**, **direction-grid hover guides**, **direction-grid export-image support**, and a **shutdown status banner**.

## 1.20.15 - 2026-04-22
- Expanded display controls with **edge FOV projection control**, a **credits link in the help menu**, **square client-area resize**, and **visibility boost controls**.

## 1.20.7 - 2026-04-19
- Expanded JPL and Horizons integration with **spacecraft overlays**, **body-search fallback**, **small-body samples**, and **persistent refresh**.

## 1.9.18 - 2026-04-17
- Fixed GUI-script ephemeris download on Windows and updated overlay label color tests.

## 1.9.17 - 2026-04-14
- Added **overture release cache checks** and adjusted the white theme background, satellite overlay opacity, and earth-guide rendering during interaction.

## 1.9.10 - 2026-04-12
- Added **`--observation-info`**, an **observation-info overlay mode**, and **auto location resolution through `ip-api.com`**.
- Added a 3-second rate limit for IP-API access.

## 1.9.5 - 2026-04-12
- Added **decorated window mode** and **octilinear earth-guide preview**, expanding the earth-guide presentation.
- Lowered the minimum jump/search target view center to `-5` degrees.

## 1.8.11 - 2026-04-05
- Added a **building-top observer option** and tuned theme text, hover text, compass label placement, and the night-style status line.

## 1.8.8 - 2026-04-02
- Added **tower review candidates** and improved cloud stripe rendering and two-shell-height blending.

## 1.8.2 - 2026-03-30
- Added **guideline opacity control** and refined the observation overlay and window-edge chrome.

## 1.7.4 - 2026-03-29
- Fixed light-theme text, stale clouds during view changes, and splash behavior on location errors.
- Also started caching cloud samplers per source.

## 1.7.0 - 2026-03-27
- Added **export CLI cache-dir support**, **CLI version and overlay color constants**, and **famous mountain and tower viewpoints**.

## 1.6.5 - 2026-03-25
- Added **export metadata summaries**, a **satellite-element cache layer**, **OpenSky aircraft snapshot caching**, and expanded **station satellite overlays**.
- Packaging: Simplified GUI launch entry points and moved `zstarview` to `zstarview.gui.viewer:main`.

## 1.4.1 - 2026-03-22
- Added the **transparent theme** and tuned transparent-sky and faint-star rendering.
- Also skipped incomplete Himawari slots and removed the `pyresample` runtime dependency.

## 1.3.3 - 2026-03-21
- Added the **export-image CLI**, **sixel export**, and **PNG-to-stdout export**.
- Switched cloud ingestion to direct Himawari ISatSS tile input.
- Packaging: Added `zstarview-export-image` and cleaned up related entry points.

## 1.2.2 - 2026-03-20
- Added a **shared content-FOV option** and fixed stale disc renders, ocean DEM horizon fallback, and cloud hatch compatibility.

## 1.1.1 - 2026-03-19
- Added **OpenSky aircraft overlays** and an **aircraft-layer opacity toggle**.
- Also tuned cloud horizon cutoff and aircraft overlay motion.

## 1.0.6 - 2026-03-17
- Added **skyscraper urban outlines** and **tile seed data**, and consolidated Overture building layers.
- Also corrected urban-outline minimum-height handling at draw time.
- Packaging: Split `astropy` dependencies by Python version.

## 1.0.4 - 2026-03-16
- Restricted urban-outline cache loading and refreshed related documentation images.

## 1.0.3 - 2026-03-15
- Fixed off-screen planet labels and hamburger-menu visibility.
- Added an urban-outline screenshot gallery and cleaned up unused data.
- Packaging: Updated the development-status classifier from Beta to Production/Stable.

## 1.0.1 - 2026-03-15
- Added a **PLATEAU CityGML zip importer** and **import utility**, and switched urban outlines to an Overture-based pipeline.
- Also introduced height-based urban-outline opacity scaling.
- Packaging: Added the `zstarview-import-overture-buildings` entry point.

## 0.29.1 - 2026-03-14
- Added **Nominatim-based startup place search** and reordered overlay location info.
- Also continued cleaning up the legacy observer-height fallback.

## 0.28.4 - 2026-03-14
- Updated terrain and asterism stroke scaling, viewpoint-height overlay display, and cross-version data-tool compatibility.

## 0.28.0 - 2026-03-12
- Added **urban-outline layer controls**, an **urban debug layer**, and **derived tile indexes**.

## 0.27.1 - 2026-03-10
- Added **bundled mountain viewpoints**, **explicit viewpoint prefixes**, and **mountain dataset query/tooling support**.

## 0.24.4 - 2026-03-08
- Added **keyboard shortcuts for menu actions** and refined the view menu, sky-disc fallback, and terrain-opacity behavior.

## 0.24.1 - 2026-03-08
- Added a **terrain-horizon overlay**.
- Packaging: Tightened and reorganized dependency version constraints.

## 0.23.3 - 2026-03-07
- Added **boolean startup visibility options**, **hover-rotating asterism overlays**, and **great-circle asterism drawing**.

## 0.22.1 - 2026-03-05
- Added **DSO catalog generation**, **DSO overlays**, **dual hover labels**, and a **DSO visibility toggle**.
- Also added a CLI option for cloud missing-tint opacity.
- Packaging: Added `pyongc` and `mypy` as development dependencies.

## 0.20.2 - 2026-03-04
- Added the **`-w` alias for render width**, **expected-width star scaling**, and **window-geometry startup options**.
- Reduced the default sky update interval from 180 seconds to 60 seconds.

## 0.18.2 - 2026-03-02
- Introduced a 2nd-magnitude star baseline and diamond overlay, and smoothed the appearance of single-pixel stars.
- Also continued cleanup and documentation work around sky layers and helpers.

## 0.18.1 - 2026-03-02
- Precomputed star-render cache size and color factors, and tuned tiny-star culling outside the background FOV and hover-label conditions.

## 0.17.3 - 2026-03-01
- Added **extra10 split catalogs**, a **merged star-catalog generation pipeline**, a **named-star search menu**, and **famous-star jump support**.

## 0.15.0 - 2026-03-01
- Improved planet visibility and bloom rendering, and toned down excessive planet bloom.

## 0.14.3 - 2026-03-01
- Refined sky-disc layout rules for wide windows.

## 0.14.2 - 2026-03-01
- Added **altitude keys**, a **nadir marker**, **continuous star flare**, **below-horizon celestial rendering**, and **planet gauge markers**.
- Packaging: Added `pytest` as a development dependency.

## 0.12.4 - 2026-02-28
- Prepared PyPI-facing `pyproject` metadata and updated install guidance to use `pipx install zstarview`.

## 0.12.3 - 2026-02-28
- Introduced the **satellite-based cloud-disc overlay** as an experimental feature.
- Also added `--enlarge-moon`, `--sky-opacity`, and black-theme tuning.
