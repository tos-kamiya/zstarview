# zstarview Release Notes

## 1.33.0 - 2026-07-05
- Changed tall-window sky-disc geometry so portrait layouts keep the disc centered, allow the edge FOV to extend beyond the left/right window edges, and blend toward a content-FOV height fit only by `1:2` and taller aspect ratios.
- Kept GUI rendering, export-image rendering, worker payload validation, and hover lookup on the same edge/content FOV geometry.
- Normalized `content_fov_deg` so it never falls below `edge_fov_deg`, preserving the projection scale invariant across direct viewer/projection construction paths.

## 1.32.13 - 2026-07-04
- Kept aircraft labels hidden in fast mode so the rapid refresh path stays uncluttered during interaction.
- Tuned the halftone cloud overlay to use a slightly larger dot base, proportional grid spacing, and an expanded content FOV so edge clouds stay visible a bit more consistently.

## 1.32.12 - 2026-07-03
- Reworked aircraft overlays into local-plane ribbon polygons, with a fallback polyline when the ribbon collapses, so nearby aircraft tracks read more like swept shapes than thick lines.
- Simplified the aircraft ribbon stroke width to a fixed thin outline, keeping the filled ribbon visually dominant.
- Kept the aircraft and palette sample colors aligned with the current runtime accent set while the aircraft presentation was being tuned.

## 1.32.11 - 2026-07-02
- Added `zstarview-install-overturemaps-exe-cli` as a copy-only staging helper for Windows `overturemaps` release executables, and taught the Overture import path to prefer a staged cache executable before falling back to `PATH`.
- Clarified the Arm64 Windows installation flow in the English, Japanese, and PyPI-facing READMEs so users can stage a Windows x64 `overturemaps` executable when Arm64 wheels are unavailable.
- Reduced wait-log frequency in `zstarview-export-image-schedule-runner` so it logs once when it starts waiting and then every 15 minutes instead of every minute.

## 1.32.10 - 2026-06-30
- Fixed the GOES CMI loading path so the worker no longer crashes while opening downloaded GOES files with CF time decoding enabled.
- Added `zstarview-diagnose-cloud-source` as a worker-oriented diagnostic CLI for cloud-source failures, with isolated output-dir handling, structured diagnostics, and `--source-file` support for already-downloaded GOES files.

## 1.32.9 - 2026-06-29
- Reworked internal GUI render-state transport so the live view and off-screen image rendering share a clearer render-input flow, with no intended visual behavior change.

## 1.32.8 - 2026-06-29
- Kept aircraft callsign labels visible in the GUI fast interaction frame so the live view no longer drops them while the window is refreshing.
- Saved the `aircraft-ready` debug screenshot from the final composited frame so the capture now includes the same overlays and labels shown on screen.

## 1.32.7 - 2026-06-29
- Added debug-level GOES cloud-download tracing so intermittent worker failures can be isolated more easily without changing the normal status display.

## 1.32.6 - 2026-06-27
- Restored the default `--night-light-opacity` value to `0.04` and `--ridge-glow-opacity` to `0.03`, with the CLI, GUI, render defaults, tests, and docs aligned to the runtime behavior.

## 1.32.5 - 2026-06-27
- Raised the default `--ridge-glow-opacity` value to `0.06` and aligned the CLI, GUI, render defaults, and documentation with the new default.

## 1.32.4 - 2026-06-27
- Added the `Square Window` menu toggle so the client area can be corrected to a square using the shorter side after resize interactions.
- Kept the resize interaction responsive by deferring the square correction until the fast-mode to normal-render transition, and suppressing the correction while a resize drag is still in progress.
- Updated the menu wording and hierarchy in the GUI docs to match the current File / Search / Layers / View Direction layout.

## 1.32.3 - 2026-06-27
- Simplified the GUI redraw flow so cloud source-ready no longer triggers an immediate repaint, which keeps the visible update tied to the later projection-ready step instead.
- Applied the same silent source-ready policy to the Geo-satellite path, and allowed it to reuse a cached raw image from a previous launch before rebuilding the visible overlay.

## 1.32.2 - 2026-06-27
- Aligned the asterism and DSO hover stroke widths at `2.2` so both hover outlines read a bit lighter and match each other more closely.

## 1.32.1 - 2026-06-25
- Tightened the halftone cloud overlay so the grid spacing now has a 20 px minimum, keeping small windows from turning the cloud field into a dense blob.
- Kept the halftone outline pass removed, so the cloud rendering stays cleaner while preserving the cloud positions.
- Aligned the GUI cloud-stripe startup default with the CLI default so both now start in `halftone,40,1.7`.

## 1.32.0 - 2026-06-25
- Added PNG metadata to `zstarview-export-image` output so exported PNGs now carry the app version and HUD-related information, with optional `--place` and `--search` resolution details available to tools such as `exiftool`.
- Kept the export-image metadata format versioned as `zstarview.export-image-metadata.v1` so future readers can distinguish compatible payloads from later schema changes.
- Added a halftone cloud-stripe mode and made `--cloud-stripe halftone` the default so cloud overlays now render as quantized halftone circles and chains unless overridden.

## 1.31.18 - 2026-06-25
- Added a 120-second default timeout for Overture Maps building downloads used by the urban outline layer. It is exposed as `--urban-outline-download-timeout-seconds` on the `zstarview`, `zstarview-gui`, `zstarview-debug`, and `zstarview-export-image` CLIs, and as `download_timeout_s` in `UrbanOutlineController`. This prevents the GUI/CLI from waiting indefinitely when the `overturemaps` subprocess stalls.

## 1.31.17 - 2026-06-24
- Fixed `zstarview-export-image` so aircraft and satellite overlays now render in exported images; the internal `time_obj` field was left unset, causing both `draw_aircraft_overlay` and `draw_satellite_overlay` to silently skip rendering.
- Changed the frameless-window "Fit to Screen" action to scale the client area to 90% of the available screen geometry instead of 100%, leaving a small margin around the window while still preserving the aspect ratio.

## 1.31.16 - 2026-06-21
- Made the night-light and ridge-glow spread shrink with distance so nearby glow stays broad while far glow tightens more gradually.
- Increased the altitude sampling resolution for night-light profile generation and added draw-time cropping of inactive altitude rows in the compositor to keep fullscreen glow rendering faster.

## 1.31.15 - 2026-06-21
- Raised the default opacity for `--night-light-opacity` to `0.06` and
  `--ridge-glow-opacity` to `0.03`, keeping the CLI defaults aligned with the
  GUI startup values.
- Kept night-light opacity flowing through the GUI startup and export paths so
  the renderer now treats a fully disabled night-light layer as disabled during
  startup work as well.

## 1.31.14 - 2026-06-20
- Added a dedicated `--ridge-glow-opacity` control and kept ridge glow wired through the GUI, export, and startup render paths as a separate layer from night light.
- Tuned the glow defaults to `--night-light-opacity 0.04` and `--ridge-glow-opacity 0.02`, and kept the shared glow warmup path enabled when either visible layer is active.
- Hid the tropical cyclone overlay for time-shifted views so live cloud, aircraft, satellite, and cyclone overlays now follow the same real-time gate.

## 1.31.12 - 2026-06-20
- Reduced the cost of night-light and ridge-glow profile generation by precomputing terrain visibility thresholds and reusing the azimuth glow kernel across accumulation passes.
- Kept the compositor behavior unchanged while cutting the profile-generation path from about one minute to about six seconds in the reported case.

## 1.31.11 - 2026-06-20
- Removed the background mouse-press simplified-view override, so simplified view now switches only through `Space`.
- Kept the three-state simplified view cycle in place, including the label-free and label-enabled modes.

## 1.31.10 - 2026-06-20
- Separated the night-light and ridge-glow opacity controls so the base night-light layer and the additive ridge glow can be tuned independently.
- Cleaned up the compositor wiring around the glow layers and updated the public docs to match the current split-layer behavior.

## 1.31.9 - 2026-06-17
- Split the night-light and ridge-glow mask cache away from the full-frame compositing cache, so glow reuse is handled independently from unrelated frame invalidations.
- Kept the full composited frame cache in place, preserving the existing reuse behavior for the rest of the scene.

## 1.31.8 - 2026-06-16
- Softened the cloud warm-threshold blend so local brightness-temperature samples only nudge the equatorial reference, reducing overly aggressive cloud attenuation in scenes where the local view is cooler.
- Kept the warm-threshold update shared across the common cloud render path, so the gentler blend applies consistently to current cloud overlays.

## 1.31.7 - 2026-06-16
- Recolored Sun, Moon, and planet labels with body-derived hues blended toward white, keeping the object identity while reducing the visual clutter of the raw marker colors.
- Applied the same white-blended treatment to simplified-view named-star labels shown in the Space-toggled label mode.

## 1.31.6 - 2026-06-16
- Reworked simplified view into a three-state cycle with `normal`, `simplified view (no labels)`, and `simplified view (labels)`, while keeping background mouse press as a temporary return-to-normal gesture from the simplified states.
- Added fixed-position simplified labels for named stars and artificial satellites, and kept the hover-only satellite label path separate so the two label styles do not duplicate each other.

## 1.31.5 - 2026-06-15
- Added a small screen-fixed noise modulation to the GlowMask alpha path so the night-light glow and ridge glow read as slightly rougher instead of perfectly smooth.
- Kept the noise deterministic at screen resolution so the texture stays stable from frame to frame.

## 1.31.4 - 2026-06-14
- Widened the cloud warm-threshold reference from a single equator line to a `+/-5°` equatorial belt, which makes cloud amount estimation more robust and reduces overbright cloud renders in scenes like Uyuni.
- Kept the Himawari tile selection in sync with the wider reference belt so the warm-threshold sampling and source selection stay aligned.

## 1.31.3 - 2026-06-14
- Forced the `overturemaps` child process to run with `PYTHONUTF8=1` so Windows builds are less likely to fail on non-ASCII building properties during urban-outline downloads.

## 1.31.2 - 2026-06-13
- Restored the `--night-light-opacity` default to `0.07`.
- Made the urban outline overlay fill closed roofline loops with a faint interior fill, and suppressed fills for clipped or very small projected fragments.

## 1.31.1 - 2026-06-13
- Fixed a startup crash in the GUI cloud toggle setup by initializing the cloud-disc handle before the toggle-support check runs.

## 1.31.0 - 2026-06-13
- Added a press-pending simplified display mode that hides clouds, night-light glow, Earth guide, secondary ridges, water, and urban outlines while keeping hover labels available.
- Kept the main terrain horizon visible in press-pending mode and reduced its line width to better match the fast interaction view.
- Replaced the earlier timed `star-focus` wording with the new press-state-driven behavior and fixed a Qt drag-exclusion warning in the window drag path.

## 1.30.8 - 2026-06-13
- Raised the default opacity for the `--night-light-opacity` and `--ridge-glow-opacity` CLI options so the night-light glow and ridge glow are easier to see at startup.
- Updated the CLI docs and startup assertions to match the new defaults.

## 1.30.7 - 2026-06-13
- Clarified tropical cyclone empty-data handling so the overlay now shows `TC none` when the service returns an empty collection, while keeping the expected no-observation warning at `WARNING` level.

## 1.30.6 - 2026-06-12
- Cleaned up internal ridge-glow, window-render, and test harness plumbing around unused arguments and dead helpers. No user-visible behavior change.

## 1.30.5 - 2026-06-12
- Treated an empty tropical-cyclone observed-position feed as a normal no-storm state, so the overlay now clears quietly instead of logging a warning.

## 1.30.4 - 2026-06-12
- Made `zstarview-export-image` warn and continue when cloud rendering is unsupported for the selected export region, so Europe-band observers can still export a scene without an unexpected hard stop.
- Kept the Geo-satellite cloud path available for Europe-band exports when `--geo-satellite true` is explicitly enabled.
- Avoided inserting unsupported cloud imagery implicitly when the Geo-satellite flag is omitted.

## 1.30.3 - 2026-06-12
- Blended the sky-glow color from the sky-disc horizon sample and the night-light base tone so the glow stays visible even when the night sky disc is nearly black.
- Kept the closest night-light glow band anchored at ground level so the nearest glow ring remains visible instead of collapsing away.
- Added a sky-glow color override hook so the night-light renderer can accept a custom glow color source without duplicating the band drawing logic.
- Updated the public specification and design notes to reflect the current sky-glow color source and versioned user-agent strings.

## 1.30.2 - 2026-06-09
- Unified outbound `User-Agent` strings so external HTTP requests now identify `zstarview/1.30.2` with a short service suffix.
- Updated the public specification and internal design notes to list the current outbound API identifiers and version-reading rule.

## 1.30.1 - 2026-06-07
- Aligned `zstarview-export-image` water-overlay generation with the GUI path so inland water no longer gets refit against a fresh DEM during export.
- Tightened the export-image terrain-horizon scan density to reduce narrow-sector ridge misses in export output.

## 1.30.0 - 2026-06-06
- Simplified tropical cyclone summary text so multiple storms now show only the count, keeping the overlay label compact.

## 1.29.14 - 2026-06-06
- Inland water polygons without explicit `ele` or `water_level` now use DEM-backed ground height in `zstarview-export-image`, so lakes and rivers follow the sampled terrain height instead of the observer-ground fallback.

## 1.29.13 - 2026-06-06
- Fixed the `--window-frame window` startup path so the standard decorated window now builds its shared `Help` menu correctly.
- Added an `INFO` shutdown log message for user-initiated exits, including the File menu Exit action, the window close button, and the frameless hamburger menu Exit action.

## 1.29.12 - 2026-06-06
- Added six more Horizons-backed spacecraft targets, fixed the search and shortcut labels to use `Spacecraft`, and kept ISS on the separate satellite path.
- Restored the shared GUI worker pool model and raised its concurrency to `2` so long-running background tasks block each other less often.
- Updated the JPL search expectations for the new `sct` lookup pass and cleaned up an unused render import.

## 1.29.11 - 2026-06-05
- Refined the water-surface overlay to use fixed `2.0°` sampling, scan-distance-based opacity, and small unfilled ellipses that flatten faster toward the horizon while staying less visually dominant.
- Refreshed the bundled Europe Geo-satellite gray-common mask with an updated IR frame set and thresholds, improving the packaged mask used by the Geo-satellite cloud path.

## 1.29.10 - 2026-06-04
- Added a `Geo-satellite` startup checkbox in the GUI startup dialog under `General`, matching the CLI option grouping.
- Kept the GUI `Geo-satellite` menu item disabled until the resolved observer location is known, then grey it out only when the observer is outside the Europe band.
- Made the tropical cyclone far marker use the same highlight-width pattern as the satellite crosshair: `1.0` normally and `2.0` when highlighted.

## 1.29.9 - 2026-06-03
- Made the cloud HUD show concrete GOES satellite names such as `GOES G18` and `GOES G19` so partial coverage and source anomalies are easier to distinguish at a glance.
- Kept the existing partial-cloud `?` and failure states unchanged while improving the source label shown in the status line.

## 1.29.8 - 2026-06-03
- Expanded the tropical cyclone overlay to fetch and render every active storm returned by the ArcGIS service, instead of selecting only the latest observed storm.
- Versioned the tropical cyclone cache, discarded stale cache files immediately, and re-fetched fresh data when the schema changes.
- Reduced clutter in the tropical cyclone overlay by showing out-of-range storm-name labels only on hover, while keeping the far marker visible.
- Raised the default tropical cyclone opacity to `0.7` so the overlay stands out more clearly by default.

## 1.29.7 - 2026-06-02
- Refined the water overlay presentation so water markers now use a fixed coarse `4-degree` sampling across all resolutions and a flatter screen-space shape near the horizon.
- Simplified asterism rendering to a single line with a stronger hover alpha, and simplified the tropical cyclone tether to a single thin stroke.

## 1.29.6 - 2026-06-02
- Moved GUI cloud-source fetching to a one-shot subprocess so Himawari and GOES native work no longer runs in the main window process.
- Kept `zstarview-export-image` on the existing one-shot in-process path, since it already runs as a single-shot export command.

## 1.29.5 - 2026-06-02
- Changed the cloud HUD to show `?` for partial Himawari slots so incomplete-but-usable satellite data no longer looks like an internal error.
- Kept the terminal INFO logs detailed, including tile counts and slot selection rationale, so partial-cloud diagnosis still has enough context.

## 1.29.4 - 2026-06-01
- Moved tropical cyclone rendering into the fast overlay path so cyclone markers update consistently with other dynamic overlays.
- Replaced the near-cyclone cylinder with a ground-to-`15km` centerline and added a filled base marker for storms within `400km`.
- Hid aircraft, satellite, and cyclone dynamic overlays during fast-mode resize/viewport interaction frames to avoid stale pre-resize positions.
- Kept fast mode active until refreshed sky data arrives, preventing normal rendering from combining a new view center with old star/DSO subsets.

## 1.29.3 - 2026-05-31
- Simplified the tropical cyclone fallback presentation to a fixed `15km`-tall, `500m`-radius cylinder when wind polygons are unavailable.
- Removed the wind-speed-dependent top cap radius from the fallback cyclone body so the schematic marker stays deterministic.
- Updated the cyclone specification and design notes to match the fixed-cylinder fallback behavior.

## 1.29.2 - 2026-05-31
- Scaled the water-surface point density by render surface size, with `2.0° / 1.0° / 0.5°` tiers and thresholds at `1200px` and `2400px`.
- Kept 4K-class surfaces on the `1.0°` tier so the water dots stay visually separated instead of merging into a band.
- Updated the water-surface specification, design notes, cache scope handling, and export-image path to follow the new density tiers.

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
