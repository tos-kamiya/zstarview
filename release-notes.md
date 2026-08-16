# zstarview Release Notes

## 1.54.7 - 2026-08-17

- Smoothed star motion between sky-data updates with a local mesh transform,
  keeping bright stars, asterisms, labels, and twinkle effects aligned.
- Softened the outer edges of default halftone cloud dots with an
  approximately 2px alpha fade, without adding a separate outline.

## 1.54.6 - 2026-08-15

- Stopped artificial-satellite fetch work promptly when the application shuts
  down, including cancellation of in-progress Horizons sequences.

## 1.54.5 - 2026-08-15

- Removed Nominatim search terms, place names, and exact coordinates from
  operational logs while preserving them in the cache and saved location data.

## 1.54.4 - 2026-08-15

- Reported Overture Maps catalog and download failures as concise application
  warnings without exposing the subprocess traceback.
- Added a preflight check for PLATEAU temporary-filesystem space and simplified
  PLATEAU cache listing to use `--list` with external text filtering.

## 1.54.3 - 2026-08-15

- Kept Moon and Sun hover images in front of their cross markers for
  consistent hover rendering.

## 1.54.2 - 2026-08-15

- Keep the Sun hover image visible with a subdued solar-marker color when
  below-horizon atmospheric transmission would otherwise make it black.

## 1.54.1 - 2026-08-15

- Replaced the Sun hover image with time-addressable SDO/HMI Continuum
  imagery so historical display times show sunspots.

## 1.54.0 - 2026-08-15

- **Added asynchronous SDO/AIA 193 EUV imagery to the 5x Sun hover display,
  including observation-time labeling, caching, and stale-cache fallback.**
- Aligned the north-up solar imagery with the projected solar north direction
  for the observer and active screen center.
- Made the black background outside the solar disc transparent while retaining
  faint coronal emission, with a smooth transition across the solar limb.

## 1.53.0 - 2026-08-15

- **Added asynchronous NASA Dial-A-Moon imagery to the 5x Moon hover display,
  with disk masking, caching, and the existing rendered phase as a fallback.**
- Aligned NASA North Up Moon images with celestial north for the observer and
  active screen center.
- Refined the enlarged rendered Moon outline so its illuminated and dark limbs
  remain distinguishable across lunar phases.
- Fixed startup rendering getting stuck when the window was resized while the
  startup log overlay was visible; the resize is now applied after initial
  data loading completes.

## 1.52.6 - 2026-08-15

- Swapped the `A` and `P` GUI shortcuts so `A` toggles Aircraft, matching the
  CLI `-a` option, and `P` toggles Asterisms.

## 1.52.5 - 2026-08-15

- Synchronized dynamic-layer redraws to even calendar seconds, eliminating
  extra repaint bursts from independently updated planets, aircraft,
  satellites, and twinkle.
- Precomputed planet positions for the next display tick and stopped applying
  the fixed-star interpolation transform to directly calculated planets.
- Aligned the regular scheduler to one-second calendar boundaries and changed
  fast-mode recovery to an independent 1000 ms idle timer.

## 1.52.4 - 2026-08-15

- Added a session-only `Twinkle` item to the Layers menu, allowing star
  twinkle to be stopped and restarted without changing the startup candidate
  count; the item is unavailable when the configured count is `0`.
- Aligned the Layers menu, startup settings dialog, CLI help groups, and
  English/Japanese reference documentation around the Celestial, Atmosphere,
  Ground, and General categories.

## 1.52.3 - 2026-08-15

- Kept star-surface caches reusable across twinkle updates and batched twinkle
  target lookup, reducing representative 30-target lookup time from about
  `20.6 ms` to `0.8 ms` while keeping twinkle out of fast-mode rendering.
- Hid model-forecast precipitation lines during fast-mode viewport rendering.
- Lowered the default GMN meteor trail opacity from `0.7` to `0.5` and reduced
  the maximum trail body half-width from `2.0 px` to `1.75 px`.

## 1.52.2 - 2026-08-15

- Raised the default star twinkle candidate count from `10` to `30` and
  renamed the option to `--twinkle-count`.
- Renamed the internal twinkle implementation and aligned the CLI option
  groups and README reference files with the startup settings dialog.

## 1.52.1 - 2026-08-15

- **Added display-only star scintillation to the normal GUI, with configurable
  per-update candidate count via `--twinkle-count` (default: 10).**
- Kept scintillation disabled for Atlas and `zstarview-export-image`.
- Restyled GMN meteor trails as single tapered diamond-shaped marks with matching
  body and age-label colors; Atlas uses a near-black gray body for its white
  background.
- Raised the default GMN meteor trail opacity from `0.4` to `0.7` and increased
  the trail body width for better visibility.
- Restyled precipitation marks as dark-blue solid lines with thinner,
  light-blue square-ended dashes for better visibility across display sizes.

## 1.51.5 - 2026-08-14

- Added meteor-trail screenshots and descriptions to the English, Japanese, and
  PyPI README files, including context about the Perseid meteor shower.

## 1.51.4 - 2026-08-13

- Report the age range of the meteor trails actually displayed after the
  geographic filter and display-count limit.
- Raised the default GMN meteor trail display limit from 100 to 200 trails.

## 1.51.3 - 2026-08-13

- Fixed GMN meteor cache refreshes to use each file's `fetched_at_utc`
  six-hour TTL regardless of the date encoded in the GMN filename.

## 1.51.2 - 2026-08-13

- Moved `Inverted City` activation from the Layers menu into the `Space`
  display-mode cycle, skipping it when Urban Outline is unavailable.
- Added an `Inverted City [Space]` status indicator for the active mode.

## 1.51.1 - 2026-08-13

- Corrected `Inverted City` for elevated viewpoints so buildings already below
  the horizon remain below it instead of being reflected into the sky.

## 1.51.0 - 2026-08-13

- **Added the `Inverted City` GUI layer mode, which reflects building outlines,
  roof fills, and distant skyscraper outlines below the geometric horizon.**
- Kept `Inverted City` as a temporary GUI display setting; it is not persisted
  and is not available in `zstarview-export-image`.
- Kept all inverted-city building points at or below `alt=0°`, including when
  the observer is above the surrounding buildings, such as on a tower.

## 1.50.1 - 2026-08-13

- Lightened the GMN meteor trail and age-label color to a pale green while
  keeping the building-outline color unchanged.

## 1.50.0 - 2026-08-13

- Changed GMN meteor trails to use the latest available 24-hour observation
  window at or before the display time, and added its relative age range to
  the status line.
- Changed GMN meteor trails to retain and display the observation-time Alt/Az
  direction instead of following the fixed-celestial position as the display
  time changes.
- Limited the default GMN meteor trail display to the newest 100 observations,
  dropping older trails when a 24-hour window contains more records.
- Added `--meteor-trails-max-candidates` to adjust that display limit; `0`
  disables the limit.
- Added compact observation-age labels such as `-32h` at meteor trail starts;
  these labels intentionally do not participate in collision avoidance.
- Removed horizon-based filtering from GMN meteor trails; the geographic
  candidate radius is now `500 km` for a less horizon-heavy display.

## 1.49.1 - 2026-08-13

- Limited Open-Meteo forecast precipitation to real-time views, matching the
  existing cloud and tropical-cyclone time gating.
- Hid GMN meteor trails during fast-mode viewport rendering to keep interaction
  lightweight, while preserving the configured layer state.

## 1.49.0 - 2026-08-13

- Added Global Meteor Network observation trails, including cached daily-data
  loading, geographic and horizon filtering, fixed-celestial projection, and
  age-based fading over the previous 24 hours.
- Added the `--meteor-trails-opacity` option with a default of `0.4` and a
  Celestial-layer menu toggle. Setting the option to `0` disables GMN fetching
  and rendering for that run.
- Added GMN attribution and CC BY 4.0 licensing information to the bundled
  license and data-source summary.

## 1.48.5 - 2026-08-12

- Kept AKARI IR bands at the configured maximum opacity in simplified view,
  without the normal Sun-altitude fade.

## 1.48.4 - 2026-08-12

- Limited Road Lights processing to the nearest configured road candidates,
  added `--road-light-max-candidates`, and raised the default road-light
  opacity to `0.12` after visual review.

## 1.48.3 - 2026-08-12

- Clarified README screenshot descriptions and added common CLI startup
  examples for place names, view centers, aircraft, and forecast precipitation.

## 1.48.2 - 2026-08-12

- Added the observer's forecast precipitation to the center of the viewport as
  a distinct, bounded rain marker, while retaining the surrounding 48-point
  forecast display.
- Extended the precipitation streak intensity scale so heavy forecast rain
  retains useful visual differences up to 62 mm/h, while light rain remains
  visible with a single streak.
- Clarified that the Open-Meteo precipitation layer is a model forecast and may
  differ from the precipitation that actually occurs.

## 1.48.1 - 2026-08-12

- Split the viewer's bottom status display into fixed and dynamic rows, ordered
  the dynamic row as cloud, precipitation, tropical cyclone, satellite, and
  aircraft, shortened the Open-Meteo forecast label, and moved bottom-positioned
  observation info above the status rows to prevent overlap.

## 1.48.0 - 2026-08-12

- Added an opt-in Open-Meteo forecast precipitation overlay to the interactive
  viewer and image exporter. Rain is rendered as distance-faded blue streaks,
  with versioned confirmation of the Free API's non-commercial terms,
  Open-Meteo attribution, and the `-P` / `--precipitation-opacity` option.
- Reorganized `zstarview-gui` into `Observing Conditions`, `Celestial`,
  `Atmosphere`, `Ground`, and `General` tabs. Location, view, and time now use
  collapsible sections, and the City field is constrained to about three lines.
- Added a selectable, copyable `Licenses and Data Sources` dialog containing
  version-matched offline information, while retaining external links to the
  providers' current terms.

## 1.47.0 - 2026-08-11

- Softened the normal-viewer urban-outline blue-white from `(200, 225, 255)`
  to `(207, 229, 255)`, halfway back toward the previous palette.
- Reused fresh normal 2.5km or dedicated 0.15km Overture building caches for
  `--use-building-top`, downloading only missing small-cache feature types.
- Added optional display-only black/white tone compensation to `zstarview`,
  with shared numbered calibration patches in the CLI and GUI launcher,
  cached final-frame conversion, and uncorrected screenshots and volatile UI.
  Kept the device-specific options out of `zstarview-export-image` so exported
  images remain portable.
- Discarded scheduled image-export tasks that are at least three minutes late,
  with a warning that identifies the task, scheduled time, and delay instead
  of running overdue work as catch-up output.

## 1.46.6 - 2026-08-10

- Replaced local-clock visual profiles with Sun-altitude transitions for night
  lights, AKARI IR bands, the sky ambient floor, closed urban-outline roof
  fills, and Road Lights. Urban roof fills and continuous road light now also
  decrease after the Sun's minimum altitude, while outline and lamp-point
  visibility remain independent.
- Increased closed urban-outline roof-fill brightness by 20% without changing
  its Sun-altitude or post-minimum activity curves.
- Kept a subtle 5% floor for closed urban-outline roof fills during daylight.

## 1.46.5 - 2026-08-10

- Fixed fresh Road Lights API responses so successfully downloaded road data
  are returned and rendered immediately.
- Placed Road Lights against Copernicus DEM terrain in both the GUI and image
  exporter, using `DEM + 0.5m` for vehicle-light strokes and `DEM + 8m` for
  lamp points.
- Shifted the normal `zstarview` urban-outline color to a cooler LED-like
  blue-white `(200, 225, 255)` while retaining the Atlas gray palette.
- Reduced the density of orange Road Lights lamp points by increasing their
  spacing from approximately 120m to 240m.
- Softened sunset sky colors by restoring the atmospheric model's green-channel
  balance while retaining its slight red reduction.

## 1.46.4 - 2026-08-10

- Added 30-day Road Lights cache expiry, 5 km fallback after 10 km Overpass
  failures, and stale-cache reuse when refreshes fail.
- Updated Road Lights menu and startup-dialog regression tests.

## 1.46.3 - 2026-08-10

- Added the `-B` short option for `--use-building-top` to the main and
  `zstarview-export-image` command-line interfaces.
- Improved fast-mode interaction rendering by omitting water-boundary lines
  while retaining water points, and by preventing duplicate planet rendering.

## 1.46.2 - 2026-08-09

- Added Road Lights rendering to `zstarview-export-image`, reusing the GUI
  cache, geometry preparation, and rendering behavior. Road-light fetch or
  projection failures can be skipped when partial export output is allowed.

## 1.46.1 - 2026-08-09

- **Added the `Road Lights` overlay, which derives a stylized continuous road
  glow and approximately 120m-spaced lamp points from cached OSM major-road
  centerlines within the 0.5km–10km annulus. Added the
  `--road-light-opacity` setting, startup-dialog control, Layers-menu toggle,
  background loading, and cache/API status reporting.**
- Tuned Road Lights preparation and rendering: shared coordinate transforms,
  non-blocking shutdown, distance attenuation for distant roads, and a smooth
  lamp-point twilight fade from Sun altitude `0°` to `-4°`.
- Restored Sun, Moon, and planet cross gauges in `zstarview-export-image`
  output.

## 1.46.0 - 2026-08-09

- Moved Sun, Moon, and planet markers out of the base/present frame caches so
  the small solar-system marker layer is redrawn dynamically.
- Prevented the normal-size Moon marker from remaining underneath the hovered
  5x Moon, and kept the gauge cross visible for both enlarged display paths.
- Added a faint, thin white lunar limb around enlarged Moon rendering to make
  the new-Moon phase easier to see.

## 1.45.4 - 2026-08-09

- Added independent cloud-altitude shell fields so sunset tinting can be
  composited separately for the 3 km, 5 km, and 7 km cloud layers.
- Used the modelled solar-horizon sky colour for shell-specific sunset cloud
  tinting and reduced excessive green in low sunset sky colours.
- Updated sky-disc overlay tests for the antialiased circular clipping.

## 1.45.3 - 2026-08-07

- Added a phase-aware lunar outline for ordinary `--bright-bodies outline`
  rendering, drawing only the illuminated outer limb and terminator arc.
- Kept the enlarged and hover Moon rendering as a detailed filled phase image
  and prevented duplicate hover crosshairs.

## 1.45.2 - 2026-08-07

- Refreshed the sky-color disc independently near twilight while keeping star
  updates on the normal sky-update interval.
- Reduced sky-disc rendering cost and memory use with lower internal
  resolution, fewer solar-path samples, and chunked NumPy processing.
- Smoothed the final sky-disc edge with an antialiased circular clip.

## 1.45.1 - 2026-08-06

- Fixed the RGB scattering coefficients used by the sky model to the
  representative-wavelength formulation documented in the rendering design.

## 1.45.0 - 2026-08-06

- **Bundled a compact global CAMS EAC4 monthly AOD550 climatology for
  regional and seasonal aerosol effects without runtime API credentials or
  network access.**
- **Applied the bundled aerosol values to the Rayleigh/Mie sky-scattering
  pipeline** and used a 0-to-10-degree solar-horizon color average for the
  time-of-day marker.
- Added CAMS source, licence, attribution, and modified-data notices to the
  README files and bundled asset metadata.

## 1.44.1 - 2026-08-05

- Reduced the PyPI source distribution by excluding repository-only `docs/`
  content from the sdist while retaining the package README files.

## 1.44.0 - 2026-08-05

- **Made the OpenSky aircraft overlay opt-in: it is disabled by default and
  requires an explicit positive `-a` / `--aircraft-opacity` value to enable
  requests and rendering.**
- Documented the OpenSky request scope, conservative five-minute refresh
  interval, short-lived local caching, and identifying HTTP User-Agent in the
  English and Japanese README files.

## 1.43.4 - 2026-08-05

- Tuned the AKARI infrared display for the refreshed screenshot examples,
  including the creative-Hubble palette, gamma and highlight compression.
- Moved AKARI display tuning constants into a shared constants module and
  reduced the night-sky ambient background contribution.
- Hardened cached overlay updates and saved aircraft diagnostic frames at
  one-minute intervals.
- Removed the star-pixel Gaussian blur and refreshed the screenshot gallery.
- Added localized README timelapse examples and aligned AKARI band descriptions.

## 1.43.3 - 2026-08-02

- Faded water surfaces hidden behind nearer, higher terrain during normal
  rendering. Kept terrain-dependent water alpha in runtime memory and skipped
  the extra calculation during fast viewport-interaction rendering.

## 1.43.2 - 2026-08-02

- **Improved offline handling for satellite, aircraft, tropical cyclone, DEM,
  and water overlays.** Reused valid caches where possible, hid unavailable
  layers without traceback-style UI output, and retained partial DEM caches
  when adjacent tiles were unavailable.

## 1.43.1 - 2026-07-31

- Applied the local-time artificial-light activity curve to AKARI IR bands,
  with a default effective opacity range of `0.05` to `0.10`.

## 1.43.0 - 2026-07-31

- **Refreshed the bundled GeoNames and IAU star-name snapshots.**
- **Regenerated the bundled star catalogs from Hipparcos I/311 with the
  optional Tycho-2 I/259 supplement;** retained the legacy I/239 input for
  catalog-generation compatibility.

## 1.42.0 - 2026-07-30

- **Added `--temp-dir` to `zstarview-download-plateau-buildings`** so large
  CityGML downloads and extraction can use a filesystem with sufficient free
  space instead of the system temporary directory.
- Added actionable error guidance for temporary-filesystem exhaustion and
  documented the recovery command in the English and Japanese README
  troubleshooting sections.

## 1.41.3 - 2026-07-30

- **Added a smooth local-time deep-night profile for artificial-light layers,
  reducing night-light intensity around 02:00–04:00 while keeping edge/ridge
  glow and building outlines independently readable.**
- Made closed roof-polygon fills follow the night-light reduction; the fills
  may disappear in deep night while building outlines remain visible.
- Added time-dependent AKARI IR bands opacity, with a default effective range
  of `0.10` during ordinary night to `0.15` around 02:00–04:00. The
  `--akari-ir-bands-opacity` value now represents the deep-night maximum.
- Added a matching sky-disc ambient floor transition from RGB `[1, 2, 5]`
  to `[2, 4, 10]` and tuned the AKARI display gamma to `0.35`.
- Added regression coverage for the time-dependent sky, night-light, AKARI,
  and roof-fill behavior.

## 1.41.2 - 2026-07-30

- Changed the time-of-day marker's solar-direction sky-color sample altitude
  from 20 degrees to 0 degrees.

## 1.41.1 - 2026-07-29

- Clarified transparent-window border behavior so the resize grip remains
  owned and rendered by the dedicated resize-grip widget.

## 1.41.0 - 2026-07-29

- **Replaced empirical daytime sky-color tuning with spherical-Earth
  Rayleigh and aerosol Mie scattering, including a broader sunset
  transition.**
- **Added the creative Hubble palette for the molecular-cloud overlay** and
  documented the AKARI source attribution.
- Moved the exported time-of-day marker to the upper-left corner.
- Fixed the duplicate resize-grip line in the lower-right corner.

## 1.40.2 - 2026-07-29

- Switched the AKARI IR bands overlay to the JWST-inspired display palette.
- Increased the AKARI overlay gamma lift to `0.5` while retaining the HSV V-value compression.

## 1.40.1 - 2026-07-28

- Reduced AKARI IR bands rendering work by generating the overlay at one
  quarter of the final width and height, then applying smooth upscaling for
  the final image and `zstarview-export-image` output.

## 1.40.0 - 2026-07-28

- **Added a reversible AKARI infrared color palette switch, with the AKARI
  blue/green/red mapping as the default and a JWST-inspired palette available
  for comparison.**
- **Raised the default `--akari-ir-bands-opacity` from `0.1` to `0.15`.**

## 1.39.3 - 2026-07-28

- Extended the sunset transition from `-2°..2°` to `-1°..4°` and added a
  restrained yellow tint to the solar-direction glow as the Sun rises through
  `4°..12°`.
- Blended the solar glow as 70% additive light and 30% target-color replacement
  so sunset hues remain visible at stronger sky opacity settings.
- Made the weak all-direction horizon haze shift toward the sunset RGB as the
  Sun descends, while keeping its base strength unchanged.

## 1.39.2 - 2026-07-28

- Tuned the solar glow toward a whiter additive light and strengthened the
  sunset color near the horizon.

## 1.39.1 - 2026-07-28

- Added the solar-direction sky-color marker to the regular viewer, including
  twilight color handling and HUD-aware corner placement.
- Removed the decorative outer border from the frameless window chrome while
  retaining the menu control.

## 1.39.0 - 2026-07-26

- **Added `--list-schedule` (with `--list` as an alias) to the export-image
  schedule runner. It lists the next occurrence of each task in time order and
  exits without starting the scheduler.**
- **Retired the urban-outline feature-type and minimum-building-height tuning
  controls from the documented CLI and the initial GUI settings dialog.
  Existing CLI invocations remain compatible but now warn and use the active
  defaults.** Use `--urban-outline-max-candidates` for outline performance
  tuning.

## 1.38.2 - 2026-07-26

- Set the minimum sampling distance for inland-water dots, including rivers and
  lakes, to 5m.

## 1.38.1 - 2026-07-26

- Extended `zstarview-download-coastline --all` to download the optional global
  25m water-mask ZIP in addition to all coastline columns.
- Added `zstarview-download-coastline --water-25m` for explicit 25m water-mask
  downloads, with manifest and SHA-256 verification and a versioned cache.
- **Added optional 25m water-mask sampling for the nearest 250m, with a 50m
  minimum sample distance and fallback to the bundled 125m tiles.**

## 1.38.0 - 2026-07-26

- Fixed coastline cache replacement after an interrupted download leaves a non-empty column directory.
- Prevented the coastline reader from using a cache without a completed `READY` marker.
- Added per-column download progress and extraction status messages to `zstarview-download-coastline`.
- **Added `zstarview-download-coastline` for selecting and downloading
  coastline vector data by longitude range from the GitHub Release cache.**
- Added manifest, SHA-256, safe ZIP path, and atomic per-column cache installation checks.

## 1.37.6 - 2026-07-26

- Refreshed the bundled 125m, 250m, and 500m water-mask tiles from the updated OSM Water Polygons source.
- Added preparation tooling for the separately distributed coastline vector Release assets.

## 1.37.5 - 2026-07-25

- Fixed night-light altitude placement at elevated observation sites by using the observer's absolute elevation instead of treating the added observer height as an elevation above sea level.
- Raised the default `--night-light-opacity` from `0.1` to `0.16` so the overlay is easier to distinguish at startup.
- Added the `-s` short option for configuring sky opacity.
- Refined sky-color gradients and added warmer low-horizon haze for more natural dawn and dusk views.
- Preserved Geo-satellite cloud overlays when the viewing direction changes.
- Improved night-light rendering at elevated observation sites and added a night-light boundary scanner.
- Covered the sky-disc field-of-view edge to avoid visible gaps.

## 1.37.4 - 2026-07-20
- Bumped the package version to `1.37.4`.

## 1.37.3 - 2026-07-20
- Added read-only `--list` support to `zstarview-download-plateau-buildings` for listing valid PLATEAU caches by municipality code, dataset year, and saved path.
- Added `--jsonl` output for detailed cache metadata and `--city-code` filtering in list mode.
- Documented the cache listing workflow in the English and Japanese READMEs.

## 1.37.2 - 2026-07-19
- Disabled star-position interpolation when `--sky-update-interval` exceeds 90 seconds.
- Aligned the interpolation midpoint and range with the configured update interval, so 90-second updates use a midpoint at 45 seconds and a `+/-45` second range.

## 1.37.1 - 2026-07-19
- Removed unused production helpers and obsolete internal cache code identified during the production-function audit.

## 1.37.0 - 2026-07-19
- Set the regular `zstarview` projected edge FOV default to `90` degrees and the shared content FOV default to `115` degrees.
- Kept Atlas at its dedicated `90`-degree edge FOV and `110`-degree content FOV defaults.

## 1.36.6 - 2026-07-18
- Darkened the below-horizon ground reset color to `(18, 18, 18)` consistently across themes while preserving theme-specific opacity.

## 1.36.5 - 2026-07-18
- Corrected `zstarview-atlas --help` to show the Atlas-specific Vmag default (`4.0`) and catalog limit (`6.0`), while keeping the regular `zstarview` help unchanged.

## 1.36.4 - 2026-07-18
- Added Atlas bright-star diamond markers with subdued dark contrast and tuned colored-marker opacity for clearer white-background presentation.
- Added dark outlines around Atlas planet discs, changed the Atlas Sun marker to the same full-size cross gauge as `zstarview`, and kept the Moon phase rendering unchanged.
- Improved `zstarview-download-plateau-buildings` handling when a requested batch has no available building tiles, and clarified PLATEAU data availability and disclosure documentation.
- Updated the ClickPy download badges and links in the English and Japanese README files.

## 1.36.3 - 2026-07-17
- Tuned the sky-color disc toward a cooler, bluer appearance at high Sun altitudes and a stronger red-orange appearance near the horizon.
- Reduced haze blending, softened the sunset tint in the exact solar direction, and documented the sky-disc color factors and processing order.

## 1.36.2 - 2026-07-16
- Improved `zstarview-download-plateau-buildings` handling of HTTP 404 catalog responses so the CLI reports that PLATEAU building data may not be available for the municipality instead of showing a traceback.

## 1.36.1 - 2026-07-16
- Added the urban-outline source and final outline count to `zstarview-export-image` PNG metadata under `extra.urban_outline`.

## 1.36.0 - 2026-07-16
- **Added optional PLATEAU building-cache preparation and runtime source
  selection, with Overture Maps fallback when no completed PLATEAU cache
  covers the observation area.**
- Extended `zstarview-download-plateau-buildings` to accept municipality code ranges and comma-separated code lists.
- Added catalog-based PLATEAU cache update detection when the preparation CLI is run.
- **Added PLATEAU LOD1 building roof outlines and experimental LOD2
  `RoofSurface` outlines for more recognizable three-dimensional
  silhouettes.**
- Merged shared edges from nearby-height LOD2 roof surfaces and increased the merging tolerance with distance, reaching 3m at 3km and beyond.
- Added PLATEAU data-source labels to the urban-outline status display.

## 1.35.3 - 2026-07-14
- Restored `width,50,0.85` as the default cloud-stripe style for the regular viewer and CLI.
- Raised the default cloud opacity from `0.06` to `0.10` for clearer visibility across displays.
- Added the optional `width-quantized` cloud-stripe mode with five width levels, gaps at width transitions, and round line caps.

## 1.35.2 - 2026-07-14
- Improved Atlas compass-direction labels with a restrained gray color and medium font weight so they remain readable near the horizon line.
- Improved `zstarview-atlas` text readability with white outlines on labels, the fixed HUD, and the status line.
- Added a sun-altitude-based blue/orange/dark-gray right-isosceles triangle marker at the HUD-free left corner in `zstarview-atlas`, with an orange horizon transition for sunset-like readability.
- Thickened the artificial-satellite cross cursor in `zstarview-atlas` to match the visual weight of the aircraft overlay lines.
- Made Atlas simplified views hide the ground fill and urban outlines, and added named-star labels to the labels-enabled simplified view.
- Clarified the simplified-view status text as `Simplified: no labels [Space]` and `Simplified: with labels [Space]`.
- Split the Atlas and regular zstarview render pipelines while keeping shared rendering primitives in the common pipeline.
- Restored `width,50,0.85` as the default cloud-stripe style so cloud visibility does not depend as strongly on display contrast.
- Raised the default cloud opacity from `0.06` to `0.10` for clearer visibility across displays.

## 1.35.1 - 2026-07-12
- Raised the default cloud, night-light, and ridge-glow opacities slightly so these faint overlays remain visible on displays with limited color reproduction.

## 1.35.0 - 2026-07-11
- **Adopted the 2025 annual VIIRS Nighttime Lights VNL v2.2 product from the
  Earth Observation Group (EOG) for the night-light overlay.**
- **Documented distribution of the converted GeoTIFF through GitHub Releases
  instead of bundling the large raster data in the PyPI package.**
- Added the EOG attribution, change notice, source citation, and CC BY 4.0 licensing requirements for the redistributed night-light product.
- Added manifest and per-tile SHA-256 verification, atomic cache installation, and lazy download of only the tiles needed around the observation area.

## 1.34.0 - 2026-07-11
- **Added `zstarview-atlas`, a parchment-inspired white-background GUI entry
  for clouds, bright stars, planets, aircraft, satellites, and geographic
  helper layers.**
- **Added the Atlas presentation with black labels, sky-color rendering
  disabled by default, blue-gray cloud tinting on white, and a light-background
  star path** that skips stars below 1 px apparent diameter while drawing star
  outlines before colored bodies.
- Added Atlas-specific rendering for DSO and asterism outlines, Earth guides, cloud stripes, observation HUD information, and opaque startup-log backgrounds for improved readability on the white map surface.
- Added a shared OpenSky aircraft fetch lock and rate-limit marker so multiple GUI instances do not multiply aircraft requests across different observation areas.
- Kept `zstarview-export-image` outside the GUI shared rate-limit skip so explicit single-shot captures can still fetch aircraft data, while respecting the shared fetch lock.
- Refactored night-light profile inputs and processing helpers without changing the existing user-facing layer controls.

## 1.33.0 - 2026-07-05
- **Changed tall-window sky-disc geometry so portrait layouts keep the disc
  centered, allow the edge FOV to extend beyond the left/right window edges,
  and blend toward a content-FOV height fit only by `1:2` and taller aspect
  ratios.**
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
- **Added a halftone cloud-stripe mode and made `--cloud-stripe halftone` the
  default** so cloud overlays now render as quantized halftone circles and
  chains unless overridden.

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
- **Added a dedicated `--ridge-glow-opacity` control and kept ridge glow
  wired through the GUI, export, and startup render paths as a separate layer
  from night light.**
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
- **Added the initial tropical cyclone overlay for displaying nearby storm
  positions and projected paths.**
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
