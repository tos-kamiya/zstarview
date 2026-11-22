# General

| Option | Description | Default |
| :----- | :---------- | :------ |
| `-h`, `--help` | Show this help message and exit. | |
| `--window-geometry restore\|X,Y,W,H` | Set initial window geometry. Use `restore` to load the last saved position/size, or `X,Y,W,H` to specify explicit integers. Note: on Wayland, window position restore is not available (size restore works). | |
| `--window-frame {frameless,window}` | Choose window decorations. `frameless` keeps the current borderless presentation; `window` uses the platform title bar and frame. | `frameless` |
| `--observation-info auto\|top\|bottom\|off` | Startup mode for the observation-info block. | `bottom` |
| `--initial-display-mode MODE` | Choose the GUI's initial display mode and starting position in the Space display-mode cycle. Values: `normal` (`default`), `inverted-city` (`inverted`), `simple-no-labels` (`simple`), `simple-labels` (`labels`), or `landscape`. `Inverted City` falls back to `normal` when Urban Outline is unavailable. This option is available only for `zstarview`; it is not supported by `zstarview-export-image`. | `normal` |
| `--set-cache-retention-hours HOURS` | Save cache history retention in `config.json` and exit without starting the GUI. Accepts 1 through 720 hours; the saved value applies to later launches. | |
| `--include-direction-grid` | `zstarview-export-image` only. Include the direction grid in exported images, with 30-degree major lines and 10-degree intersection crosses. | |
| `-t`, `--theme {night,transparent,transparent-10..90}` | Theme preset for the window background and sky-disc opacity. `night` is the default dark theme. `transparent` is the `transparent-40` alias; `transparent-10` through `transparent-90` are 10-step transparent presets. | `night` |
| `--visibility-boost MULTIPLIER` | Visibility boost for faint support layers. Values above `1.0` raise opacity for layers such as the terrain horizon, Earth guide, urban outline, sky disc, cloud disc, satellites, aircraft, and asterisms. | `1.0` |
| `--overlay-font-size POINTS` | Base font size for window-drawn labels and HUD text only. This does not affect status line text, menus, dialogs, or standard Qt widgets. | `11` |
| `--offline` | Start without network access. Cached long-lived data (terrain, urban outlines, night lights, road lights, water) is used when present and skipped otherwise. Aircraft, satellites, forecast precipitation, satellite clouds, tropical cyclones, meteor trails, and Moon/Sun hover images are not used. An omitted location uses the last location, then `Tokyo`; `--place` and search use only cached or local data. Also available in `zstarview-export-image` and as `Offline` in the startup dialog. | |
| `--clear-long-lived-cache` | Troubleshooting option. Delete long-lived DEM and urban-outline caches before startup. If used again within 3 days, startup is refused and the app tells you when retry is allowed. | |

#### Footnotes

\*3 Geo-satellite remains experimental because the upstream imagery is a display-oriented product rather than raw cloud data. Coastlines, borders, and other masked regions can have gaps that require inpainting or other fallback handling, so the path is useful but not yet as robust as the standard GOES/Himawari flow.
