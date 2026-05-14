# General CLI Options

| Option | Description | Default |
| :----- | :---------- | :------ |
| `-h`, `--help` | Show this help message and exit. | |
| `--window-geometry restore\|X,Y,W,H` | Set initial window geometry. Use `restore` to load the last saved position/size, or `X,Y,W,H` to specify explicit integers. Note: on Wayland, window position restore is not available (size restore works). | |
| `--window-frame {frameless,window}` | Choose window decorations. `frameless` keeps the current borderless presentation; `window` uses the platform title bar and frame. | `frameless` |
| `--observation-info auto\|top\|bottom\|off` | Startup mode for the observation-info block. | `bottom` |
| `--include-direction-grid` | `zstarview-export-image` only. Include the direction grid in exported images, with 30-degree major lines and 10-degree intersection crosses. | |
| `-t`, `--theme {night,day,white,black,transparent}` | Theme preset for background and star contrast. | `night` |
| `--visibility-boost MULTIPLIER` | Visibility boost for faint support layers. Values above `1.0` raise opacity for layers such as the terrain horizon, Earth guide, urban outline, sky disc, cloud disc, satellites, aircraft, and ground tint. | `1.0` |
| `--clear-long-lived-cache` | Troubleshooting option. Delete long-lived DEM and urban-outline caches before startup. If used again within 3 days, startup is refused and the app tells you when retry is allowed. | |

\*1 When using non-realtime sky options (`--hours`, `--days`, `--datetime`), cloud, aircraft, and artificial satellite overlays are not shown.

\*2 Cloud rendering uses infrared data from meteorological satellites (**Himawari** and **NOAA GOES** series), retrieved from their public S3 buckets.
   See Troubleshooting for tips on slow networks or offline use (e.g., disabling clouds with `-c 0`).

\*3 The brightest-magnitude multiplier cannot exceed the classical Pogson value of \(100^{1/5}\approx2.512\).

\*4 Terrain horizon rendering downloads Copernicus DEM tiles on first use and reuses the cached DEM later. When enabled, the terrain profile also becomes the boundary for the ground-color fill inside the disc. The Earth guide is a separate layer with its own opacity toggle, and both layers currently share the dark olive-brown ground tone from the palette.

\*5 `--place` uses the public OpenStreetMap Nominatim search service. It sends a single search request with a User-Agent and Accept-Language. See the Nominatim usage policy if you plan to rely on this option heavily or from automation.
