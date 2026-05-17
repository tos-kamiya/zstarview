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
