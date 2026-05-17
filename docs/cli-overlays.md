# Overlays

| Option | Description | Default |
| :----- | :---------- | :------ |
| `-c`, `--cloud-opacity CLOUD_OPACITY` | Opacity of cloud rendering (0.0–1.0). Use 0.0 to disable. \*2 | `0.07` |
| `--cloud-stripe MODE[,COUNT[,WIDTH]]` | Cloud stripe style. `width` draws centered symmetric stripes whose visible width varies with cloud amount; `alpha` keeps width fixed and varies stripe alpha. `COUNT` is treated as the stripe density for the default 600x600 star render surface, and the effective count is scaled to match the star layer's downsampled surface size. `width` expands to `width,50,0.85`; `alpha` expands to `alpha,50,0.25`. If count or width is `0`, cloud rendering is disabled. | `width,50,0.85` |
| `--cloud-missing-tint-opacity OPACITY` | Opacity of missing-cloud-data yellow tint (0.0–1.0). | `0.176` |
| `--night-light-opacity OPACITY` | Opacity of the NASA night lights overlay (0.0–1.0). Use 0.0 to disable the on-demand Black Marble download and drawing for that run. | `0.022` |
| `--water-surface-opacity OPACITY` | Opacity of the water-surface dots (0.0–1.0). Use 0.0 to disable the on-demand water-surface download and drawing for that run. | `0.12` |
| `-a`, `--aircraft-opacity OPACITY` | Opacity of the aircraft overlay (0.0–1.0). Use 0.0 to disable aircraft queries and drawing for that run. | `0.5` |
| `--satellite-opacity OPACITY` | Opacity of the artificial satellite overlay (0.0–1.0). Use 0.0 to disable satellite element fetch and drawing for that run. | `0.5` |
| `--show-guidelines-initial true\|false` | Whether guideline overlays are shown at startup. This controls the geometric horizon, celestial equator, ecliptic, never-rises circle, direction labels, and zenith marker. | `show` |
| `-d`, `--terrain-horizon-opacity OPACITY` | Opacity of the terrain horizon polyline (0.0–1.0). Use 0.0 to disable DEM download, terrain-horizon calculation, terrain-horizon drawing, and the Earth guide. \*4 | `0.003` |
| `-e`, `--earth-guide-opacity OPACITY` | Opacity of the below-horizon Earth guide line layer (0.0–1.0). Use 0.0 to disable Earth guide drawing for that run. \*4 | `0.028` |
| `--ground-tint-opacity OPACITY` | Strength of the ground-color fill below the geometric/terrain horizon (0.0–1.0). | `0.1` |
| `-u`, `--urban-outline-opacity OPACITY` | Opacity of the urban outline overlay (0.0–1.0). Use 0.0 to disable it for that run. | `0.2` |
| `--urban-outline-feature-type {both,building}` | Overture cache mode for the urban outline. `both` combines `building` and `building_part`, preferring parts when available. | `both` |
| `-r`, `--urban-outline-radius-km RADIUS_KM` | Fetch and render urban-outline buildings within this radius from the observer location. The value is also part of the cache key. | `2.5` |
| `--urban-outline-skyscraper-radius-km RADIUS_KM` | Outer radius of the far-range skyscraper helper layer. Use `0` to disable skyscraper-tile lookup for that run; otherwise the value must be greater than or equal to `--urban-outline-radius-km`. | `60.0` |
| `-b`, `--urban-outline-min-building-height-m METERS` | Ignore buildings lower than this height when fetching/caching the urban outline. The value is also part of the cache key. | `0.0` |

#### About Water Surface

Water-surface rendering uses two different data paths: sea dots come from local sea-mask tiles derived from [OSM Water Polygons](https://osmdata.openstreetmap.de/data/water-polygons.html), while inland water dots come from OpenStreetMap features fetched through the [Overpass API](https://overpass-api.de/).

Some coastal, higher viewpoints can make a single Overpass fetch return around 9 MB of inland water data. The raw footprint is cached by viewpoint, but repeatedly changing viewpoints can accumulate download volume and may approach the public Overpass instance's rough 1 GB/day safety guideline.

#### Footnotes

\*2 Cloud rendering uses infrared data from meteorological satellites (**Himawari** and **NOAA GOES** series), retrieved from their public S3 buckets. See Troubleshooting for tips on slow networks or offline use (for example, disabling clouds with `-c 0`).

\*4 Terrain horizon rendering downloads Copernicus DEM tiles on first use and reuses the cached DEM later. When enabled, the terrain profile also becomes the boundary for the ground-color fill inside the disc. The Earth guide is a separate layer with its own opacity toggle, and both layers currently share the dark olive-brown ground tone from the palette.
