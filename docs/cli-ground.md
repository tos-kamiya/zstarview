# Ground

| Option | Description | Default |
| :----- | :---------- | :------ |
| `-d`, `--terrain-horizon-opacity OPACITY` | Opacity of the terrain horizon polyline (0.0–1.0). Use 0.0 to disable DEM download, terrain-horizon calculation, terrain-horizon drawing, and the Earth guide. \*4 | `0.003` |
| `-e`, `--earth-guide-opacity OPACITY` | Opacity of the below-horizon Earth guide line layer (0.0–1.0). Use 0.0 to disable Earth guide drawing for that run. \*4 | `0.028` |
| `--ground-tint-opacity OPACITY` | Strength of the ground-color fill below the geometric/terrain horizon (0.0–1.0). | `0.025` |
| `--water-surface-opacity OPACITY` | Opacity of the water-surface dots (0.0–1.0). Use 0.0 to disable the on-demand water-surface download and drawing for that run. \*5 | `0.4` |
| `--night-light-opacity OPACITY` | Opacity of the street-light portion of the NASA night lights overlay (0.0–1.0). Use 0.0 to disable the on-demand Black Marble download and street-light drawing for that run. | `0.14` |
| `--ridge-glow-opacity OPACITY` | Opacity of the ridge glow layer derived from the night-light profile (0.0–1.0). Use 0.0 to disable ridge glow rendering for that run. | `0.04` |
| `--road-light-opacity OPACITY` | Opacity of the OSM road lights layer (0.0–1.0). Use 0.0 to disable road light rendering for that run. | `0.12` |
| `--road-light-max-candidates N` | Keep at most `N` road `way` candidates before expensive processing. Use `0` to disable the layer. | `5000` |
| `-u`, `--urban-outline-opacity OPACITY` | Opacity of the urban outline overlay (0.0–1.0). Use 0.0 to disable it for that run. | `0.2` |
| `-r`, `--urban-outline-radius-km RADIUS_KM` | Fetch and render urban-outline buildings within this radius from the observer location. | `2.5` |
| `--urban-outline-skyscraper-radius-km RADIUS_KM` | Outer radius of the far-range skyscraper helper layer. Use `0` to disable skyscraper-tile lookup; otherwise the value must be at least `--urban-outline-radius-km`. | `60.0` |
| `--urban-outline-max-candidates N` | Keep at most `N` urban-outline ring candidates before expensive ring sampling. Use `0` to disable the layer. | `5000` |
| `--urban-outline-skyscraper-only` | Draw only the far-range skyscraper helper layer and skip the normal near-range urban outline. | off |
| `--urban-outline-download-timeout-seconds SECONDS` | Maximum time allowed for each Overture urban-outline download subprocess. | `120.0` |

<a id="about-water-surface"></a>

#### Footnotes

\*4 Terrain horizon rendering downloads Copernicus DEM tiles on first use and reuses the cached DEM later. When enabled, the terrain profile also becomes the boundary for the ground-color fill inside the disc. The Earth guide is a separate layer.

\*5 Water-surface rendering uses local sea-mask tiles for sea dots and OpenStreetMap features fetched through the Overpass API for inland water dots. The raw footprint is cached by viewpoint.
