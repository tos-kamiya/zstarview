# Overlays

| Option | Description | Default |
| :----- | :---------- | :------ |
| `-c`, `--cloud-opacity CLOUD_OPACITY` | Opacity of cloud rendering (0.0–1.0). Use 0.0 to disable for the session, even when `--geo-satellite true` is enabled. At night, the effective opacity is smoothly raised by up to 30% according to solar altitude to keep clouds visible. \*2 | `0.05` |
| `--cloud-stripe MODE[,COUNT[,WIDTH]]` | Cloud stripe style. `width` draws centered symmetric stripes whose visible width varies continuously with cloud amount; `width-quantized` is its five-level variant with gaps at width transitions and round line caps; `alpha` keeps width fixed and varies stripe alpha. `COUNT` is the absolute number of stripes drawn across the disc; it is no longer scaled with window or star-layer surface size. Stripe spacing therefore grows with the window, so the visual density decreases on larger windows while the stripe count stays fixed. `halftone` expands to `halftone,30,1.7`; `width` and `width-quantized` expand to `width,50,0.85` and `width-quantized,50,0.85`; `alpha` expands to `alpha,50,0.25`. If count or width is `0`, cloud rendering is disabled. | `halftone,30,1.7` |
| `--cloud-missing-tint-opacity OPACITY` | Opacity of missing-cloud-data yellow tint (0.0–1.0). | `0.176` |
| `--night-light-opacity OPACITY` | Opacity of the street-light portion of the NASA night lights overlay (0.0–1.0). Use 0.0 to disable the on-demand Black Marble download and street-light drawing for that run. | `0.14` |
| `--ridge-glow-opacity OPACITY` | Opacity of the ridge glow layer derived from the night-light profile (0.0–1.0). Use 0.0 to disable ridge glow rendering for that run. | `0.04` |
| `--water-surface-opacity OPACITY` | Opacity of the water-surface dots (0.0–1.0). Use 0.0 to disable the on-demand water-surface download and drawing for that run. \*5 | `0.4` |
| `-a`, `--aircraft-opacity OPACITY` | Opacity of the aircraft overlay (0.0–1.0). Use a positive value to explicitly enable aircraft queries and drawing; 0.0 disables them. | `0.0` |
| `--satellite-opacity OPACITY` | Opacity of the artificial satellite overlay (0.0–1.0). Use 0.0 to disable satellite element fetch and drawing for that run. | `0.5` |
| `--meteor-trails-opacity OPACITY` | Opacity of GMN meteor trails (0.0–1.0). A positive value enables loading and drawing; 0.0 disables fetching, drawing, and menu re-enabling for that run. | `0.0` |
| `--tropical-cyclone-opacity OPACITY` | Opacity of the tropical cyclone overlay (0.0–1.0). Use 0.0 to disable cyclone API fetch and drawing for that run. The overlay is also hidden automatically for time-shifted views. | `0.7` |
| `--show-guidelines-initial true\|false` | Whether guideline overlays are shown at startup. This controls the geometric horizon, celestial equator, ecliptic, never-rises circle, direction labels, zenith marker, and celestial pole markers. | `show` |
| `-d`, `--terrain-horizon-opacity OPACITY` | Opacity of the terrain horizon polyline (0.0–1.0). Use 0.0 to disable DEM download, terrain-horizon calculation, terrain-horizon drawing, and the Earth guide. \*4 | `0.003` |
| `-e`, `--earth-guide-opacity OPACITY` | Opacity of the below-horizon Earth guide line layer (0.0–1.0). Use 0.0 to disable Earth guide drawing for that run. \*4 | `0.028` |
| `--ground-tint-opacity OPACITY` | Strength of the ground-color fill below the geometric/terrain horizon (0.0–1.0). | `0.025` |
| `--overlay-font-size POINTS` | Base font size for window-drawn labels and HUD text only. This affects canvas text such as overlay info, planet labels, DSO labels, asterism names, search labels, and direction labels, but not status line text, menus, dialogs, or standard Qt widgets. Decimal values are accepted, and Qt may round the requested precision slightly depending on platform. | `11` |
| `-u`, `--urban-outline-opacity OPACITY` | Opacity of the urban outline overlay (0.0–1.0). Use 0.0 to disable it for that run. | `0.2` |
| `-r`, `--urban-outline-radius-km RADIUS_KM` | Fetch and render urban-outline buildings within this radius from the observer location. | `2.5` |
| `--urban-outline-skyscraper-radius-km RADIUS_KM` | Outer radius of the far-range skyscraper helper layer. Use `0` to disable skyscraper-tile lookup for that run; otherwise the value must be greater than or equal to `--urban-outline-radius-km`. | `60.0` |
| `--urban-outline-max-candidates N` | Keep at most `N` urban-outline ring candidates before expensive ring sampling. This is the main performance knob for outline pruning; use `0` to disable the layer. | `5000` |
| `--road-light-max-candidates N` | Keep at most `N` road `way` candidates before expensive simplification, terrain sampling, and projection. The limit does not truncate the road cache; use `0` to disable the Road Lights layer. Candidate counts are not shown in the normal status line. | `5000` |

<a id="about-water-surface"></a>

#### Footnotes

\*2 Cloud rendering uses infrared data from meteorological satellites (**Himawari** and **NOAA GOES** series), retrieved from their public S3 buckets. See Troubleshooting for tips on slow networks or offline use (for example, disabling clouds with `-c 0`). If Geo-satellite is enabled, `-c 0` still keeps cloud rendering disabled until the user re-enables it manually.

\*5 Water-surface rendering uses two data paths: sea dots come from local sea-mask tiles derived from [OSM Water Polygons](https://osmdata.openstreetmap.de/data/water-polygons.html), while inland water dots come from OpenStreetMap features fetched through the [Overpass API](https://overpass-api.de/). When inland water polygons do not carry an explicit `ele` or `water_level`, the renderer can sample DEM-backed ground height for the projected water surface. Some coastal, higher viewpoints can make a single Overpass fetch return around 9 MB of inland water data. The raw footprint is cached by viewpoint, but repeatedly changing viewpoints can accumulate download volume and may approach the public Overpass instance's rough 1 GB/day safety guideline.

\*4 Terrain horizon rendering downloads Copernicus DEM tiles on first use and reuses the cached DEM later. When enabled, the terrain profile also becomes the boundary for the ground-color fill inside the disc. The Earth guide is a separate layer with its own opacity toggle, and both layers currently share the dark olive-brown ground tone from the palette.
