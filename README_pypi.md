# zstarview

**Zenith Star View** is a desktop sky viewer for your chosen location.

It renders an all-sky view with stars, the Sun, Moon, planets, deep-sky objects, and guide overlays.
When enabled, it can also add real-time cloud imagery, terrain horizon, urban outlines, night lights, nearby aircraft, and the ISS/JWST/Voyager 1/Voyager 2/Parker/Europa Clipper/Lucy/Psyche/JUICE/Solar Orbiter/BepiColombo artificial satellite overlays.
Locations can be set by city or viewpoint name, direct coordinates, online place search, or supported Google Maps URLs.

## Screen descriptions

<p>This view shows the night sky over a Japanese city, displayed with <code>-p "Matsue Station" -A5 -Anw</code>. The mouse is hovering near Dubhe, so the asterism it belongs to, the Big Dipper, is shown. Buildings are displayed as an <em>urban outline</em>.</p>
<p>Building data is generally obtained from Overture Maps, but <a href="https://github.com/tos-kamiya/zstarview#plateau-building-data-preparation">PLATEAU data</a> is also available for some Japanese cities. This screenshot uses PLATEAU data.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot1.png" width="60%" alt="Night sky over Matsue Station with the Big Dipper asterism and urban outline"></p>

<p>This view shows Atlanta airport and the busy airspace above it, with more than ten aircraft visible. Aircraft trails are rendered as purple ribbons, and the ellipse below the horizon marks the <em>never-rises</em> region: the part of the celestial sphere that never comes above the horizon. The celestial equator and ecliptic are shown as dashed reference lines.</p>
<p>The celestial equator is a great circle whose normal axis is the Earth's rotation axis, the line connecting the north and south celestial poles. The never-rises boundary is a smaller circle around the same axis.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot4.png" width="60%" alt="Atlanta airport and busy airspace with aircraft trails and the never-rises region"></p>

<p>This view looks straight up from Salar de Uyuni, known as one of the flattest places in the world. The horizon forms an almost complete circle around the view because the difference in elevation along the horizon is very small. The <code>-V10.5</code> option displays stars down to visual magnitude 10.5, while <code>-s5</code> makes the stars slightly larger than the default so that the smaller stars stand out.</p>
<p>Note: higher magnitude limits increase rendering time. See <a href="https://github.com/tos-kamiya/zstarview#about-magnitude-limit">about magnitude limit</a>.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot3.png" width="60%" alt="View of the sky and nearly circular horizon from Salar de Uyuni"></p>

<p>This view is from a tower 108 meters above the ground, with the viewing altitude lowered slightly to <code>-A-5</code> degrees. Some towers, such as Kobe Port Tower, have their locations and heights registered in the internal database. This view also uses <em>night lights</em>, a dataset of the brightness of the ground as seen from space, which makes the ground glow differently over urban and sea areas.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot9.png" width="60%" alt="Night-light ground glow viewed from a 108-meter tower at -5 degrees altitude"></p>

<p>This view shows the Moon enlarged to 5x its normal apparent size after the mouse is moved over it. The enlargement makes the Moon's lunar phase—the shape of its illuminated and dark portions—easier to see. Unlike stars and planets, which are drawn according to visual magnitude, the Moon is rendered as a disk using its apparent angular size. Cloud cover is visible over roughly the eastern third of the sky, rendered as a halftone pattern of circular dots.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot11.png" width="60%" alt="Moon enlarged by mouse hover with partial cloud cover over Matsue"></p>

<p>This view shows a mountainous region of Switzerland, where steep terrain is visible along the horizon. Because the Sun is above the horizon, the sky is colored rather than fully dark. Terrain ridges are drawn in olive. Clouds are rendered from cloud amounts estimated from satellite data such as GOES and Himawari. Most of Europe is outside their usual coverage, but experimental support is available through <code>--geo-satellite true</code>, which processes cloud imagery from a geostationary satellite.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot10.png" width="60%" alt="Mountain terrain, colored sky, and halftone cloud overlay in Switzerland"></p>

<p>This view is a star-field image generated with <code>zstarview-export-image</code>, rather than a GUI application screenshot. The object search option <code>--search "Torifune"</code> is used to show the position of the minor body. Osaka Castle, a Japanese building, is visible on the right.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot6.png" width="60%" alt="zstarview-export-image output showing Torifune and Osaka Castle"></p>

<p>This view uses <code>-A45</code> to change the altitude of the viewing direction and shows a view looking slightly up into the sky. The field of view reaches 90 degrees at the edge of the screen, producing a subtle fisheye-lens effect.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot8.png" width="60%" alt="A slightly upward-looking sky view rendered with -A45"></p>

## Install

Recommended:

```bash
pipx install zstarview
```

Or with pip:

```bash
pip install zstarview
```

> Note: The previous Windows on Arm64 installation blocker has been removed.
> Installation is now possible there, but Windows Security may still block
> Python extension modules during startup on some systems. If that happens,
> see the troubleshooting notes in the main project README.

First run check:

```bash
zstarview-gui
```

This opens the startup dialog first. If you select `City` as the location
source and press `Auto Search`, the dialog fills in your current location
automatically.

Prerequisite for the urban outline overlay:

For non-Arm64 platforms, install the `overturemaps-py` package with `pipx`:

```bash
pipx install overturemaps-py
```

On Windows Arm64, stage an x64 `overturemaps` v1.0.1 or newer release executable into the zstarview cache after installing `zstarview` with `pipx`; see [Windows on Arm64](https://github.com/tos-kamiya/zstarview#windows-on-arm64) for the full steps and the dependency-wheel note.

## Highlights

- Stars: the sky view shows stars from the selected catalog, with asterisms and other overlays layered on top.
- Solar-system bodies: supports Sun, Moon, and major planets.
- Deep-sky objects: named galaxies/open clusters/globular clusters are shown as soft blue extents.
- Asterism overlay: popular line patterns are shown as dim guide lines, alongside other sky guides such as the never-rises region, which uses a distinct neutral gray, plus the zenith and celestial pole markers.
- Satellite cloud imagery and sky-color disc: real-time Himawari/GOES satellite data are downloaded and rendered as a stylized hatched overlay, with the sky-color disc still visible beneath the clouds. The default `width` preset uses 50 stripes whose visible width varies with cloud amount. An experimental Geo-satellite path is available for Europe-band observers when `--geo-satellite true` is set.
- Aircraft and artificial satellite overlays: nearby aircraft from OpenSky can be drawn on the sky view, and ISS, JWST, Voyager 1, Voyager 2, Parker, Europa Clipper, Lucy, Psyche, JUICE, Solar Orbiter, and BepiColombo can be drawn as small purple markers between the planet and aircraft layers. ISS uses `wheretheiss.at` with CelesTrak fallback, while the other spacecraft use JPL Horizons.
- Tropical cyclone overlay: active hurricanes / typhoons from a public ArcGIS `Active_Hurricanes_v1` FeatureServer can be shown as small markers with projected current-position tracking and a distance cutoff.
- Fresh current caches are reused for up to 24 hours for both the ISS cache and the Horizons-backed spacecraft cache.
- Terrain horizon and earth guide: Copernicus DEM data can be used to render the local terrain skyline and ground region, with a separate below-horizon continental hatch layer for orientation.
- Urban outline overlay: major rooflines are drawn for the current viewpoint, with optional distant skyscrapers in dense urban areas.
- Night lights overlay: The 2025 annual VIIRS Nighttime Lights VNL v2.2 GeoTIFF is downloaded from a GitHub Release on demand, cached locally, and rendered as a separate glow layer above the horizon and terrain ridges.
- Flexible location input: start from a city, tower, mountain, lat/lon, or online place/station search.
- Adjustable view center: change the view center from the CLI or with the arrow keys.
- Python support: routinely tested on CPython 3.10, 3.11, 3.12, 3.13, and 3.14.

## Entry Points

This package provides three entry points:

- `zstarview-gui` for the interactive startup dialog and GUI launcher, with its own startup profile separate from the legacy CLI config.
- `zstarview` for GUI launches driven by command-line options, sharing the same option families described in the main README and the CLI reference docs.
- `zstarview-export-image` for headless one-shot image export, using the same option families as `zstarview` where applicable.

## Code, Data Licenses, and Credits

- Code: MIT License. See `LICENSE`.
- Bundled and runtime-fetched data may be subject to their own licenses, attribution rules, or service terms.
- See the main project README for the full credits and third-party data notes.

## Links

- Homepage: https://github.com/tos-kamiya/zstarview
- Documentation: https://github.com/tos-kamiya/zstarview#readme
- Issues: https://github.com/tos-kamiya/zstarview/issues
