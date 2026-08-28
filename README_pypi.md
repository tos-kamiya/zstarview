# zstarview

See the sky from anywhere, at any time -- even when it is cloudy or the Sun is up.

**Zenith Star View** is a desktop sky-scene simulator that combines the celestial sphere with atmospheric conditions, terrain, cityscapes, and other features around the horizon for a chosen location.

It renders an all-sky view with stars, the Sun, Moon, planets, deep-sky objects, and guide overlays.
Locations can be set by city or viewpoint name, direct coordinates, online place search, or supported Google Maps URLs.
When enabled, it can also add real-time cloud imagery, terrain horizon, urban outlines, night lights, nearby aircraft, and the ISS/JWST/Voyager 1/Voyager 2/Parker/Europa Clipper/Lucy/Psyche/JUICE/Solar Orbiter/BepiColombo artificial satellite overlays.
It also supports small blue-dot water surfaces and an optional false-color AKARI IR dust-map layer. The default AKARI display palette combines the 90 and 140 micrometre maps; the preparation command caches all four available bands by default.

https://github.com/user-attachments/assets/b0a4e340-1089-4256-9c48-b795d5c7b200

<p align="center">
  <video controls width="600" aria-label="Timelapse of the sky viewed from Matsue Station">
    <source src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/timelapse-matsueeki.mp4" type="video/mp4" />
  </video>
</p>

<p align="center"><em>Timelapse of the sky viewed from Matsue Station.</em></p>

## Screen descriptions

<p>This view shows the night sky over a Japanese city, displayed with <code>zstarview -p "Matsue Station" -A5 -Znw -a0.4 -P0.4</code>. The mouse is hovering near Vega, so the asterism it belongs to, the Summer Triangle, is shown. Buildings are displayed as an <em>urban outline</em>.</p>
<p>Clouds are visible on the right side of the image, rendered as a halftone pattern of circular dots. Clouds are rendered from cloud amounts estimated from satellite data such as GOES and Himawari.</p>
<p>Building data is generally obtained from Overture Maps, but <a href="https://github.com/tos-kamiya/zstarview#plateau-building-data-preparation">PLATEAU data</a> is also available for some Japanese cities. This screenshot uses PLATEAU data.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot1.png" width="60%" alt="Night sky over Matsue Station with the Summer Triangle asterism, clouds on the right, and urban outline"></p>

<p>This view shows Atlanta airport and the busy airspace above it, with more than ten aircraft visible. Aircraft trails are rendered as purple ribbons, and the ellipse below the horizon marks the <em>never-rises</em> region: the part of the celestial sphere that never comes above the horizon. The celestial equator and ecliptic are shown as dashed reference lines.</p>
<p>The celestial equator is a great circle whose normal axis is the Earth's rotation axis, the line connecting the north and south celestial poles. The never-rises boundary is a smaller circle around the same axis.</p>
<p>Note: To display aircraft as ribbons, specify <code>-a 0.4</code> or another opacity value in the CLI. In <code>zstarview-gui</code>, set <code>Aircraft opacity</code> to 0.4 or another value in the <code>Overlays</code> tab of the startup dialog.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot4.png" width="60%" alt="Atlanta airport and busy airspace with aircraft trails and the never-rises region"></p>

<p>This view looks straight up from Salar de Uyuni, known as one of the flattest places in the world. The horizon forms an almost complete circle around the view because the difference in elevation along the horizon is very small. The <code>-V9</code> option displays stars down to visual magnitude 9, while <code>-s5</code> makes the stars slightly larger than the default so that the smaller stars stand out.</p>
<p>Note: higher magnitude limits increase rendering time. See <a href="https://github.com/tos-kamiya/zstarview#about-magnitude-limit">about magnitude limit</a>.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot3.png" width="60%" alt="View of the sky and nearly circular horizon from Salar de Uyuni"></p>

<p>This view is from a tower 108 meters above the ground, with the viewing altitude lowered slightly to <code>-A-5</code> degrees. Some towers, such as Kobe Port Tower, have their locations and heights registered in the internal database. Use <code>--list-viewpoints t</code> to list the available tower names, or <code>--list-viewpoints m</code> to list the available mountain names. This view also uses <em>night lights</em>, a dataset of the brightness of the ground as seen from space, which makes the ground glow differently over urban and sea areas.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot9.png" width="60%" alt="Night-light ground glow viewed from a 108-meter tower at -5 degrees altitude"></p>

<p>This view uses <code>-A40</code> to change the altitude of the viewing direction and shows a view looking upward into the sky. The field of view reaches 90 degrees at the edge of the screen, producing a subtle fisheye-lens effect. The Moon is displayed in the upper-left part of the image.</p><p>The Moon is enlarged to 5x its normal apparent size after the mouse is moved over it. During normal operation, a simulated Moon image for the display time is fetched from NASA Scientific Visualization Studio's Dial-A-Moon API. The simulation includes the Moon's phase as the shadowed portion of the image.</p><p>In the default bright-body mode, the normal Moon marker uses a compact phase-aware outline. Unlike stars and planets, which are drawn according to visual magnitude, the Moon is rendered as a disk using its apparent angular size.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot11.png" width="60%" alt="Moon enlarged by mouse hover in the upper-left of an upward-looking view over Matsue"></p>

<p>This view shows a mountainous region of Switzerland, where steep terrain is visible along the horizon. Because the Sun is above the horizon, the sky is colored rather than fully dark. Terrain ridges are drawn in olive.</p><p>For cloud rendering, most of Europe is outside the coverage of GOES and Himawari, so clouds are not drawn by default. Experimental support is available through <code>--geo-satellite true</code>, which processes cloud imagery from a geostationary satellite.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot10.png" width="60%" alt="Mountain terrain, colored sky, and halftone cloud overlay in Switzerland"></p>

<p>This view is a star-field image generated with <code>zstarview-export-image</code>, rather than a GUI application screenshot. The object search option <code>--search "Torifune"</code> is used to show the position of the minor body. Osaka Castle, a Japanese building, is visible on the right.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot6.png" width="60%" alt="zstarview-export-image output showing Torifune and Osaka Castle"></p>

<p>This view shows a typhoon moving near Tokyo. The red marker shows the typhoon, and the blue diagonal lines show precipitation predicted by a forecast model. Unlike observation-based layers such as clouds, this forecast may differ from the precipitation that actually occurs. The number of lines represents the forecast amount of precipitation.</p>
<p>Forecast precipitation is enabled by setting <code>--precipitation-opacity</code> to a positive value. The first time it is enabled, you must agree to the Open-Meteo Free API terms.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot12.png" width="60%" alt="Typhoon and forecast precipitation viewed from Tokyo Tower"></p>

<p>This view shows observed Global Meteor Network meteor trails over Budapest. The trails are drawn in the directions recorded at observation time, with brighter tapered marks for the meteor paths and compact age labels such as <code>-56h</code> and <code>-67h</code>. The image was captured during the Perseid meteor shower, so an especially large number of trails is visible.</p>
<p>Meteor trails are enabled with a positive <code>--meteor-trails-opacity</code> value. Use <code>--meteor-trails-max-candidates</code> to limit the number of displayed trails.</p>
<p align="center"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot14.png" width="60%" alt="Global Meteor Network trails viewed from Budapest"></p>

## Installation (Recommended: `pipx`)

Zstarview is intended to be installed using [`pipx`](https://pypa.github.io/pipx/).

```bash
pipx install zstarview
```

Upgrade:

```bash
pipx upgrade zstarview
```

> Note: Linux x86_64 is the primary tested platform. The cloud-disc path no
> longer depends on `pyresample`, so its previous Windows on Arm64 installation
> blocker has been removed.

First run check:

```bash
zstarview-gui
```

This opens the startup dialog first. Use `zstarview-gui` for the default GUI
launch flow after installation. If you select `City` as the location source
and press `Auto Search`, the startup dialog fills in your current location
automatically.

- Optional coastline data: `zstarview-download-coastline --all` ([details](https://github.com/tos-kamiya/zstarview#1-optional-coastline-data))
- Optional AKARI IR bands data: `zstarview-download-akari-ir-bands` ([details](https://github.com/tos-kamiya/zstarview#2-optional-akari-ir-bands-data))
- Optional urban outline data: `pipxu install -f overturemaps`; Windows Arm64 users should see the [full instructions](https://github.com/tos-kamiya/zstarview#3-optional-urban-outline-data).

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
