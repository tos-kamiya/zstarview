# zstarview 🌌

See the starry sky, even when it's cloudy or the sun is out.

**Zenith Star View** is a desktop sky viewer for your chosen location.
The name emphasizes the *zenith*, the point directly overhead.

<table align="right">
  <tr>
    <td align="center">
      <img src="docs/images/clickpy-top25-20260523.png" alt="clickpy Top 25 medal" width="160" />
      <br />
      <sub>clickpy Top 25 as of 2026-05-23</sub>
      <br />
      <sub><a href="https://clickpy.clickhouse.com/dashboard/zstarview">View on clickpy</a></sub>
    </td>
  </tr>
</table>

It renders an all-sky view with stars, the Sun, Moon, planets, deep-sky objects, and guide overlays.
When enabled, it can also add real-time cloud imagery, terrain horizon, urban outlines, night lights, and nearby aircraft.
Locations can be set by city or viewpoint name, direct coordinates, online place search, or supported Google Maps URLs.

**Features:**

- **Stars**: the sky view shows stars from the selected catalog, with asterisms and other overlays layered on top.
- **Solar-system bodies**: supports Sun, Moon, and major planets. Minor planets (asteroids) are not displayed yet.
- **Search**: search named stars, asterisms, places, known artificial satellites, and JPL bodies from one dialog. When local matches are not found, known artificial satellites use the app's current position first; if a target is recognized as a known artificial satellite but its current position cannot be obtained, the search fails instead of falling back to JPL. Persistent marker mode keeps both the marker and label visible on the selected target.
- **Overlay groups**: the features below are grouped by display layer, and each overlay can be toggled individually from the menu or with its keyboard shortcut.
- **Flexible location and view center**: specify the observer location with a city name, tower name, mountain name, direct latitude/longitude input, supported Google Maps coordinate URLs, or online place/station search via Nominatim. Adjust the view center with `-A` / `-Z` or the arrow keys. The HUD also shows the observer location and the current `alt/az` view center.
- **Terminal image export**: `zstarview-export-image` can render the sky headlessly and write it to a file or display it directly in sixel-capable terminals.
- **Python support**: the project is routinely tested on CPython 3.10, 3.11, 3.12, 3.13, and 3.14.

**Celestial overlays:**

- **Deep-sky objects**: named galaxies/open clusters/globular clusters are shown as soft blue extents.
- **Asterism overlay**: popular line patterns rather than formal IAU constellation boundaries are shown as dim ambient lines. Mouse-hovering a star in an asterism brightens the matching pattern and shows its label, with 3-second rotation when multiple asterisms share that star.
- **Sky Guides**: guide overlays include the never-rises region as a guide-line style solid circle, and the celestial equator as a dashed line with longer on-segments in the same neutral gray, along with direction labels around the horizon, a zenith marker, and celestial pole markers.

**Atmospheric and man-made overlays:**

- **Sky Color**: the sky-color disc can be shown as a gradient sky disc or as a flat dark-disc fallback.
- **Clouds**: real-time Himawari/GOES satellite data are downloaded and rendered as a stylized hatched (striped) overlay, and the sky-color disc remains visible beneath the clouds. Missing regions are shown in faint yellow when satellite coverage is partial. See [an example with partial coverage and yellow missing-data tint](docs/images/screenshot5.png).
- **Artificial satellites**: ISS, JWST, Voyager 1, Voyager 2, and Parker can be drawn as small purple markers between the planet and aircraft layers.
- **Aircraft**: nearby aircraft from OpenSky can be drawn as purple predicted-motion polylines.

**Building and ground guide overlays:**

- **Earth guide**: an independent layer draws a simplified continental hatch pattern below the horizon in the same ground tone to help with orientation.
- **Terrain horizon**: Copernicus DEM data can be downloaded to render the local terrain skyline. The terrain overlay shows banded ridge lines in the same horizon color, with nearby bands drawn thicker and distant bands drawn thinner. Blue-tinted ridge lines mark the parts visible from the observer, not hidden by nearer ridges. The disc is filled with the same ground tone below the terrain horizon, or below the geometric horizon when terrain is disabled.
- **Urban outline**: major rooflines are drawn as a white urban outline overlay for the current viewpoint. In some skyscraper-heavy cities, distant skyscrapers can also be added from within a 60km radius.
- **Water surface**: nearby water bodies are rendered as small blue dots. Sea points come from OSM Water Polygons sea-mask tiles, while inland water points come from OpenStreetMap features fetched via Overpass API. See [About Water Surface](docs/cli-overlays.md#about-water-surface).
- **Night lights**: NASA Earth at Night / Black Marble 2016 Grayscale 500m GeoTIFF tiles are downloaded on demand, cached locally, and rendered as a separate glow layer above the horizon and terrain ridges.

## Screenshots

The first screenshot shows the asterism overlay and also serves as a terrain-horizon example.
The second screenshot shows the aircraft overlay together with the never-rises region.
The third screenshot shows a denser star field rendered with `-V10.5 -s4.5`.
The fourth screenshot shows `zstarview-export-image --search "C/1861 G1" --sixel` displaying the Thatcher comet (`C/1861 G1 Thatcher`) in a sixel terminal.

  <p align="center">
    <img src="docs/images/screenshot1.png" alt="Screenshot showing the asterism overlay and a terrain horizon example" width="49%" />
    <img src="docs/images/screenshot4.png" alt="Screenshot showing the aircraft overlay and the never-rises region" width="49%" />
  </p>

  <p align="center">
    <img src="docs/images/screenshot3.png" alt="Screenshot showing a denser star field rendered with -V10.5 -s4.5" width="49%" />
    <img src="docs/images/screenshot6.png" alt="Screenshot showing zstarview-export-image --search C/1861 G1 --sixel displaying the Thatcher comet (C/1861 G1 Thatcher) in a sixel terminal" width="49%" />
  </p>

Note: higher magnitude limits increase rendering time. See [About magnitude limit](docs/cli-sky-and-stars.md#about-magnitude-limit).

Urban outline and terrain horizon examples from several cities worldwide:

<table>
  <tr>
    <td align="center" width="25%"><img src="docs/images/screenshot-sydney.png" alt="Circular Quay, Sydney" width="100%" /></td>
    <td align="center" width="25%"><img src="docs/images/screenshot-tokyotower.png" alt="Near Tokyo Tower, Tokyo" width="100%" /></td>
    <td align="center" width="25%"><img src="docs/images/screenshot-mtfuji.png" alt="View of Mt. Fuji, Japan" width="100%" /></td>
    <td align="center" width="25%"><img src="docs/images/screenshot-marinabay.png" alt="Marina Bay, Singapore" width="100%" /></td>
  </tr>
  <tr>
    <td align="center"><sub>Circular Quay, Sydney</sub></td>
    <td align="center"><sub>Near Tokyo Tower, Tokyo</sub></td>
    <td align="center"><sub>View of Mt. Fuji, Japan</sub></td>
    <td align="center"><sub>Marina Bay, Singapore</sub></td>
  </tr>
</table>

<div style="height: 0.8em;"></div>

<table>
  <tr>
    <td align="center" width="25%"><img src="docs/images/screenshot-burjkhalifa.png" alt="Near Burj Khalifa, Dubai" width="100%" /></td>
    <td align="center" width="25%"><img src="docs/images/screenshot-jungfraujoch.png" alt="Jungfraujoch, Switzerland" width="100%" /></td>
    <td align="center" width="25%"><img src="docs/images/screenshot-manarola.png" alt="Manarola, Italy" width="100%" /></td>
    <td align="center" width="25%"><img src="docs/images/screenshot-sagradafamilia.png" alt="Sagrada Familia, Barcelona" width="100%" /></td>
  </tr>
  <tr>
    <td align="center"><sub>Near Burj Khalifa, Dubai</sub></td>
    <td align="center"><sub>Jungfraujoch, Switzerland</sub></td>
    <td align="center"><sub>Manarola, Italy</sub></td>
    <td align="center"><sub>Sagrada Familia, Barcelona</sub></td>
  </tr>
</table>

## Installation (Recommended: `pipx`)

Prerequisite for the urban outline overlay: install the `overturemaps` CLI separately.
Installation: <https://pypi.org/project/overturemaps/>

Confirm it with:

```bash
overturemaps --help
```

Zstarview is intended to be installed using [`pipx`](https://pypa.github.io/pipx/).

> Note: Linux x86_64 is the primary tested platform. The cloud-disc path no
> longer depends on `pyresample`, so its previous Windows on Arm64 installation
> blocker has been removed.

```bash
pipx install zstarview
```

Upgrade:

```bash
pipx upgrade zstarview
```

First run check:

```bash
zstarview-gui
```

This opens the startup dialog first. Use `zstarview-gui` for the default GUI
launch flow after installation.
If you select `City` as the location source and press `Auto Search`, the
startup dialog fills in your current location automatically.

> Note: Troubleshooting tips (including X11 libraries and slow network) are summarized below.

## Usage

Use `zstarview-gui` for the interactive startup dialog, `zstarview` for GUI
launches driven by CLI options, and `zstarview-export-image` for headless
one-shot image export.

After launching `zstarview-gui`, the `Auto` button in the Location tab uses
network connection information to estimate the current observer location and
fills it into the dialog.

The CLI reference below documents `zstarview` and `zstarview-export-image`.

```bash
zstarview [options] [location]
```

> Note (Ubuntu/Wayland, GNOME): If the taskbar icon does not appear when launching from a terminal, follow the steps in [Tools](#tools) under [Generating a `.desktop` launcher (GNOME only)](#generating-a-desktop-launcher-gnome-only).

Quick examples:

```bash
zstarview Tokyo
zstarview auto
zstarview "Tokyo Skytree"
zstarview "35.68;139.76" --datetime "2025-09-12 21 JST"
zstarview --place "Matsue Station" --place-countrycode jp
zstarview --search Ceres
```

### CLI Reference

#### Argument

| Format | Description | Default |
| :----- | :---------- | :------ |
| City Name | City name (e.g., `Tokyo`) | Last run location (or `Tokyo`) |
| Tower/Mountain | Tower (e.g., `t/Tokyo Skytree`) or Mountain (e.g., `m/Mount Fuji`) | |
| Coordinates | Direct coordinates (e.g., `35.68;139.76`, `@35.68,139.76`) or Google Maps URL | |
| `auto` | Automatically detect current location via IP | |

The links below cover the detailed option groups, and the linked docs files contain the detailed option tables, footnotes, and examples.

- [Observing Location and Time](docs/cli-observing-location-and-time.md)
- [Viewpoint Dataset Queries for Observing Locations](docs/cli-viewpoint-dataset-queries.md)
- [Search Objects at startup](docs/cli-search-objects.md)
- [Sky and Stars](docs/cli-sky-and-stars.md)
- [Overlays](docs/cli-overlays.md)
- [General](docs/cli-general.md)

<details>
  <summary>Tools</summary>

## Tools

### `zstarview-export-image`

Use `zstarview-export-image` for one-shot image export without starting the GUI.

```bash
zstarview-export-image Matsue -o matsue.png
```

`zstarview-export-image` writes the usual location/time/view/vmag summary to `stderr` after rendering, or immediately before terminal image output when `--sixel` is used.

For troubleshooting or manual cache maintenance, you can print the cache root directory without rendering:

```bash
zstarview-export-image --print-cache-dir
```

If you really need to bypass the `--clear-long-lived-cache` cooldown, first run `zstarview-export-image --print-cache-dir`, then remove these subdirectories under that cache root before startup:

- `copernicus-dem`
- `overture_buildings`
- `overture_skyscrapers`

### Generating a `.desktop` launcher (GNOME only)

On GNOME-based environments (including Ubuntu Dock and DockToPanel),
a `.desktop` file is required for the correct icon to appear in the taskbar.

This application includes a helper command to generate it:

```bash
# Create zstarview-gui.desktop in the current directory
zstarview-make-desktop-file

# Install to ~/.local/share/applications
zstarview-make-desktop-file --write
```

* Without `--write`, the file `zstarview-gui.desktop` is created in the current directory.
* With `--write`, it is installed to `~/.local/share/applications` and registered with the desktop database.

> **Note:** This launcher integration is only intended for GNOME-based environments.
> It is not required on other desktop environments, and may not work as intended elsewhere.

</details>

The GUI supports direct keyboard and menu-based navigation, search, and overlay toggles.

<details>
  <summary>GUI operations</summary>

### GUI

#### Key Operations

* **← / →**: Rotate view azimuth by ±5°
* **↑ / ↓**: Change view altitude by ±5° (clamped to -45°..90°)
* **Shift + arrow keys**: Fine-tune view direction by 1°
  While arrow-key input continues, the app keeps a simplified viewport-interaction mode for about 0.7 seconds after the last input. In this mode, it shows stars up to `Vmag <= 4.0`, the Sun, Moon, planets, the celestial equator, ecliptic, horizon, terrain horizon, direction labels, the zenith marker, and the celestial pole markers; full star density, sky-color disc, clouds, night lights, DSO, asterisms, and urban outlines are temporarily hidden.
* **M**: Toggle moon enlarged to 5x size
* **D**: Toggle DSO overlays
* **A**: Toggle asterism overlays
* **G**: Toggle sky guides
* **S**: Toggle sky color visibility
* **C**: Toggle cloud overlays
* **L**: Toggle night lights overlay
* **P**: Toggle aircraft overlay
* **I**: Toggle artificial satellite overlay
* **T**: Toggle terrain horizon overlay
* **E**: Toggle earth guide overlay
* **U**: Toggle urban outline overlay
* **Ctrl+J**: Open Jump to Named Star
* **Ctrl+F**: Open Search Objects
* **F11**: Toggle fullscreen display
* **ESC**: Exit fullscreen
* **Q**: Quit

#### Menu Operations

From the hamburger menu (`☰`), you can use:

* **Jump to Named Star...**: Choose from representative named stars (`Vmag <= 2.0`), grouped into North / Equatorial / South, then jump the view center to that star.
* **Search Objects...**: Search across named stars, supported asterisms, places, known artificial satellites, and JPL bodies, then jump to the selected target. If the local star/asterism search finds nothing, known artificial satellites use the app-side current position first; if a target is recognized as a known artificial satellite but its current position cannot be obtained, the search fails instead of falling back to JPL. Use `Keep marker` to keep both the marker and label visible after the jump.
* **Search Places...**: Open a separate place-search dialog backed by OpenStreetMap Nominatim, list matching places/stations/facilities, and jump toward the selected ground location.
* **Set View Center...**: Open a direct Alt/Az dialog with the current values prefilled, then apply the entered view center immediately.
* **Enlarge Moon**: Toggle moon enlarged to 5x size.
* **DSO**: Toggle deep-sky object overlays on/off.
* **Asterisms**: Toggle asterism overlays on/off (when enabled, dim overlays stay visible; hovering a member star brightens the matching asterism and shows its label).
* **Guidelines**: Toggle the geometric horizon, celestial equator, ecliptic, never-rises solid circle, direction labels, zenith marker, and celestial pole markers on/off. The celestial equator uses a longer dashed stroke in the same neutral gray as the never-rises circle.
* **Observation Info**: Toggle the observation-info block on/off. When shown, it stays at the bottom-left by default; `auto` keeps the older hover-avoid placement behavior.
* **Sky Color**: Switch between the full sky-color gradient and the flat dark-disc fallback.
* **Clouds**: Toggle real-time cloud overlays on/off.
* **Night Lights**: Toggle the NASA Earth at Night / Black Marble overlay on/off. If disabled from the CLI with `--night-light-opacity 0`, the menu item cannot re-enable it for that run.
* **Aircraft**: Toggle the OpenSky-based aircraft overlay on/off. If disabled from the CLI with `-a 0` / `--aircraft-opacity 0`, the menu item cannot re-enable it for that run.
* **Satellites**: Toggle the ISS artificial satellite overlay on/off. If disabled from the CLI with `--satellite-opacity 0`, the menu item cannot re-enable it for that run.
* **Terrain Horizon**: Toggle the terrain skyline overlay on/off. If disabled from the CLI with `-d 0` / `--terrain-horizon-opacity 0`, the menu item cannot re-enable it for that run.
* **Earth Guide**: Toggle the below-horizon earth-guide overlay on/off. If disabled from the CLI with `-e 0` / `--earth-guide-opacity 0`, the menu item cannot re-enable it for that run.
* **Urban Outline**: Toggle the Overture-derived urban roofline overlay on/off. If disabled from the CLI with `-u 0` / `--urban-outline-opacity 0`, the menu item cannot re-enable it for that run.
* **Fullscreen**: Toggle fullscreen display.
* **Exit**: Quit the application.

After a jump/search, the selected star is highlighted for about 3 seconds using the same UI style as mouse hover (circle marker + name label).

The same simplified viewport-interaction mode is also used during window resize to keep interaction responsive while heavier layers catch up after idle.

</details>

<details>
  <summary>Urban Outline Data</summary>

`zstarview` now fetches urban-outline source data on demand from Overture Maps and
caches the derived building tiles under the app cache directory. The first launch
for a new viewpoint/radius/height combination may take a few seconds while the
download finishes; the outline appears automatically after the cache is ready.

This path requires the `overturemaps` CLI to be installed separately. Confirm it
with:

```bash
overturemaps --help
```

Useful startup options:

```bash
zstarview "Tokyo Tower" -r 2.5 -b 0
zstarview -p "Matsue Station" -r 2.0 -b 20
```

- `-r`, `--urban-outline-radius-km`: fetch radius in kilometers
- `-b`, `--urban-outline-min-building-height-m`: minimum building height in meters
- `--urban-outline-feature-type`: choose which urban-outline data is used for display; default `both`

</details>

<details>
  <summary>Troubleshooting and platform notes</summary>

## Troubleshooting

### X11 (Ubuntu/Debian)
Qt's xcb platform plugin may require `libxcb-cursor0` at runtime.
If you're not watching for X11 vs Wayland differences, this can be confusing — running from a terminal may show errors like:

```sh
$ zstarview
qt.qpa.plugin: From 6.5.0, xcb-cursor0 or libxcb-cursor0 is needed to load the Qt xcb platform plugin.
qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.

Available platform plugins are: eglfs, offscreen, wayland-egl, linuxfb, wayland, minimal, xcb, vkkhrdisplay, minimalegl, vnc.
```

Install the missing `libxcb-cursor0` package with:

`sudo apt install libxcb-cursor0`

### Wayland Window Shadows
On some Wayland desktops, a normal framed `zstarview --window-frame window` window may appear without the usual outer shadow.
This is usually caused by the Wayland decoration/compositor path rather than by zstarview's own window settings.

If you prefer the X11-style shadowed window appearance, a practical workaround is to launch the app through XWayland:

```sh
QT_QPA_PLATFORM=xcb zstarview --window-frame window
```

If the shadow appears with `QT_QPA_PLATFORM=xcb`, that confirms the difference is in the Wayland vs X11 window-decoration path on your desktop environment.

### Wayland Flicker With `QT_QPA_PLATFORM=xcb`
On some Wayland systems, forcing `QT_QPA_PLATFORM=xcb` can make the frameless window flicker, especially when maximized.
The visible symptom is that the desktop or windows behind zstarview briefly show through during repaint.

If this happens, do not force `QT_QPA_PLATFORM=xcb` for the normal frameless UI.
Launch zstarview without that variable, or use the standard framed window mode if you specifically need XWayland:

```bash
zstarview
zstarview --window-frame window
```

### Slow or Unstable Network / Offline Use
1. Planetary ephemeris data

   On the very first launch, the app downloads a planetary ephemeris file (`de442s.bsp`).
   This requires network connectivity once. After it is cached, the app can run offline.

2. Cloud satellite imagery

   Cloud rendering downloads satellite imagery from public S3 buckets (Himawari / NOAA GOES) and relies on heavy dependencies.
   If your network is slow or unavailable, disable clouds with `-c 0`.
   You can still explore stars/planets and sky colors without cloud overlays.

3. Terrain horizon

   Terrain horizon rendering downloads Copernicus DEM tiles once and then reuses the local cache.
   If your network is slow or unavailable, disable terrain horizon rendering with `-d 0` or `--terrain-horizon-opacity 0`.
   You can still explore stars/planets and sky colors without terrain overlays.

4. Night lights data

   Night lights use NASA Earth at Night / Black Marble 2016 Grayscale 500m GeoTIFF tiles.
   The app downloads the tiles on demand, caches them locally, and reuses the cache on later launches.
   If your network is slow or unavailable, disable the layer with `--night-light-opacity 0`.
   If the cache is already present, the app can keep showing the night lights overlay without network access.

5. Artificial satellite data

   The artificial satellite overlay fetches ISS orbital data at runtime, using `wheretheiss.at` as the primary source and CelesTrak as a fallback, and fetches JWST, Voyager 1, Voyager 2, and Parker from JPL Horizons. Fresh current caches are reused for up to 24 hours for both the ISS cache and the Horizons-backed spacecraft cache.
   The layer is available only for realtime views; time-shifted views do not fetch or display artificial satellites.
   If your network is slow or unavailable, disable the layer with `--satellite-opacity 0`.
   If a fresh cache is already present, the app can keep showing the satellite overlay without network access.

6. Aircraft data

   The aircraft overlay fetches OpenSky Network state data at runtime.
   By default it refreshes once every 5 minutes. This interval is intentionally conservative so the app keeps practical headroom for free-tier use, temporary failures, and retries rather than polling more aggressively.
   If you want to avoid OpenSky queries entirely, disable the layer with `-a 0`.

Cloud-related status text uses `idle` / `downloading` / `partial`:
- `downloading`: fetching source imagery from S3
- `partial`: rendered with available data only; missing regions are tinted faint yellow

After the GOES-East refresh to GOES-19, some places that previously showed only a generic "satellite not covering this region" result can now render as partial coverage instead. This now includes some locations in Europe. In those cases, covered parts of the sky show cloud imagery and uncovered parts show the faint yellow missing-data tint. See [screenshot5](docs/images/screenshot5.png) for an example around 77-87% coverage.

### Sky Update Interval and CPU Load
Frequent sky updates can be CPU‑intensive on lower‑end machines. Increase the interval to reduce load (e.g., `-i 300` for every 5 minutes). Lower it only if your machine can keep up.

### Viewing Logs
Launching from a terminal as `$ zstarview` shows startup messages and errors.
Logs are also written to a file (platform‑dependent). Examples:
- Linux: `~/.cache/zstarview/logs/app.log`
- macOS: `~/Library/Logs/zstarview/app.log`
- Windows: `%LOCALAPPDATA%/tos-kamiya/zstarview/Logs/app.log`

If normal startup closes too quickly to inspect, try `zstarview-debug` from a terminal.
It is mainly for Windows troubleshooting. On Linux, `zstarview-debug` behaves the same as `zstarview`.

On Windows, Windows Security may block loading Python extension modules and the app may stop during startup.
If that happens, changing the Smart App Control setting under Windows Security `App & browser control` may help.
However, this weakens security, so it is not recommended outside a trusted environment.
See [this Smart App Control screenshot](docs/images/windows-smart-app-control.png).

</details>

<details>
  <summary>Code, Data Licenses, and Credits</summary>

## Code, Data Licenses, and Credits

This software is provided under the [MIT](LICENSE.txt) License.

However, the **included data** is redistributed according to their respective licenses.

All paths below are relative to `src/zstarview/data/`.

| File                                                           | Content                                          | Source                                                             | License                                                                                                                             |
| -------------------------------------------------------------- | ------------------------------------------------ | ------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------- |
| `cities1000.txt`, `admin1CodesASCII.txt`                       | List of cities with a population of 1000 or more | [GeoNames](https://download.geonames.org/export/dump/)             | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)                                                                           |
| `viewpoints/tower_viewpoints.json`                             | Tower/viewpoint dataset packaged for tower-name startup resolution (derived and normalized from Wikidata) | [Wikidata](https://www.wikidata.org/) via local normalization/query workflow documented in `docs/developer/viewpoint-dataset-generation.md` | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) (Wikidata data) |
| `viewpoints/mountain_viewpoints.json`                          | Mountain/viewpoint dataset packaged for mountain-name startup resolution (Wikipedia-curated candidates normalized with Wikidata metadata) | [Wikipedia](https://www.wikipedia.org/) candidate collection plus [Wikidata](https://www.wikidata.org/) normalization workflow documented in `docs/developer/viewpoint-dataset-generation.md` | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) (Wikidata data) |
| `earth_guide_land_110m.json`                                   | Simplified land geometry used to generate the below-horizon Earth guide hatch pattern, derived from Natural Earth 1:110m land polygons | [Natural Earth](https://www.naturalearthdata.com/) | [Public domain](https://www.naturalearthdata.com/about/terms-of-use/) |
| Runtime `--place` geocoding requests sent to OpenStreetMap Nominatim | Online place-name geocoding used only when `--place` is requested | [OpenStreetMap Nominatim](https://nominatim.openstreetmap.org/) | [Nominatim Usage Policy](https://operations.osmfoundation.org/policies/nominatim/) |
| Runtime IP geolocation requests sent to `ip-api.com` | IP-based location lookup used when `auto` is requested | [ip-api.com](https://ip-api.com/) | [ip-api.com Terms of Service / Privacy Policy](https://ip-api.com/docs/legal) |
| Runtime water-surface overlay data fetched via Overpass API | Water-surface points derived from OpenStreetMap inland water features for the optional river/lake/pond layer; sea-mask tiles are derived from `https://osmdata.openstreetmap.de/data/water-polygons.html` | [OpenStreetMap](https://www.openstreetmap.org/), [Overpass API](https://overpass-api.de/), [OSM Water Polygons](https://osmdata.openstreetmap.de/data/water-polygons.html) | [ODbL 1.0](https://opendatacommons.org/licenses/odbl/) |
| On-demand urban-outline cache under the app cache directory | Derived building tiles and `tile_index.json` files produced from downloaded Overture building data | [Overture Maps Buildings](https://docs.overturemaps.org/guides/buildings/) downloaded at runtime via the `overturemaps` CLI | [ODbL 1.0](https://opendatacommons.org/licenses/odbl/) |
| On-demand night lights cache under the app cache directory | NASA Earth at Night / Black Marble 2016 Grayscale 500m GeoTIFF tiles used for the optional night lights overlay | [NASA Earth at Night / Black Marble maps](https://science.nasa.gov/earth/earth-observatory/earth-at-night/maps/) | See [NASA Earth at Night / Black Marble maps](https://science.nasa.gov/earth/earth-observatory/earth-at-night/maps/) for the published NASA data use terms. |
| Runtime aircraft overlay data fetched from OpenSky Network | Aircraft state vectors used for the optional nearby-aircraft overlay | [OpenSky Network REST API](https://openskynetwork.github.io/opensky-api/rest.html) | [OpenSky Network Terms of Use](https://opensky-network.org/about/terms-of-use) |
| Runtime JPL Horizons / Small-Body Database requests | Search / ephemeris data used for celestial-body lookup and the JWST / Voyager / Parker spacecraft cache | [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/), [JPL Small-Body Database](https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html) | See the JPL/JPL SSD sites for current usage terms and data notes |
| Runtime artificial satellite overlay data fetched from `wheretheiss.at` with CelesTrak fallback | Orbital-element data used for the optional ISS overlay | [wheretheiss.at](https://wheretheiss.at/w/developer), [CelesTrak](https://celestrak.org/) | See each source site for current terms and licensing details |
| Runtime Geo-satellite cloud overlay data fetched from MET Norway | Geo-satellite cloud imagery used for the optional Europe workflow cloud overlay | [MET Norway Geo-satellite documentation](https://api.met.no/weatherapi/geosatellite/1.4/documentation), [MET Norway Licensing and crediting](https://www.met.no/en/free-meteorological-data/Licensing-and-crediting) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| `dso.csv`                                                       | Deep-sky object catalog (named galaxies/open clusters/globular clusters; generated from OpenNGC) | [OpenNGC](https://github.com/mattiaverga/OpenNGC) via [PyOngc](https://github.com/mattiaverga/PyOngc) | [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) (OpenNGC database) |
| On-demand terrain DEM cache under the app cache directory | Terrain horizon source data (Copernicus DEM GLO-90) | [Copernicus DEM / Copernicus Data Space Ecosystem](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM) via public AWS S3 distribution used by the app | Copernicus DEM GLO-90 access terms as described by Copernicus Data Space Ecosystem, including "Licence for COP-DEM-GLO-90-F Global 90m Full, Free & Open" / "Licence for the use of the Copernicus WorldDEM™-90" |
| `stars/IAU-Catalog of Star Names (always up to date).csv`      | IAU WGSN catalog of approved star names          | [exopla.net](https://exopla.net/star-names/modern-iau-star-names/) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)                                                                           |
| `Noto_Sans/*`                                                   | Font for displaying text                          | [Google Fonts](https://fonts.google.com/)                          | [SIL Open Font License 1.1](https://openfontlicense.org)                                                                            |

## Credits

* Astronomical data provided by CDS Strasbourg and the ESA Hipparcos Mission.
* City data based on GeoNames.
* Tower/viewpoint startup data are derived from Wikidata and redistributed under Wikidata's CC0 data terms.
* Mountain/viewpoint startup data are curated from Wikipedia candidates and normalized with Wikidata metadata; redistributed here under Wikidata's CC0 data terms.
* Earth guide land geometry is derived from **Natural Earth** 1:110m land polygons. Natural Earth treats the data as public domain; credit is optional, but we note the source here.
* Urban outline source data are downloaded on demand from **Overture Maps Buildings** and converted into cached derived building tiles for runtime use.
* Night lights source data are downloaded on demand from **NASA** Earth at Night / Black Marble and cached locally as GeoTIFF tiles for runtime use.
* Star proper names provided by the **IAU** Working Group on Star Names (**WGSN**) (via exopla.net).
* Cloud data are based on infrared observations from the **Himawari** satellite (provided by **JMA**) and the **NOAA GOES** series (provided by **NOAA/NESDIS**), retrieved from their public S3 buckets.
* Geo-satellite cloud imagery is provided by **MET Norway** (The Norwegian Meteorological Institute) and used under **CC BY 4.0** terms, with attribution to MET Norway as the data source.
* Aircraft overlay data are fetched from **OpenSky Network** at runtime and are subject to the OpenSky Network Terms of Use.
* Celestial-body search uses **JPL Horizons** and the **JPL Small-Body Database** at runtime to resolve matches and observer ephemerides. Search results and ephemerides are subject to the current JPL/JPL SSD usage terms and data notes.
* Orbital data (TLE/OMM) for the artificial satellite overlay are fetched from **wheretheiss.at** with **CelesTrak** as a fallback source.
* Terrain horizon data are based on **Copernicus DEM GLO-90**, managed by **ESA** on behalf of the European Commission and obtained by the app through its public AWS distribution/cache flow.
* Place/station search via `--place` uses the public **OpenStreetMap Nominatim** service and is subject to the Nominatim Usage Policy.
* Automatic IP-based location lookup uses **ip-api.com** and is subject to the ip-api.com Terms of Service / Privacy Policy, including its non-commercial-use restriction and 45 requests per minute limit.
* Water-surface overlay data for rivers, lakes, and ponds are derived from **OpenStreetMap** inland water features fetched through **Overpass API** and are attributed to **OpenStreetMap contributors** under **ODbL 1.0**.
* Sea-mask tiles for ocean water surfaces come from [OSM Water Polygons](https://osmdata.openstreetmap.de/data/water-polygons.html) and are also covered by the **OpenStreetMap** attribution / **ODbL 1.0** terms used by that dataset.
* Thanks to **Overture Maps** and its source data contributors for making large-scale building data available.
* Thanks to **AWS** and dataset providers for making the public S3 distribution/mirror endpoints available for cloud imagery and terrain DEM access.
* Fonts provided by the **Google Noto Project**.
* The window title "Zenith Star View" was suggested by **ChatGPT**.
* Specification discussions, code generation, and debugging were greatly assisted by **Gemini** and **ChatGPT**.

</details>

## Appendix

→ [Developer Docs](docs/developer/README.md)

→ [Specification](docs/specification.md), [Design](docs/design.md)

→ [Release Notes](release-notes.md)

→ [Lunar eclipses in 2026-2028, Solar eclipses 2026-2028](docs/appendix-eclipses.md)
