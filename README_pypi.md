# zstarview

Transparent desktop sky viewer with stars, planets, eclipses, optional real-time satellite cloud overlay, optional terrain horizon overlay, optional urban outline overlay, and an optional nearby-aircraft overlay.

## Screenshots

The first screenshot shows the asterism overlay together with the never-rises region.
The second screenshot shows the aircraft overlay.
The third screenshot shows a denser star field rendered with `-V10 -s4.5`.
The fourth screenshot shows terminal output via sixel using `zstarview-export-image`.

<p align="center">
  <img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot1.png" width="49%" alt="Screenshot showing the asterism overlay and the never-rises region">
  <img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot4.png" width="49%" alt="Screenshot showing the aircraft overlay">
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot3.png" width="49%" alt="Screenshot showing a denser star field rendered with -V10 -s4.5">
  <img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot6.png" width="49%" alt="Screenshot showing sixel terminal output from zstarview-export-image">
</p>

Note: higher magnitude limits increase rendering time. See [about magnitude limit](https://github.com/tos-kamiya/zstarview#about-magnitude-limit).

Urban outline examples from several cities worldwide:

<table>
  <tr>
    <td align="center" width="25%"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot-tokyotower.png" alt="Near Tokyo Tower, Tokyo" width="100%" /></td>
    <td align="center" width="25%"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot-downtowndubai.png" alt="Downtown Dubai" width="100%" /></td>
    <td align="center" width="25%"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot-marinabay.png" alt="Marina Bay, Singapore" width="100%" /></td>
    <td align="center" width="25%"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot-sydney.png" alt="Circular Quay, Sydney" width="100%" /></td>
  </tr>
  <tr>
    <td align="center"><sub>Near Tokyo Tower, Tokyo</sub></td>
    <td align="center"><sub>Downtown Dubai</sub></td>
    <td align="center"><sub>Marina Bay, Singapore</sub></td>
    <td align="center"><sub>Circular Quay, Sydney</sub></td>
  </tr>
</table>

## Install

Recommended:

Prerequisite for the urban outline overlay: install the `overturemaps` CLI separately.
Installation: <https://pypi.org/project/overturemaps/>

Confirm it with:

```bash
overturemaps --help
```

```bash
pipx install zstarview
```

Or with pip:

```bash
pip install zstarview
```

> Note: Windows on Arm64 is currently not supported for installation.
> As of 2026-03-15, native dependencies such as `shapely` can fail there
> because they may require a source build.

## Quick Start

```bash
zstarview [options] [location]
```

Examples:

```bash
zstarview Tokyo
zstarview "Tokyo Skytree"
zstarview "35.68;139.76"
zstarview --place "Matsue Station" --place-countrycode jp
zstarview -Z E -A 25 Tokyo
```

## Highlights

- Deep-sky objects: named galaxies/open clusters/globular clusters are shown as soft blue extents.
- Asterism overlay: popular line patterns rather than formal IAU constellation boundaries are shown as dim ambient lines.
- Solar-system bodies: supports Sun, Moon, and major planets.
- Flexible location input: specify the observer location through the CLI argument using a city name, tower name, mountain name, direct latitude/longitude input, or online place/station search via Nominatim.
- Adjustable view center: adjust the view center with CLI options `-A` and `-Z`, or with the arrow keys.
- Never-rises region: the celestial region that never rises above the horizon for the observer's latitude is shown in a red tint.
- Satellite cloud imagery: real-time Himawari/GOES satellite data are downloaded and rendered as a stylized hatched overlay.
- Terrain horizon and ground fill: Copernicus DEM data can be downloaded to render the local terrain skyline and ground region below the horizon.
- Urban outline overlay: major rooflines are drawn as a white overlay for the current viewpoint. In some skyscraper-heavy cities, distant skyscrapers can also be added from within a 10km radius.
- Aircraft overlay: nearby aircraft from OpenSky can be drawn as purple predicted-motion polylines, with startup opacity control and cached reuse when the layer is hidden and shown again.
- Python support: routinely tested on CPython 3.10, 3.11, 3.12, and 3.13.

## Common Options

- `--place QUERY`
- `--place-countrycode CODE`
- `--place-lang LANG`
- `--sky-opacity 0.0..1.0`
- `--cloud-opacity 0.0..1.0`
- `--cloud-missing-tint-opacity 0.0..1.0`
- `--terrain-horizon-opacity 0.0..1.0`
- `--ground-tint-opacity 0.0..1.0`
- `--urban-outline-opacity 0.0..1.0`
- `-a, --aircraft-opacity 0.0..1.0`
- `--observer-height-m METERS`
- `--datetime "YYYY-MM-DD HH[:MM[:SS]] [TZ]"`

Notes:

- `--place` uses the public OpenStreetMap Nominatim search service and sends a single request with a User-Agent and `Accept-Language`.
- Satellite cloud rendering downloads Himawari/GOES data from public S3 buckets.
- Terrain horizon rendering downloads Copernicus DEM tiles on first use and reuses cached data later.
- Aircraft rendering uses OpenSky state data when enabled; `-a 0` disables both aircraft queries and drawing for that run.

## Code, Data Licenses, and Credits

- Code: MIT License. See `LICENSE`.
- Bundled and runtime-fetched data may be subject to their own licenses, attribution rules, or service terms.
- See the main project README for the full credits and third-party data notes.

## Links

- Homepage: https://github.com/tos-kamiya/zstarview
- Documentation: https://github.com/tos-kamiya/zstarview#readme
- Issues: https://github.com/tos-kamiya/zstarview/issues
