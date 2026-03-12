# zstarview

Transparent desktop sky viewer with stars, planets, eclipses, optional real-time satellite cloud overlay, optional terrain horizon overlay, and optional urban outline overlay.

<img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot1.png" width="640" alt="zstarview screenshot">

## Install

Recommended:

```bash
pipx install zstarview
```

Or with pip:

```bash
pip install zstarview
```

## Quick Start

```bash
zstarview [options] [location]
```

Examples:

```bash
zstarview Tokyo
zstarview "Tokyo Skytree"
zstarview "35.68;139.76"
zstarview -Z E -A 25 Tokyo
```

## Highlights

- Deep-sky objects: named galaxies/open clusters/globular clusters are shown as soft blue extents.
- Asterism overlay: popular line patterns rather than formal IAU constellation boundaries are shown as dim ambient lines.
- Solar-system bodies: supports Sun, Moon, and major planets.
- Flexible location input: specify the observer location through the CLI argument using a city name, tower name, mountain name, or direct latitude/longitude input.
- Adjustable view center: adjust the view center with CLI options `-A` and `-Z`, or with the arrow keys.
- Satellite cloud imagery: real-time Himawari/GOES satellite data are downloaded and rendered as a stylized hatched overlay.
- Terrain horizon and ground fill: Copernicus DEM data can be downloaded to render the local terrain skyline and ground region below the horizon.
- Urban outline overlay: where bundled PLATEAU-derived building tiles are available, major rooflines are drawn as a white overlay; very narrow roof spans are simplified to thick horizontal strokes.
- Never-rises region: the celestial region that never rises above the horizon for the observer's latitude is shown in a red tint.
- Python support: routinely tested on CPython 3.10, 3.11, 3.12, and 3.13.

## Common Options

- `--cloud-opacity 0.0..1.0`
- `--terrain-horizon-opacity 0.0..1.0`
- `--urban-outline-opacity 0.0..1.0`
- `--observer-height-m METERS`
- `--datetime "YYYY-MM-DD HH[:MM[:SS]] [TZ]"`

## Links

- Homepage: https://github.com/tos-kamiya/zstarview
- Documentation: https://github.com/tos-kamiya/zstarview#readme
- Issues: https://github.com/tos-kamiya/zstarview/issues
