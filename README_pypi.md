# zstarview

**Zenith Star View** is a desktop sky viewer for your chosen location.

It renders an all-sky view with stars, the Sun, Moon, planets, deep-sky objects, and guide overlays.
When enabled, it can also add real-time cloud imagery, terrain horizon, urban outlines, nearby aircraft, and the ISS/JWST/Voyager 1/Voyager 2/Parker artificial satellite overlays.
Locations can be set by city or viewpoint name, direct coordinates, online place search, or supported Google Maps URLs.

## Screenshots

The first screenshot shows the asterism overlay and also serves as a terrain-horizon example.
The second screenshot shows the aircraft overlay together with the never-rises region.
The third screenshot shows a denser star field rendered with `-V10.5 -s4.5`.
The fourth screenshot shows `zstarview-export-image` searching for and displaying `PANSTARRS (C/2026 R3)` in a sixel terminal.

<p align="center">
  <img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot1.png" width="49%" alt="Screenshot showing the asterism overlay and a terrain horizon example">
  <img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot4.png" width="49%" alt="Screenshot showing the aircraft overlay and the never-rises region">
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot3.png" width="49%" alt="Screenshot showing a denser star field rendered with -V10.5 -s4.5">
  <img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot6.png" width="49%" alt="Screenshot showing zstarview-export-image searching for and displaying PANSTARRS (C/2026 R3) in a sixel terminal">
</p>

Note: higher magnitude limits increase rendering time. See [about magnitude limit](https://github.com/tos-kamiya/zstarview#about-magnitude-limit).

Urban outline examples from several cities worldwide:

<table>
  <tr>
    <td align="center" width="25%"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot-tokyotower.png" alt="Near Tokyo Tower, Tokyo" width="100%" /></td>
    <td align="center" width="25%"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot-burjkhalifa.png" alt="Near Burj Khalifa, Dubai" width="100%" /></td>
    <td align="center" width="25%"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot-marinabay.png" alt="Marina Bay, Singapore" width="100%" /></td>
    <td align="center" width="25%"><img src="https://raw.githubusercontent.com/tos-kamiya/zstarview/main/docs/images/screenshot-sydney.png" alt="Circular Quay, Sydney" width="100%" /></td>
  </tr>
  <tr>
    <td align="center"><sub>Near Tokyo Tower, Tokyo</sub></td>
    <td align="center"><sub>Near Burj Khalifa, Dubai</sub></td>
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

> Note: The previous Windows on Arm64 installation blocker has been removed.
> Installation is now possible there, but Windows Security may still block
> Python extension modules during startup on some systems. If that happens,
> see the troubleshooting notes in the main project README.

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

- Stars: the sky view shows stars from the selected catalog, with asterisms and other overlays layered on top.
- Solar-system bodies: supports Sun, Moon, and major planets.
- Deep-sky objects: named galaxies/open clusters/globular clusters are shown as soft blue extents.
- Asterism overlay: popular line patterns are shown as dim guide lines, alongside other sky guides such as the never-rises region.
- Satellite cloud imagery and sky-color disc: real-time Himawari/GOES satellite data are downloaded and rendered as a stylized hatched overlay, with the sky-color disc still visible beneath the clouds.
- Aircraft and artificial satellite overlays: nearby aircraft from OpenSky can be drawn on the sky view, and ISS, JWST, Voyager 1, Voyager 2, and Parker can be drawn as small purple markers between the planet and aircraft layers. ISS uses `wheretheiss.at` with CelesTrak fallback, while the other four spacecraft use JPL Horizons.
- Fresh current caches are reused for up to 24 hours for both the ISS cache and the Horizons-backed spacecraft cache.
- Terrain horizon and earth guide: Copernicus DEM data can be used to render the local terrain skyline and ground region, with a separate below-horizon continental guide layer for orientation.
- Urban outline overlay: major rooflines are drawn for the current viewpoint, with optional distant skyscrapers in dense urban areas.
- Flexible location input: start from a city, tower, mountain, lat/lon, or online place/station search.
- Adjustable view center: change the view center from the CLI or with the arrow keys.
- Python support: routinely tested on CPython 3.10, 3.11, 3.12, 3.13, and 3.14.

## Common Options

- `--place QUERY`
- `--place-countrycode CODE`
- `--place-lang LANG`
- `--datetime "YYYY-MM-DD HH[:MM[:SS]] [TZ]"`
- `-Z, --view-center-az VIEW_CENTER_AZ`
- `-A, --view-center-alt VIEW_CENTER_ALT`
- `--observer-height-m METERS`
- `-V, --vmag-limit V_MAG_LIMIT`
- `--theme {night,day,white,black}`
- `-o, --output PATH` for `zstarview-export-image`

Notes:

- `--place` uses the public OpenStreetMap Nominatim search service and sends a single request with a User-Agent and `Accept-Language`.
- Satellite cloud rendering downloads Himawari/GOES data from public S3 buckets.
- Terrain horizon rendering downloads Copernicus DEM tiles on first use and reuses cached data later.
- Detailed layer-tuning options such as per-layer opacity remain available in the main README.

## Code, Data Licenses, and Credits

- Code: MIT License. See `LICENSE`.
- Bundled and runtime-fetched data may be subject to their own licenses, attribution rules, or service terms.
- See the main project README for the full credits and third-party data notes.

## Links

- Homepage: https://github.com/tos-kamiya/zstarview
- Documentation: https://github.com/tos-kamiya/zstarview#readme
- Issues: https://github.com/tos-kamiya/zstarview/issues
