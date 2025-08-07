# zstarview 🌌

See the starry sky, even when it's cloudy or the sun is out.

Zenith Star View is an application that displays the starry sky from any city on Earth.

- Renders bright stars, planets, the celestial equator, and the ecliptic in real-time.
- Specify location by city name (from GeoNames).

## Installation (Recommended: `pipx`)

It is intended to be installed using [`pipx`](https://pypa.github.io/pipx/).

```bash
pipx install git+https://github.com/yourusername/zstarview.git
```

## Usage

```bash
zstarview [options] city
```

### Argument

| Argument | Description | Default |
|:---|:---|:---|
| `city` | Specify the city name to display. | `Tokyo` |

### Options

| Option | Description | Default |
|:---|:---|:---|
| `-h`, `--help` | Show this help message and exit. | |
| `-H`, `--hours HOURS` | Number of hours to add to the current time. | `0` |
| `-D`, `--days DAYS` | Number of days to add to the current time. | `0` |
| `-m`, `--enlarge-moon` | Show the moon in 3x size. | |
| `-s`, `--star-base-radius STAR_BASE_RADIUS` | Base size of stars. | `15.0` |
| `-Z`, `--view-center-az VIEW_CENTER_AZ` | Viewing azimuth angle [deg] (0=N, 90=E, 180=S, 270=W). | `180` |
| `-A`, `--view-center-alt VIEW_CENTER_ALT` | Viewing altitude angle [deg] (90=zenith, 0=horizon). | `90` |

By default, the view shows the sky from Tokyo (`--city Tokyo`), looking south (`-Z 180`) towards the zenith (`-A 90`). In this view, the bottom of the screen is South, and the left side is East.

### Key Operations

*   **F11**: Toggle fullscreen display
*   **ESC**: Exit fullscreen
*   **Q**: Quit

## Dependencies

*   Python 3.10+
*   [`appdirs`](https://pypi.org/project/appdirs/)
*   [`PyQt5`](https://pypi.org/project/PyQt5/)
*   [`astropy`](https://pypi.org/project/astropy/)
*   [`skyfield`](https://pypi.org/project/skyfield/)
*   [`numpy`](https://pypi.org/project/numpy/)
*   [`Pillow`](https://pypi.org/project/Pillow/)

## License

This software is provided under the [MIT](LICENSE.txt) License.

However, the **included data** is redistributed according to their respective licenses.

| File | Content | Source | License |
|---|---|---|---|
| `data/cities1000.txt` | List of cities with a population of 1000 or more | [GeoNames](https://download.geonames.org/export/dump/) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| `data/catalog` | Hipparcos Star Catalog | [CDS Strasbourg](https://cdsarc.cds.unistra.fr/viz-bin/cat/V/50) | [ODbL](https://www.data.gouv.fr/licences) or [CC BY-NC 3.0 IGO](https://creativecommons.org/licenses/by-nc/3.0/igo/) Non-commercial |
| `data/NotoSansSymbols-VariableFont_wght.ttf` | Font for displaying planetary symbols | [Google Fonts](https://fonts.google.com/noto/specimen/Noto+Sans+Symbols) | [SIL Open Font License 1.1](https://openfontlicense.org) |

## Credits

*   Astronomical data provided by CDS Strasbourg and the ESA Hipparcos Mission.
*   City data based on GeoNames.
*   Fonts provided by the Google Noto Project.
*   The window title "Zenith Star View" was suggested by ChatGPT.
*   Specification discussions, code generation, and debugging were greatly assisted by Gemini and ChatGPT.
