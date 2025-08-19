# zstarview 🌌

See the starry sky, even when it's cloudy or the sun is out.

Zenith Star View is an application that displays the starry sky from any city on Earth.

- Renders bright stars, planets, the celestial equator, and the ecliptic in real-time.
- Specify location by city name (from GeoNames).

![](docs/images/screenshot1.png)

- With the option `-A altitude`, you can also view the horizon-level sky
- To make the sky situation easier to understand, a simulated sky color is lightly overlaid (since version 0.8.2)

![](docs/images/screenshot4.png)

## Installation (Recommended: `pipx`)

It is intended to be installed using [`pipx`](https://pypa.github.io/pipx/).

```bash
pipx install git+https://github.com/tos-kamiya/zstarview.git
```

> Note (X11 on Ubuntu/Debian): On X11 sessions, Qt's xcb platform plugin may require `libxcb-cursor0` at runtime. Install it with: `sudo apt install libxcb-cursor0`.

## Usage

```bash
zstarview [options] [city]
```

### Argument

| Argument | Description                                                                                                                          | Default                            |
| :------- | :----------------------------------------------------------------------------------------------------------------------------------- | :--------------------------------- |
| `city`   | Specify the city name to display. If omitted, the city from the **last run** will be used. On the first run, it defaults to `Tokyo`. | Last run (or `Tokyo` on first run) |

### Options

### Options

| Option                                      | Description                                             | Default |
| :------------------------------------------ | :------------------------------------------------------ | :------ |
| `-h`, `--help`                              | Show this help message and exit.                        |         |
| `-H`, `--hours HOURS`                       | Number of hours to add to the current time.             | `0`     |
| `-D`, `--days DAYS`                         | Number of days to add to the current time.              | `0`     |
| `--datetime "YYYY-MM-DD HH:MM:SS [TZ]"`     | Specify an absolute date/time (UTC if no TZ given).     |         |
| `-m`, `--enlarge-moon`                      | Show the moon in 5x size.                               |         |
| `-s`, `--star-base-radius STAR_BASE_RADIUS` | Base size of stars.                                     | `8.0`   |
| `-Z`, `--view-center-az VIEW_CENTER_AZ`     | Viewing azimuth (degrees or compass points).            | `180`   |
| `-A`, `--view-center-alt VIEW_CENTER_ALT`   | Viewing altitude angle (90=zenith, 0=horizon).          | `90`    |
| `-V`, `--vmag-threshold V_MAG_THRESHOLD`    | Maximum visual magnitude of stars to display.           | `6.0`   |
| `--sky-opacity SKY_OPACITY`                 | Opacity of the simulated sky-color disc (0.0–1.0). Use 0.0 to disable sky-color rendering. | `0.2` |


**About the datetime option**

Use `--datetime "YYYY-MM-DD HH:MM:SS [TZ]"` to specify an absolute date and time.  
If no timezone (TZ) is specified, UTC is assumed.

You can specify the timezone as:
- A common abbreviation (JST, UTC, GMT, KST, HKT, AWST, ACST, AEST, NZST, NZDT, MSK, EAT)
- A full IANA timezone name (e.g., `Asia/Tokyo`, `Europe/Moscow`)
- A UTC offset (e.g., `UTC+9`, `UTC-07:30`)

Example:

```bash
zstarview --datetime "2025-08-17 21:00:00 JST" Tokyo
```

**About the view center options**

The `-Z` (azimuth) and `-A` (altitude) options specify the center of the displayed sky.

By default, `-Z 180` (facing south) and `-A 90` (zenith) are used.
In this view, the bottom of the screen is south, the left side is east, and the display is a circular view looking straight up toward the zenith.

For example, setting `-Z 90` (facing east) and `-A 10` (altitude 10°, i.e., looking 10° above the horizon)
will produce a roughly semicircular sky view.
→ This will capture the eastern sky showing the [Summer Triangle (Vega, Altair, Deneb)](docs/images/screenshot2.png).

Azimuth can be given in degrees or compass points (case-insensitive).
Examples: `-Z E`, `-Z ne`, `-Z SSW` (202.5°).
(Compass mapping: 0=N, 90=E, 180=S, 270=W; accepts N, NNE, NE, ENE, E, ESE, SE, SSE, S, SSW, SW, WSW, W, WNW, NW, NNW.)

**About magnitude threshold**

Use `-V magnitude` to limit the displayed stars to those brighter than the given magnitude.
The default is `-V 6.0`. For example, specifying 9.0 will display about 83k stars.
Note that higher values will increase rendering time.

→ Example: display up to [magnitude 9.0](docs/images/screenshot3.png).

### Key Operations

* **← / →**: Rotate view azimuth by ±5°
* **M**: Toggle moon enlarged to 5x size
* **F11**: Toggle fullscreen display
* **ESC**: Exit fullscreen
* **Q**: Quit

## Generating a `.desktop` launcher (GNOME only)

On GNOME-based environments (including Ubuntu Dock and DockToPanel),  
a `.desktop` file is required for the correct icon to appear in the taskbar.

This application includes a helper command to generate it:

```bash
# Create zstarview.desktop in the current directory
zstarview-make-desktop-file

# Install to ~/.local/share/applications
zstarview-make-desktop-file --write
```

* Without `--write`, the file `zstarview.desktop` is created in the current directory.
* With `--write`, it is installed to `~/.local/share/applications` and registered with the desktop database.

> **Note:** This launcher integration is only intended for GNOME-based environments.
> It is not required on other desktop environments, and may not work as intended elsewhere.

## License

This software is provided under the [MIT](LICENSE.txt) License.

However, the **included data** is redistributed according to their respective licenses.

| File                                         | Content                                          | Source                                                                   | License                                                                                                                             |
| -------------------------------------------- | ------------------------------------------------ | ------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------- |
| `data/cities1000.txt`                        | List of cities with a population of 1000 or more | [GeoNames](https://download.geonames.org/export/dump/)                   | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)                                                                           |
| `data/stars/hip_main.dat`                               | Hipparcos and Tycho Catalogues (ESA 1997)        | [CDS Strasbourg](https://cdsarc.cds.unistra.fr/ftp/I/239/)               | [ODbL](https://www.data.gouv.fr/licences) or [CC BY-NC 3.0 IGO](https://creativecommons.org/licenses/by-nc/3.0/igo/) Non-commercial |
| `data/stars/IAU-Catalog-of-Star-Names.csv`         | IAU Working Group on Star Names (WGSN) catalog of approved star names | [exopla.net](https://exopla.net/star-names/modern-iau-star-names/)       | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)                                                                           |
|  `data/Noto_Sans/*`, `data/Noto_Sans_Symbols/*` | Font for displaying text / planetary symbols            | [Google Fonts](https://fonts.google.com/) | [SIL Open Font License 1.1](https://openfontlicense.org)                                                                            |

## Credits

* Astronomical data provided by CDS Strasbourg and the ESA Hipparcos Mission.
* City data based on GeoNames.
* Star proper names provided by the IAU Working Group on Star Names (via [exopla.net](https://exopla.net/star-names/modern-iau-star-names/)).
* Fonts provided by the Google Noto Project.
* The window title "Zenith Star View" was suggested by ChatGPT.
* Specification discussions, code generation, and debugging were greatly assisted by Gemini and ChatGPT.

## Appendix: Lunar eclipses in 2025

You can preview the two total lunar eclipses of 2025 by specifying the datetime and city.

### March 13-14, 2025 (Total Lunar Eclipse, visible in the Americas)

```bash
zstarview --datetime "2025-03-14 02:58:43 America/New_York" "New York City"
zstarview --datetime "2025-03-13 23:58:43 America/Los_Angeles" "US/Los Angeles"
zstarview --datetime "2025-03-13 20:58:43 HST" "Honolulu"
zstarview --datetime "2025-03-14 00:58:43 America/Mexico_City" "Mexico City"
zstarview --datetime "2025-03-14 03:58:43 America/Sao_Paulo" "BR/São Paulo"
```

### September 7-8, 2025 (Total Lunar Eclipse, visible in Asia, Europe, Africa, Oceania)

```bash
zstarview --datetime "2025-09-08 03:12:00 JST" "Tokyo"
zstarview --datetime "2025-09-08 03:11:47 KST" "Seoul"
zstarview --datetime "2025-09-08 02:11:47 Asia/Shanghai" "Beijing"
zstarview --datetime "2025-09-08 01:11:47 Asia/Bangkok" "Bangkok"
zstarview --datetime "2025-09-07 23:41:47 Asia/Kolkata" "New Delhi"
zstarview --datetime "2025-09-07 21:11:47 Europe/Istanbul" "Istanbul"
zstarview --datetime "2025-09-07 21:11:47 EAT" "Nairobi"
zstarview --datetime "2025-09-08 04:11:47 Australia/Sydney" "Sydney"
zstarview --datetime "2025-09-08 06:11:47 Pacific/Auckland" "Auckland"
```
