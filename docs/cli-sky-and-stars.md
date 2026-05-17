# Sky and Stars

| Option | Description | Default |
| :----- | :---------- | :------ |
| `--sky-opacity SKY_OPACITY` | Opacity of the simulated sky-color disc (0.0–1.0). Use 0.0 to disable. | `0.1` |
| `--sky-disc-altaz-rings {off,dimalt,altaz}` | Always-visible sky-disc alt/az overlay. `dimalt` shows subtle altitude rings; `altaz` shows the full grid. | `dimalt` |
| `--sky-disc-altaz-rings-hover {off,dimalt,altaz}` | Hover-time sky-disc alt/az overlay. Same meanings as above. | `altaz` |
| `--bright-bodies {outline,fill}` | Bright bodies rendering mode. `outline` renders bright stars as diamond outlines, planets as outlines, and the Moon as outline-only except for enlarged moon / hover moon views. `fill` keeps the normal filled rendering. | `outline` |
| `-m`, `--enlarge-moon` | Show the moon in 5x size. | |
| `-s`, `--star-base-radius STAR_BASE_RADIUS` | Base size of 2nd-magnitude stars. | `4.0` |
| `-w`, `--expected-render-width EXPECTED_RENDER_WIDTH` | Expected window width for full-resolution star rendering. When celestial-disc width exceeds this, star rendering uses square-root scaling. | `600` |
| `-V`, `--vmag-limit V_MAG_LIMIT` | Maximum visual magnitude of stars to display. | `7.0` |
| `--vmag-brightness-multiplier MULTIPLIER` | Brightness multiplier per magnitude step (allowed range 1.58–2.512, default `2.5`; 2.512 is the historical Pogson ratio). \*3 | `2.5` |
| `-i`, `--sky-update-interval SECONDS` | Interval for updating stars/sky-color disc in seconds. | `60` |
| `--show-dso-initial true\|false` | Whether DSO overlays are shown at startup. | auto (`show` when catalog is available) |
| `--show-asterisms-initial true\|false` | Whether asterism overlays are shown at startup. | `show` |

#### Overlay visibility at startup

Use these options to control initial overlay states without post-launch menu operations:

```bash
# Start with DSO hidden and asterisms visible
zstarview --show-dso-initial false --show-asterisms-initial true Tokyo
```

#### About the view center options

The `-Z` (azimuth) and `-A` (altitude) options specify the center of the displayed sky.

By default, `-Z 180` (facing south) and `-A 90` (zenith) are used.
In this view, the bottom of the screen is south, the left side is east, and the display is a circular view looking straight up toward the zenith.

For example, setting `-Z 90` (facing east) and `-A 25` (altitude 25° above the horizon) produces a sky view toward the eastern sky.

Azimuth can be given in degrees or compass points (case-insensitive).
Examples: `-Z E`, `-Z ne`, `-Z SSW` (202.5°).
(Compass mapping: 0=N, 90=E, 180=S, 270=W; accepts N, NNE, NE, ENE, E, ESE, SE, SSE, S, SSW, SW, WSW, W, WNW, NW, NNW.)

While resizing the window, the same simplified viewport-interaction mode is used so that the view stays responsive.

#### About magnitude limit

Use `-V magnitude` to limit the displayed stars to those brighter than the given magnitude.
The default is `-V 7.0`. The bundled catalog currently supports up to `-V 10.5`, for which the candidate star set contains about 536,000 stars.
Note that higher values will increase rendering time.

#### About theme presets

Use `--theme` to change the background treatment and contrast style.

* `night`: default dark theme
* `black`: darker opaque background
* `transparent`: dark translucent background for already-dark desktops; use a flat, uniform alpha
* `day`: bright sky/background treatment
* `white`: brightest light theme

### Supported Asterisms

These overlays are **asterisms** (popular line patterns), not formal IAU constellation boundaries.

- Winter: `Winter Triangle`, `Orion's Belt`, `Winter Hexagon`, `Southern Cross`, `Southern Pointers`, `Diamond Cross`, `False Cross`
- Spring: `Big Dipper`, `Little Dipper`, `Spring Triangle`, `Arc to Arcturus`, `Leo Sickle`, `Southern Triangle`
- Summer: `Summer Triangle`, `Northern Cross`, `Teapot`, `Keystone`
- Autumn: `Great Square of Pegasus`, `Circlet of Pisces`, `Water Jar of Aquarius`, `Cassiopeia W`, `House of Cepheus`, `Job's Coffin`

#### Footnotes

\*3 The brightest-magnitude multiplier cannot exceed the classical Pogson value of \(100^{1/5}\approx2.512\).
