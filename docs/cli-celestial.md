# Celestial

| Option | Description | Default |
| :----- | :---------- | :------ |
| `-V`, `--vmag-limit V_MAG_LIMIT` | Maximum visual magnitude of stars to display. | `7.0` |
| `--vmag-brightness-multiplier MULTIPLIER` | Brightness multiplier per magnitude step (allowed range 1.58–2.512, default `2.5`; 2.512 is the historical Pogson ratio). \*3 | `2.5` |
| `-m`, `--enlarge-moon` | Compatibility shortcut for `--moon-style sphere --moon-scale 5`. It cannot be combined with `--moon-style` or `--moon-scale`. | |
| `--moon-style {marker,sphere,image}` | Moon rendering style. `marker` draws a compact phase-aware outline, `sphere` draws a procedural Lambert-shaded phase, and `image` uses NASA Dial-A-Moon with a flat procedural phase fallback while the image is unavailable. | `marker` |
| `--moon-scale {1,2,3,4,5,6,7,8}` | Integer Moon display scale. The selected scale applies to `marker`, `sphere`, and `image`. | `1` |
| `--bright-bodies {outline,fill}` | Rendering mode for bright stars, the Sun, and planets. `outline` emphasizes their outlines; `fill` keeps filled rendering. This option does not affect the Moon. | `outline` |
| `-s`, `--star-base-radius STAR_BASE_RADIUS` | Base size of 2nd-magnitude stars. | `4.0` |
| `-w`, `--expected-render-width EXPECTED_RENDER_WIDTH` | Expected window width for full-resolution star rendering. When celestial-disc width exceeds this, star rendering uses square-root scaling. | `600` |
| `--show-dso-initial true\|false` | Whether DSO overlays are shown at startup. | auto (`show` when catalog is available) |
| `--show-asterisms-initial true\|false` | Whether asterism overlays are shown at startup. | `show` |
| `--asterism-opacity OPACITY` | Absolute opacity of normal asterism lines (0.0–0.5). Use 0.0 to hide them. Hover emphasis remains visible and is not dimmed below the normal lines. | `0.1` (Atlas keeps its theme default) |
| `--diffuse-sky-source {akari,gaia}` | Select the diffuse all-sky layer. Gaia uses the bundled EDR3 integrated brightness/colour texture; AKARI uses the prepared far-infrared cache. | `gaia` |
| `--diffuse-sky-opacity OPACITY` | Opacity of the selected diffuse sky layer (0.0–1.0). Use 0.0 to disable it. | `0.15` for Gaia and AKARI |
| `--akari-ir-bands-opacity OPACITY` | Legacy shortcut for `--diffuse-sky-source akari --diffuse-sky-opacity OPACITY`. It cannot be combined with either `--diffuse-sky-*` option. | — |
| `--twinkle-count N` | Number of star-twinkle candidates selected per 2-second update. Use `0` to disable twinkle. This option is available only in the normal GUI. | `30` |
| `--show-guidelines-initial true\|false` | Whether geometric and celestial guideline overlays are shown at startup. | `show` |
| `-i`, `--sky-update-interval SECONDS` | Interval for updating the star/DSO sky snapshot in seconds. The sky-colour disc is independent: it updates every 15 seconds when Sun altitude is between `+15` and `-15` degrees, and every 60 seconds otherwise. Its calculation image is normally quarter-width/height, with an additional square-root reduction above a 1920-pixel disc width. | `60` |

#### Overlay visibility at startup

Use these options to control initial overlay states without post-launch menu operations:

```bash
# Start with DSO hidden and asterisms visible
zstarview --show-dso-initial false --show-asterisms-initial true Tokyo
```

#### Moon display options

Moon style and scale are independent. For example:

```bash
# A procedural shaded Moon at 3x
zstarview --moon-style sphere --moon-scale 3 Tokyo

# A Dial-A-Moon image at 5x, with a procedural fallback
zstarview --moon-style image --moon-scale 5 Tokyo
```

In the GUI, the M key and **Moon Option** menu item temporarily toggle the
configured Moon display. With the default `marker` at 1x, they switch to a 5x
procedural sphere. With any non-default style or scale, they switch that
configuration off to `marker` at 1x and restore it on the next toggle.

#### About magnitude limit

Use `-V magnitude` to limit the displayed stars to those brighter than the given magnitude.
The default is `-V 7.0`. The bundled catalog currently supports up to `-V 10.5`, for which the candidate star set contains about 536,000 stars.
Note that higher values will increase rendering time.

### Supported Asterisms

These overlays are **asterisms** (popular line patterns), not formal IAU constellation boundaries.

- Winter: `Winter Triangle`, `Orion's Belt`, `Winter Hexagon`, `Southern Cross`, `Southern Pointers`, `Diamond Cross`, `False Cross`
- Spring: `Big Dipper`, `Little Dipper`, `Spring Triangle`, `Arc to Arcturus`, `Leo Sickle`, `Southern Triangle`
- Summer: `Summer Triangle`, `Northern Cross`, `Teapot`, `Keystone`
- Autumn: `Great Square of Pegasus`, `Circlet of Pisces`, `Water Jar of Aquarius`, `Cassiopeia W`, `House of Cepheus`, `Job's Coffin`

#### Footnotes

\*3 The brightest-magnitude multiplier cannot exceed the classical Pogson value of \(100^{1/5}\approx2.512\).
