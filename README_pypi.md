# zstarview

Transparent desktop sky viewer with stars, planets, eclipses, optional real-time satellite cloud overlay, and optional terrain horizon overlay.

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
zstarview [options] [city]
```

Examples:

```bash
zstarview Tokyo
zstarview "35.68;139.76"
zstarview -Z E -A 25 Tokyo
```

## Highlights

- Real-time rendering of bright stars, planets, celestial equator, and ecliptic.
- Optional deep-sky object (DSO) overlays.
- Optional asterism overlays (line patterns, not IAU constellation boundaries).
- Real-time cloud overlays from Himawari/GOES (can be disabled with `-c 0`).
- Optional terrain horizon overlay from Copernicus DEM (can be disabled with `--terrain-horizon-opacity 0`).
- Hover interactions for object labels and contextual overlays.

## Useful Options

- `-V`, `--vmag-limit`: star magnitude limit.
- `--datetime "YYYY-MM-DD HH[:MM[:SS]] [TZ]"`: absolute sky time.
- `--show-dso-initial true|false`: DSO visibility at startup.
- `--show-asterisms-initial true|false`: asterism visibility at startup.
- `-c`, `--cloud-opacity`: cloud overlay opacity (`0` disables clouds).
- `--terrain-horizon-opacity`: terrain horizon opacity (`0` disables DEM download and terrain rendering).

## Notes

- On first launch, ephemeris data may be downloaded once and cached.
- Terrain horizon rendering downloads Copernicus DEM tiles on first use and then reuses the local cache.
- If cloud fetching is slow or unavailable, use `-c 0` for offline-friendly star/planet viewing.

## Links

- Homepage: https://github.com/tos-kamiya/zstarview
- Documentation: https://github.com/tos-kamiya/zstarview#readme
- Issues: https://github.com/tos-kamiya/zstarview/issues
