# Atmosphere

## Atmosphere

| Option | Description | Default |
| :----- | :---------- | :------ |
| `-S`, `--sky-opacity SKY_OPACITY` | Overall intensity of the simulated sky-color disc (0.0–1.0). Higher values favor luminance over proportional saturation growth. Use 0.0 to disable the disc and its bright-body contrast underlay. | `0.3` |
| `--sky-disc-altaz-rings {off,dimalt,altaz}` | Always-visible sky-disc alt/az overlay. `dimalt` shows subtle altitude rings; `altaz` shows the full grid. | `dimalt` |
| `--sky-disc-altaz-rings-hover {off,dimalt,altaz}` | Hover-time sky-disc alt/az overlay. Same meanings as above. | `altaz` |
| `-c`, `--cloud-opacity CLOUD_OPACITY` | Opacity of cloud rendering (0.0–1.0). Use 0.0 to disable for the session, even when `--geo-satellite true` is enabled. \*2 | `0.175` |
| `--geo-satellite true\|false` | Use the experimental MET Norway infrared image source inside the supported Europe workflow band. Voxel clouds estimate a 1–9 km layer allocation from the display image. | `false` |
| `--cloud-stripe holes[,RADIUS]\|cutout[,COUNT]` | Transparent marks that make clouds easier to recognize. `holes` cuts small round holes at 126 points fixed in alt/az, on concentric equal-altitude rings around the zenith (a zenith point plus rings 12 degrees apart, with points about 13 to 15 degrees apart along each ring; each ring is rotated by the golden angle so the rings form sunflower-like spiral arms); `RADIUS` is in pixels of the 513x513 source raster, so the on-screen size does not change with the field of view. `cutout` draws `COUNT` (default 16) transparent lines that slope down to the right at 45 degrees and are 1 px wide on the source raster. A value of `0` disables only the marks. The former shell-cloud styles (`halftone2`, `halftone`, `width`, `alpha`) were removed and are rejected. | `holes,2` |
| `-P`, `--precipitation-opacity OPACITY` | Opacity of the opt-in Open-Meteo forecast precipitation rain streaks (0.0–1.0). The nearest 15-minute forecast interval midpoint is selected. Outside native 15-minute model regions, values may be interpolated from hourly forecasts. A positive value requires one-time confirmation of the non-commercial Free API terms. | `0.0` |
| `--tropical-cyclone-opacity OPACITY` | Opacity of the tropical cyclone overlay (0.0–1.0). Use 0.0 to disable cyclone API fetch and drawing for that run. The overlay is also hidden automatically for time-shifted views. | `0.7` |
| `-a`, `--aircraft-opacity OPACITY` | Opacity of the aircraft overlay (0.0–1.0). Use a positive value to explicitly enable aircraft queries and drawing; 0.0 disables them. | `0.0` |
| `--satellite-opacity OPACITY` | Opacity of the artificial satellite overlay (0.0–1.0). Use 0.0 to disable satellite element fetch and drawing for that run. | `0.7` |
| `--meteor-trails-opacity OPACITY` | Opacity of GMN meteor trails (0.0–1.0). Use 0.0 to disable fetching, drawing, and menu re-enabling for that run. | `0.5` |
| `--meteor-trails-max-candidates N` | Display at most `N` of the newest GMN trails after geographic filtering. Use `0` for no limit. | `150` |
#### Footnotes

\*2 Cloud rendering uses infrared data from meteorological satellites (**Himawari** and **NOAA GOES** series), retrieved from their public S3 buckets. See Troubleshooting for tips on slow networks or offline use (for example, disabling clouds with `-c 0`). If Geo-satellite is enabled, `-c 0` still keeps cloud rendering disabled until the user re-enables it manually.
