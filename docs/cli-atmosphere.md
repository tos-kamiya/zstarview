# Atmosphere

## Atmosphere

| Option | Description | Default |
| :----- | :---------- | :------ |
| `-S`, `--sky-opacity SKY_OPACITY` | Overall intensity of the simulated sky-color disc (0.0–1.0). Higher values favor luminance over proportional saturation growth. Use 0.0 to disable the disc and its bright-body contrast underlay. | `0.3` |
| `--sky-disc-altaz-rings {off,dimalt,altaz}` | Always-visible sky-disc alt/az overlay. `dimalt` shows subtle altitude rings; `altaz` shows the full grid. | `dimalt` |
| `--sky-disc-altaz-rings-hover {off,dimalt,altaz}` | Hover-time sky-disc alt/az overlay. Same meanings as above. | `altaz` |
| `-c`, `--cloud-opacity CLOUD_OPACITY` | Opacity of cloud rendering (0.0–1.0). Use 0.0 to disable for the session, even when `--geo-satellite true` is enabled. At night, the effective opacity is smoothly raised by up to 30% according to solar altitude to keep clouds visible. \*2 | `0.3` |
| `--cloud-mode {voxel,shell}` | Use volumetric satellite-pixel voxels or project altitude shells to the display. In shell mode, `--cloud-stripe` selects the visual style. | `voxel` |
| `--geo-satellite true\|false` | Use the experimental Geo-satellite infrared cloud path inside the supported Europe workflow band. | `false` |
| `--cloud-stripe MODE[,COUNT[,WIDTH]]` | Select the appearance of shell-mode clouds; ignored in voxel mode. `halftone2` renders the current per-shell quantized round-dot appearance; `halftone` renders the legacy aggregate-cloud appearance. `width` draws centered symmetric stripes whose visible width varies continuously with cloud amount; `width-quantized` is its five-level variant with gaps at width transitions and round line caps; `alpha` keeps width fixed and varies stripe alpha. `COUNT` is the absolute number of stripes drawn across the disc; it is not scaled with window size. The mode defaults expand to `halftone,30,1.7`, `halftone2,30,1.7`, `width,50,0.85`, `width-quantized,50,0.85`, and `alpha,50,0.25`. If count or width is `0`, cloud rendering is disabled. | not selected by default |
| `--cloud-missing-tint-opacity OPACITY` | Opacity of missing-cloud-data yellow tint (0.0–1.0). | `0.176` |
| `-P`, `--precipitation-opacity OPACITY` | Opacity of the opt-in Open-Meteo forecast precipitation rain streaks (0.0–1.0). The nearest 15-minute forecast interval midpoint is selected. Outside native 15-minute model regions, values may be interpolated from hourly forecasts. A positive value requires one-time confirmation of the non-commercial Free API terms. | `0.0` |
| `--tropical-cyclone-opacity OPACITY` | Opacity of the tropical cyclone overlay (0.0–1.0). Use 0.0 to disable cyclone API fetch and drawing for that run. The overlay is also hidden automatically for time-shifted views. | `0.7` |
| `-a`, `--aircraft-opacity OPACITY` | Opacity of the aircraft overlay (0.0–1.0). Use a positive value to explicitly enable aircraft queries and drawing; 0.0 disables them. | `0.0` |
| `--satellite-opacity OPACITY` | Opacity of the artificial satellite overlay (0.0–1.0). Use 0.0 to disable satellite element fetch and drawing for that run. | `0.7` |
| `--meteor-trails-opacity OPACITY` | Opacity of GMN meteor trails (0.0–1.0). Use 0.0 to disable fetching, drawing, and menu re-enabling for that run. | `0.5` |
| `--meteor-trails-max-candidates N` | Display at most `N` of the newest GMN trails after geographic filtering. Use `0` for no limit. | `150` |
#### Footnotes

\*2 Cloud rendering uses infrared data from meteorological satellites (**Himawari** and **NOAA GOES** series), retrieved from their public S3 buckets. See Troubleshooting for tips on slow networks or offline use (for example, disabling clouds with `-c 0`). If Geo-satellite is enabled, `-c 0` still keeps cloud rendering disabled until the user re-enables it manually.
