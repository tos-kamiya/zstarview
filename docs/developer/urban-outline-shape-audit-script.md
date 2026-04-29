# Urban Outline Shape Audit Script

`scripts/audit_urban_outline_shapes.py` scans derived building tiles and reports
how much of the dataset is made up of small footprint geometries.

## What It Measures

The script measures each building footprint in meters by projecting the ring
coordinates into a local tangent plane centered on that footprint. For each
building it reports:

- the number of rings
- the number of points
- the footprint bounding-box width and height in meters
- the larger of those two dimensions
- the approximate area and perimeter

The `max_dim_m` value is the most useful quick indicator for "small shapes".

## Usage

```bash
python scripts/audit_urban_outline_shapes.py \
  --derived-dir ~/.cache/zstarview/overture_buildings \
  --threshold-m 5 \
  --threshold-m 10 \
  --threshold-m 20 \
  --threshold-m 50
```

By default the script scans both derived cache roots when they exist:

- `overture_buildings`
- `overture_skyscrapers`

Use `--json` when you want machine-readable summary output for follow-up
analysis.

Use `--csv-out <path>` to write a histogram CSV of `max_dim_m` values. The CSV
uses one row per meter bin by default:

- `bin_start_m`
- `bin_end_m`
- `count`

Set `--bin-size-m` if you want a different bin width. The histogram uses the
raw `max_dim_m` values; height-based widening is handled by the renderer, not
by this audit script.

## Reading the Output

- `min_max_dim_m` shows the smallest footprint in the scanned dataset.
- `p50_max_dim_m` is the median footprint size.
- `threshold_counts` tells you how many buildings fall below each size limit.
- The `smallest` section lists the smallest buildings, along with their source
  file and footprint metrics.

The CSV histogram is the easiest format for plotting a distribution chart in a
spreadsheet or notebook.

The script is intended for exploratory auditing, not for generating committed
assets.
