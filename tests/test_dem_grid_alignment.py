from __future__ import annotations

import numpy as np
import rasterio
from affine import Affine

from zstarview.terrain.dem import GeoTiffDem


def test_dem_sample_uses_source_grid_for_different_bounds(tmp_path) -> None:
    path = tmp_path / "dem.tif"
    data = np.arange(100, dtype=np.float32).reshape(10, 10)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=10,
        height=10,
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=Affine(0.01, 0, 10.003, 0, -0.01, 20.007),
    ) as dataset:
        dataset.write(data, 1)

    dem = GeoTiffDem((path,))
    try:
        small = dem.build_grid((10.026, 19.941, 10.075, 19.996))
        large = dem.build_grid((10.013, 19.918, 10.094, 20.002))
        lon = np.array([10.050])
        lat = np.array([19.959])
        small_value = small.sample_lonlat(lon, lat)[0]
        large_value = large.sample_lonlat(lon, lat)[0]
        assert np.isclose(small_value, 47.2)
        assert np.isclose(large_value, small_value)
        offset = (small.transform.c - 10.003) / 0.01
        assert np.isclose(offset, round(offset), atol=1e-9)
    finally:
        dem.close()
