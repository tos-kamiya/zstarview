from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

from zstarview.utils.geostationary_grid_overlay import (
    draw_geostationary_latlon_grid,
    load_lonlat_grid_from_npz,
)


def test_draw_geostationary_latlon_grid_colors_major_and_minor_lines(tmp_path: Path) -> None:
    width = 61
    height = 61
    image = Image.new("RGB", (width, height), color=(255, 255, 255))
    lon_row = np.linspace(-30.0, 30.0, width, dtype=np.float64)
    lat_col = np.linspace(30.0, -30.0, height, dtype=np.float64)
    lon_deg = np.tile(lon_row, (height, 1))
    lat_deg = np.tile(lat_col[:, None], (1, width))

    overlay = draw_geostationary_latlon_grid(image, lon_deg, lat_deg, step_deg=10, major_step_deg=30)
    pixels = np.asarray(overlay.convert("RGB"), dtype=np.uint8)

    assert tuple(pixels[30, 30]) == (255, 0, 0)
    assert tuple(pixels[20, 40]) == (0, 0, 0)
    assert tuple(pixels[40, 20]) == (0, 0, 0)

    npz_path = tmp_path / "grid.npz"
    np.savez_compressed(npz_path, lon_deg=lon_deg, lat_deg=lat_deg)
    loaded = load_lonlat_grid_from_npz(npz_path)
    assert loaded.lon_deg.shape == (height, width)
    assert loaded.lat_deg.shape == (height, width)
