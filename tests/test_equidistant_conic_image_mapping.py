from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image
from pyproj import Transformer

from zstarview.utils.geostationary_image_mapping import (
    build_equidistant_conic_projection,
    fit_equidistant_conic_projection_from_image,
)


def test_equidistant_conic_fit_recovers_known_projection(tmp_path: Path) -> None:
    image = Image.new("RGB", (240, 180), color=(0, 0, 0))
    projection = build_equidistant_conic_projection(
        longitude_of_projection_origin=10.0,
        latitude_of_projection_origin=45.0,
        standard_parallel_1=30.0,
        standard_parallel_2=60.0,
    )
    to_lonlat = Transformer.from_crs(projection, "EPSG:4326", always_xy=True)
    control_pixels = [
        ("#ff0000", (30.0, 24.0)),
        ("#cc0000", (70.0, 40.0)),
        ("#aa0000", (110.0, 68.0)),
        ("#880000", (160.0, 96.0)),
        ("#660000", (190.0, 130.0)),
        ("#440000", (90.0, 150.0)),
    ]
    lines: list[str] = []
    for color, (x_px, y_px) in control_pixels:
        x_proj = (x_px - 120.0) * 8000.0
        y_proj = (90.0 - y_px) * 8000.0
        lon_deg, lat_deg = to_lonlat.transform(x_proj, y_proj)
        lines.append(f"{float(lat_deg):.10f}, {float(lon_deg):.10f} {color}")
        image.putpixel((int(x_px), int(y_px)), tuple(int(color[i : i + 2], 16) for i in (1, 3, 5)))

    map_path = tmp_path / "latlonmap.txt"
    map_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    fitted, located = fit_equidistant_conic_projection_from_image(
        image,
        map_path,
        longitude_of_projection_origin_candidates=[0.0, 10.0, 20.0],
        latitude_of_projection_origin_candidates=[35.0, 45.0, 55.0],
        standard_parallel_1_candidates=[20.0, 30.0, 40.0],
        standard_parallel_2_candidates=[50.0, 60.0, 70.0],
    )

    assert len(located) == 6
    assert np.isclose(fitted.standard_parallel_1, 30.0)
    assert np.isclose(fitted.standard_parallel_2, 60.0)
    assert fitted.rms_pixel_error < 1e-4
    assert fitted.max_pixel_error < 1e-4
