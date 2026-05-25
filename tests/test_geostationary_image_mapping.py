from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from PIL import Image
from pyproj import Transformer

from zstarview.utils.geostationary_image_mapping import (
    build_projection,
    estimate_pixel_to_lonlat,
    fit_pixel_lonlat_polynomial_from_image,
    fit_pixel_lonlat_tps_from_image,
    fit_pixel_map_from_image,
)

pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:You will likely lose important projection information when converting to a PROJ string:UserWarning"
    ),
]


def _write_control_fixture(tmp_path: Path) -> tuple[Path, Path, dict[str, tuple[float, float]]]:
    image = Image.new("RGB", (400, 300), color=(0, 0, 0))
    projection = build_projection(longitude_of_projection_origin=0.0, perspective_point_height=35785831.0, sweep_angle_axis="x")
    to_lonlat = Transformer.from_crs(projection, "EPSG:4326", always_xy=True)

    # Affine map from projected meters to pixels:
    #   x_px = 120 + x_proj / 12000
    #   y_px = 80 - y_proj / 10000
    control_pixels = {
        "#ff0000": (60.0, 50.0),
        "#cc0000": (140.0, 70.0),
        "#aa0000": (220.0, 110.0),
        "#880000": (300.0, 140.0),
    }
    expected = {
        "#ff0000": (60.0, 50.0),
        "#cc0000": (140.0, 70.0),
        "#aa0000": (220.0, 110.0),
        "#880000": (300.0, 140.0),
    }
    lines: list[str] = []
    for color, (x_px, y_px) in control_pixels.items():
        x_proj = (x_px - 120.0) * 12000.0
        y_proj = (80.0 - y_px) * 10000.0
        lon_deg, lat_deg = to_lonlat.transform(x_proj, y_proj)
        lines.append(f"{float(lat_deg):.10f}, {float(lon_deg):.10f} {color}")
        image.putpixel((int(x_px), int(y_px)), tuple(int(color[i : i + 2], 16) for i in (1, 3, 5)))

    image_path = tmp_path / "geo.png"
    map_path = tmp_path / "latlonmap.txt"
    image.save(image_path)
    map_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return image_path, map_path, expected


def test_fit_geostationary_image_mapping_recovers_known_affine_projection(tmp_path: Path) -> None:
    image_path, map_path, expected_pixels = _write_control_fixture(tmp_path)
    assert image_path.exists()
    with Image.open(image_path) as image:
        pixel_map, located = fit_pixel_map_from_image(
            image,
            map_path,
            projection=build_projection(longitude_of_projection_origin=0.0, perspective_point_height=35785831.0, sweep_angle_axis="x"),
        )
        assert len(located) == 4
        for point in located:
            assert point.match_count == 1
            exp_x, exp_y = expected_pixels["#{:02x}{:02x}{:02x}".format(*point.rgb)]
            assert np.isclose(point.x_px, exp_x)
            assert np.isclose(point.y_px, exp_y)

        lon_deg, lat_deg = pixel_map.pixel_to_lonlat(180.0, 95.0)
        to_lonlat = Transformer.from_crs(pixel_map.projection, "EPSG:4326", always_xy=True)
        x_proj = (180.0 - 120.0) * 12000.0
        y_proj = (80.0 - 95.0) * 10000.0
        expected_lon, expected_lat = to_lonlat.transform(x_proj, y_proj)
        assert np.isclose(float(np.asarray(lon_deg).item()), float(expected_lon), atol=1e-10)
        assert np.isclose(float(np.asarray(lat_deg).item()), float(expected_lat), atol=1e-10)

        lon_grid, lat_grid = estimate_pixel_to_lonlat(
            image,
            map_path,
            projection=build_projection(longitude_of_projection_origin=0.0, perspective_point_height=35785831.0, sweep_angle_axis="x"),
        )[2:]
        assert lon_grid.shape == (300, 400)
        assert lat_grid.shape == (300, 400)


def test_quadratic_lonlat_fit_recovers_control_points(tmp_path: Path) -> None:
    image = Image.new("RGB", (200, 160), color=(0, 0, 0))
    projection = build_projection(longitude_of_projection_origin=0.0, perspective_point_height=35785831.0, sweep_angle_axis="x")
    to_lonlat = Transformer.from_crs(projection, "EPSG:4326", always_xy=True)
    control_pixels = [
        ("#ff0000", (20.0, 20.0)),
        ("#cc0000", (60.0, 40.0)),
        ("#aa0000", (100.0, 60.0)),
        ("#880000", (140.0, 80.0)),
        ("#660000", (50.0, 120.0)),
        ("#440000", (160.0, 110.0)),
    ]
    lines: list[str] = []
    for color, (x_px, y_px) in control_pixels:
        x_proj = (x_px - 100.0) * 8000.0 + (y_px - 80.0) * 1000.0
        y_proj = (80.0 - y_px) * 9000.0 + (x_px - 100.0) * 800.0
        lon_deg, lat_deg = to_lonlat.transform(x_proj, y_proj)
        lines.append(f"{float(lat_deg):.10f}, {float(lon_deg):.10f} {color}")
        image.putpixel((int(x_px), int(y_px)), tuple(int(color[i : i + 2], 16) for i in (1, 3, 5)))

    map_path = tmp_path / "latlonmap.txt"
    map_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    pixel_map, located = fit_pixel_lonlat_polynomial_from_image(image, map_path)
    assert len(located) == 6
    for point in located:
        lon_fit, lat_fit = pixel_map.pixel_to_lonlat(point.x_px, point.y_px)
        assert np.isclose(float(np.asarray(lon_fit).item()), float(point.lon_deg), atol=0.02)
        assert np.isclose(float(np.asarray(lat_fit).item()), float(point.lat_deg), atol=0.02)


def test_tps_fit_recovers_control_points(tmp_path: Path) -> None:
    image = Image.new("RGB", (320, 240), color=(0, 0, 0))
    projection = build_projection(longitude_of_projection_origin=0.0, perspective_point_height=35785831.0, sweep_angle_axis="x")
    to_lonlat = Transformer.from_crs(projection, "EPSG:4326", always_xy=True)
    control_pixels = [
        ("#ff0000", (40.0, 35.0)),
        ("#cc0000", (90.0, 50.0)),
        ("#aa0000", (150.0, 80.0)),
        ("#880000", (220.0, 120.0)),
        ("#660000", (260.0, 160.0)),
        ("#440000", (120.0, 190.0)),
    ]
    lines: list[str] = []
    for color, (x_px, y_px) in control_pixels:
        x_proj = (x_px - 160.0) * 9000.0
        y_proj = (120.0 - y_px) * 9000.0
        lon_deg, lat_deg = to_lonlat.transform(x_proj, y_proj)
        lines.append(f"{float(lat_deg):.10f}, {float(lon_deg):.10f} {color}")
        image.putpixel((int(x_px), int(y_px)), tuple(int(color[i : i + 2], 16) for i in (1, 3, 5)))

    map_path = tmp_path / "latlonmap.txt"
    map_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    fitted, located = fit_pixel_lonlat_tps_from_image(image, map_path)
    assert len(located) == 6
    for point in located:
        lon_fit, lat_fit = fitted.pixel_to_lonlat(point.x_px, point.y_px)
        assert np.isclose(float(np.asarray(lon_fit).item()), float(point.lon_deg), atol=0.02)
        assert np.isclose(float(np.asarray(lat_fit).item()), float(point.lat_deg), atol=0.02)
