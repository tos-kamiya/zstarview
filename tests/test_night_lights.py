from __future__ import annotations

import numpy as np
from PySide6.QtCore import QRectF
from PySide6.QtGui import QImage, QPainter
from PySide6.QtWidgets import QApplication

from zstarview import night_lights
from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.render.night_lights import draw_night_light_glow
from zstarview.types import ScreenGeometry

app = QApplication.instance() or QApplication([])


def test_read_or_build_manifest_parses_tile_urls(tmp_path, monkeypatch) -> None:
    html = "\n".join(
        f'<a href="https://assets.science.nasa.gov/content/dam/science/esd/eo/images/'
        f'imagerecords/144000/144897/BlackMarble_2016_{tile}_geo_gray.tif">{tile}</a>'
        for tile in night_lights.NIGHT_LIGHTS_TILE_NAMES
    )
    monkeypatch.setattr(night_lights, "_read_url", lambda *_args, **_kwargs: html)

    manifest = night_lights._read_or_build_manifest(cache_root=tmp_path)

    assert manifest["dataset_version"] == night_lights.NIGHT_LIGHTS_DATASET_VERSION
    assert manifest["source_page_url"] == night_lights.NIGHT_LIGHTS_PAGE_URL
    assert set(manifest["tile_urls"]) == set(night_lights.NIGHT_LIGHTS_TILE_NAMES)
    assert (tmp_path / night_lights.NIGHT_LIGHTS_DATASET_VERSION / "manifest.json").exists()


def test_compute_night_light_glow_profile_returns_none_for_daytime() -> None:
    profile = night_lights.compute_night_light_glow_profile(
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        sun_alt_deg=1.0,
    )
    assert profile is None


def test_compute_night_light_glow_profile_uses_mocked_sampling(tmp_path, monkeypatch) -> None:
    tile_paths = {}
    for tile in night_lights.NIGHT_LIGHTS_TILE_NAMES:
        path = tmp_path / f"{tile}.tif"
        path.write_text("dummy", encoding="utf-8")
        tile_paths[tile] = path

    monkeypatch.setattr(night_lights, "_ensure_night_light_tiles", lambda **_kwargs: tile_paths)
    monkeypatch.setattr(
        night_lights,
        "_build_azimuth_grid",
        lambda terrain_profile_altaz=None: (
            np.asarray([0.0, 90.0, 180.0], dtype=np.float64),
            np.zeros(3, dtype=np.float64),
        ),
    )
    monkeypatch.setattr(
        night_lights,
        "_sample_ray_brightness",
        lambda **kwargs: float(kwargs["azimuth_deg"]) / 180.0,
    )

    profile = night_lights.compute_night_light_glow_profile(
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        sun_alt_deg=-5.0,
    )

    assert profile is not None
    assert len(profile.samples) == 3
    strengths = [sample.strength for sample in profile.samples]
    assert strengths[0] <= strengths[1] <= strengths[2]


def test_draw_night_light_glow_smoke() -> None:
    image = QImage(200, 100, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        profile = night_lights.NightLightGlowProfile(
            samples=(
                night_lights.NightLightGlowSample(azimuth_deg=170.0, horizon_alt_deg=0.0, strength=0.4),
                night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=0.0, strength=1.0),
                night_lights.NightLightGlowSample(azimuth_deg=190.0, horizon_alt_deg=0.0, strength=0.4),
            ),
            sun_alt_deg=-5.0,
        )
        draw_night_light_glow(
            painter,
            geometry=ScreenGeometry(center=(100, 50), radius=45),
            viewport_rect=QRectF(0.0, 0.0, 200.0, 100.0),
            profile=profile,
            terrain_profile_altaz=[
                (0.0, 170.0),
                (0.0, 180.0),
                (0.0, 190.0),
            ],
            view_center=(0.0, 180.0),
            theme=THEME_STYLES_BY_PRESET["night"],
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        )
    finally:
        painter.end()

    assert any(
        image.pixelColor(x, y).alpha() > 0
        for x in range(image.width())
        for y in range(image.height())
    )
    assert any(
        (
            (color := image.pixelColor(x, y)).alpha() > 0
            and color.red() > color.green()
            and color.green() > color.blue()
        )
        for x in range(image.width())
        for y in range(image.height())
    )
