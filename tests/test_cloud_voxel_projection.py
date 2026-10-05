from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from PySide6.QtGui import QImage, QPainter

from zstarview.clouddisc.workers import cloud_voxel_worker as worker
from zstarview.gui.composite import SkyCompositorCache
from zstarview.render.geometry import get_screen_geometry
from zstarview.render.ground_mask import inverse_project_disc
from zstarview.types import ViewProjection


@pytest.mark.parametrize(
    "width,height,edge_fov,content_fov",
    [
        (660, 741, 95.0, 115.0),
        (320, 900, 95.0, 115.0),
        (1000, 600, 95.0, 115.0),
        (600, 600, 90.0, 90.0),
    ],
)
def test_content_cloud_raster_covers_viewport_and_preserves_projection(
    monkeypatch: pytest.MonkeyPatch,
    width: int,
    height: int,
    edge_fov: float,
    content_fov: float,
) -> None:
    # Encode each viewing ray's altitude in red. This lets us check generation
    # and actual QPainter placement together without satellite data or tracing.
    rendered_angles: list[float] = []

    def shade(source, lat, lon, altitudes, azimuths, sun_alt, sun_az, base, **kwargs):
        rgb = np.zeros_like(base)
        rgb[:, 0] = (altitudes + 90.0) / 180.0
        rgb[:, 2] = 1.0
        cos_angle = np.cos(np.radians(altitudes)) * np.cos(np.radians(azimuths))
        rendered_angles.append(
            float(np.degrees(np.arccos(np.clip(cos_angle, -1, 1))).max())
        )
        return rgb, {}, np.zeros(len(altitudes), dtype=np.float32)

    def run_worker(command, **kwargs):
        result = worker._worker_main(*(Path(value) for value in command[-3:]))
        return SimpleNamespace(returncode=result, stderr="")

    monkeypatch.setattr(worker, "shade_native_voxels", shade)
    monkeypatch.setattr(worker.subprocess, "run", run_worker)
    clouds = worker.render_cloud_voxels_in_subprocess(
        SimpleNamespace(source_key="projection-test"),
        lat=0.0,
        lon=0.0,
        view_center=(0.0, 0.0),
        edge_fov_deg=edge_fov,
        content_fov_deg=content_fov,
        radius_px=256,
        sun_alt_deg=-90.0,
        sun_az_deg=0.0,
        sunlight_mix=0.0,
        night_color_rgb=(0.36, 0.385, 0.44),
        request_id=1,
    )
    assert clouds.shape == (513, 513, 4)
    assert rendered_angles[0] == pytest.approx(content_fov, abs=0.01)
    assert clouds[0, 0, 3] == 0  # Outside the content disc.

    geometry = get_screen_geometry(
        width, height, 0.0, edge_fov_deg=edge_fov, content_fov_deg=content_fov
    )
    image = QImage(width, height, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        SkyCompositorCache().draw_cloud_voxel_overlay(
            painter,
            geometry=geometry,
            projection=ViewProjection(
                view_center=(0.0, 0.0),
                edge_fov_deg=edge_fov,
                content_fov_deg=content_fov,
            ),
            cloud_rgba=clouds,
            cloud_alpha=1.0,
        )
    finally:
        painter.end()

    cx, cy = geometry.center
    for x, y in [
        (cx, 0),
        (cx, 2),
        (cx, cy),
        (cx, cy - geometry.radius // 2),
        (min(width - 1, cx + geometry.radius // 2), cy),
    ]:
        altitudes, _, inside = inverse_project_disc(
            1,
            1,
            geometry,
            (0.0, 0.0),
            edge_fov_deg=edge_fov,
            content_fov_deg=content_fov,
            origin=(x, y),
        )
        color = image.pixelColor(x, y)
        if not inside[0, 0]:
            assert color.alpha() == 0
            continue
        assert color.alpha() == 255
        assert color.redF() == pytest.approx(
            (float(altitudes[0]) + 90.0) / 180.0, abs=0.01
        )
