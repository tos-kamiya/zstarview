from __future__ import annotations

import numpy as np

from zstarview import cloud_voxel_style, night_lights
from zstarview.clouddisc.workers import cloud_voxel_worker


def test_render_cloud_voxel_image_shares_appearance_and_cutout(monkeypatch) -> None:
    source = object()
    raw_image = np.zeros((3, 3, 4), dtype=np.uint8)
    styled_image = raw_image.copy()
    render_call: dict[str, object] = {}
    cutout_call: dict[str, object] = {}

    def render(source_arg, **kwargs):
        render_call["source"] = source_arg
        render_call.update(kwargs)
        return raw_image

    def cutout(image, *, count):
        cutout_call["image"] = image
        cutout_call["count"] = count
        return styled_image

    monkeypatch.setattr(cloud_voxel_worker, "render_cloud_voxels_in_subprocess", render)
    monkeypatch.setattr(night_lights, "night_light_strength_factor", lambda _alt: 0.25)
    monkeypatch.setattr(cloud_voxel_style, "apply_cloud_voxel_cutout", cutout)

    result = cloud_voxel_style.render_cloud_voxel_image(
        source,
        lat=35.0,
        lon=135.0,
        view_center=(45.0, 180.0),
        edge_fov_deg=90.0,
        content_fov_deg=110.0,
        radius_px=256,
        sun_alt_deg=-6.5,
        sun_az_deg=220.0,
        request_id=17,
        timeout_s=45.0,
        cutout_count=8,
    )

    assert result is styled_image
    assert render_call == {
        "source": source,
        "lat": 35.0,
        "lon": 135.0,
        "view_center": (45.0, 180.0),
        "edge_fov_deg": 90.0,
        "content_fov_deg": 110.0,
        "radius_px": 256,
        "sun_alt_deg": -6.5,
        "sun_az_deg": 220.0,
        "sunlight_mix": 0.75,
        "night_color_rgb": cloud_voxel_style.CLOUD_VOXEL_NIGHT_COLOR_RGB,
        "request_id": 17,
        "timeout_s": 45.0,
    }
    assert cutout_call == {"image": raw_image, "count": 8}
