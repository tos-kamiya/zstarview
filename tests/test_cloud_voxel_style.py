from __future__ import annotations

import threading

import numpy as np
import pytest

from zstarview import cloud_voxel_style, night_lights
from zstarview.clouddisc.workers import cloud_voxel_worker
from zstarview.render.ground_mask import inverse_project_disc
from zstarview.types import ScreenGeometry


def test_render_cloud_voxel_image_shares_appearance_and_marks(monkeypatch) -> None:
    source = object()
    raw_image = np.zeros((3, 3, 4), dtype=np.uint8)
    styled_image = raw_image.copy()
    render_call: dict[str, object] = {}
    mark_call: dict[str, object] = {}

    def render(source_arg, **kwargs):
        render_call["source"] = source_arg
        render_call.update(kwargs)
        return raw_image

    def marks(image, style, **kwargs):
        mark_call["image"] = image
        mark_call["style"] = style
        mark_call.update(kwargs)
        return styled_image

    monkeypatch.setattr(cloud_voxel_worker, "render_cloud_voxels_in_subprocess", render)
    monkeypatch.setattr(night_lights, "night_light_strength_factor", lambda _alt: 0.25)
    monkeypatch.setattr(cloud_voxel_style, "apply_cloud_voxel_marks", marks)
    style = cloud_voxel_style.CloudVoxelMarkStyle("cutout", 8)

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
        mark_style=style,
        moon_altaz_phase=(45.0, 120.0, 180.0),
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
        "moon_altaz_phase": (45.0, 120.0, 180.0),
        "native_work_lock": None,
    }
    assert mark_call == {
        "image": raw_image,
        "style": style,
        "view_center": (45.0, 180.0),
        "content_fov_deg": 110.0,
    }


def test_parent_style_native_work_is_locked_but_subprocess_wait_is_not(monkeypatch) -> None:
    lock = threading.Lock()
    raw = np.zeros((3, 3, 4), dtype=np.uint8)

    def render(source, **kwargs):
        assert not lock.locked()
        assert kwargs["native_work_lock"] is lock
        return raw

    def night_mix(alt):
        assert lock.locked()
        return 1.0

    def marks(image, style, **kwargs):
        assert lock.locked()
        return image

    monkeypatch.setattr(cloud_voxel_worker, "render_cloud_voxels_in_subprocess", render)
    monkeypatch.setattr(night_lights, "night_light_strength_factor", night_mix)
    monkeypatch.setattr(cloud_voxel_style, "apply_cloud_voxel_marks", marks)
    result = cloud_voxel_style.render_cloud_voxel_image(
        object(), lat=35.0, lon=135.0, view_center=(45.0, 180.0),
        edge_fov_deg=90.0, content_fov_deg=110.0, radius_px=1,
        sun_alt_deg=-20.0, sun_az_deg=220.0, request_id=1,
        native_work_lock=lock,
    )
    assert result is raw
    assert not lock.locked()


def test_cloud_hole_points_default_rings() -> None:
    alt, az = cloud_voxel_style.cloud_hole_altaz_points(126)

    rings = sorted({round(float(value), 6) for value in alt}, reverse=True)
    assert rings == [90.0, 78.0, 66.0, 54.0, 42.0, 30.0, 18.0, 6.0]
    counts = [int(np.count_nonzero(np.isclose(alt, ring))) for ring in rings]
    assert counts == [1, 5, 11, 15, 20, 23, 25, 26]
    first_ring = np.sort(az[np.isclose(alt, 78.0)])
    assert np.allclose(first_ring, [0.0, 72.0, 144.0, 216.0, 288.0])
    golden = 180.0 * (3.0 - np.sqrt(5.0))
    for ring_index, ring in enumerate(rings[1:]):
        ring_az = np.sort(az[np.isclose(alt, ring)])
        step = 360.0 / len(ring_az)
        assert ring_az[0] == pytest.approx((ring_index * golden) % step)


@pytest.mark.parametrize("count", [1, 30, 126, 250])
def test_cloud_hole_points_form_even_concentric_rings(count) -> None:
    alt, az = cloud_voxel_style.cloud_hole_altaz_points(count)

    assert alt.shape == az.shape == (count,)
    assert alt[0] == 90.0
    assert np.all((alt > 0.0) & (alt <= 90.0))
    if count == 1:
        return
    rings = np.unique(alt)[::-1]
    step = 90.0 - rings[1]
    assert np.allclose(-np.diff(rings), step)
    assert rings[-1] == pytest.approx(step / 2.0)
    for ring in rings[1:]:
        ring_az = np.sort(az[alt == ring])
        spacing = 360.0 * np.cos(np.radians(ring)) / len(ring_az)
        assert np.allclose(np.diff(ring_az), 360.0 / len(ring_az))
        # Along-ring spacing stays close to the ring step.
        assert 0.75 * step < spacing < 1.4 * step


@pytest.mark.parametrize(
    ("view_center", "content_fov_deg"),
    [((45.0, 180.0), 110.0), ((20.0, 73.0), 40.0), ((80.0, 300.0), 120.0)],
)
def test_project_altaz_to_cloud_image_inverts_worker_projection(
    view_center, content_fov_deg
) -> None:
    radius = 64
    size = radius * 2 + 1
    alt, az, inside = inverse_project_disc(
        size,
        size,
        ScreenGeometry((radius, radius), radius),
        view_center,
        edge_fov_deg=content_fov_deg,
        content_fov_deg=content_fov_deg,
    )
    ys, xs = np.nonzero(inside)

    x, y, _ = cloud_voxel_style.project_altaz_to_cloud_image(
        alt,
        az,
        width=size,
        height=size,
        view_center=view_center,
        content_fov_deg=content_fov_deg,
    )

    assert np.max(np.abs(x - xs)) < 0.01
    assert np.max(np.abs(y - ys)) < 0.01


def test_apply_cloud_voxel_holes_clears_alpha_at_fixed_altaz_points(monkeypatch) -> None:
    monkeypatch.setattr(cloud_voxel_style, "CLOUD_VOXEL_HOLE_COUNT", 46)
    size = 513
    image = np.full((size, size, 4), 200, dtype=np.uint8)
    view_center = (90.0 - 1.0e-3, 0.0)

    result = cloud_voxel_style.apply_cloud_voxel_holes(
        image, radius_px=4, view_center=view_center, content_fov_deg=100.0
    )

    alt, az = cloud_voxel_style.cloud_hole_altaz_points(46)
    xs, ys, _ = cloud_voxel_style.project_altaz_to_cloud_image(
        alt, az, width=size, height=size, view_center=view_center, content_fov_deg=100.0
    )
    above = alt > 0.0
    for x, y in zip(xs[above], ys[above]):
        assert result[int(round(y)), int(round(x)), 3] == 0
    assert np.array_equal(result[..., :3], image[..., :3])
    count = int(np.count_nonzero(above))
    cleared = int(np.count_nonzero(result[..., 3] < 200))
    # Holes of radius 4 px, plus a one-pixel anti-aliased rim.
    assert count * np.pi * 3.5**2 < cleared < count * np.pi * 5.0**2
    assert image[0, 0, 3] == 200  # input is not modified


def test_apply_cloud_voxel_holes_follows_view_direction() -> None:
    size = 257
    image = np.full((size, size, 4), 255, dtype=np.uint8)
    alt, az = cloud_voxel_style.cloud_hole_altaz_points(
        cloud_voxel_style.CLOUD_VOXEL_HOLE_COUNT
    )
    candidates = np.flatnonzero((alt > 20.0) & (alt < 60.0))
    target = int(candidates[0])
    view_center = (float(alt[target]), float(az[target]))

    result = cloud_voxel_style.apply_cloud_voxel_holes(
        image, radius_px=4, view_center=view_center, content_fov_deg=60.0
    )

    assert result[128, 128, 3] == 0


def test_apply_cloud_voxel_marks_zero_value_disables_marks() -> None:
    image = np.full((65, 65, 4), 255, dtype=np.uint8)
    for mode in cloud_voxel_style.CLOUD_VOXEL_MARK_MODES:
        result = cloud_voxel_style.apply_cloud_voxel_marks(
            image,
            cloud_voxel_style.CloudVoxelMarkStyle(mode, 0),
            view_center=(90.0, 0.0),
            content_fov_deg=100.0,
        )
        assert np.array_equal(result, image)


def test_apply_cloud_voxel_marks_dispatches_cutout() -> None:
    image = np.full((65, 65, 4), 255, dtype=np.uint8)
    result = cloud_voxel_style.apply_cloud_voxel_marks(
        image,
        cloud_voxel_style.CloudVoxelMarkStyle("cutout", 4),
        view_center=(90.0, 0.0),
        content_fov_deg=100.0,
    )
    assert np.array_equal(
        result, cloud_voxel_style.apply_cloud_voxel_cutout(image, count=4)
    )
