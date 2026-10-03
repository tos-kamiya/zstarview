import threading

import pytest

from zstarview import water_surface_rings as rings
from zstarview.clouddisc.types import DownloadCancelledError


def point(az, radius=1000, category="sea", index=None, altitude=-2):
    return rings.WaterOverlayPoint(
        "water",
        altitude,
        az,
        radius / 1000,
        scan_distance_m=radius,
        water_category=category,
        scan_azimuth_index=index,
    )


def test_connects_all_distance_rows_and_preserves_original_points():
    points = [
        point(az, radius, altitude=radius / 1000 - 5)
        for radius in (125, 1000, 10000, 64000)
        for az in (10, 12, 14)
    ]
    lines = rings.build_water_surface_ring_polylines(points)
    assert len(lines) == 4
    assert sum(len(line.points) for line in lines) == len(points)
    assert all(
        any(vertex is original for original in points)
        for line in lines
        for vertex in line.points
    )


def test_missing_samples_split_water_runs():
    points = [point(az, index=az // 2) for az in (0, 2, 8, 10, 16)]
    lines = rings.build_water_surface_ring_polylines(points)
    assert [[p.az_deg for p in line.points] for line in lines] == [[0, 2], [8, 10]]


def test_does_not_join_different_water_categories():
    lines = rings.build_water_surface_ring_polylines(
        [point(0), point(2), point(4, category="lake"), point(6, category="lake")]
    )
    assert [line.water_category for line in lines] == ["lake", "sea"]


def test_connects_adjacent_samples_across_azimuth_seam():
    lines = rings.build_water_surface_ring_polylines([point(358), point(0), point(2)])
    assert [[p.az_deg for p in line.points] for line in lines] == [[0, 2], [358, 0]]


def test_cancel_before_grouping():
    event = threading.Event()
    event.set()
    with pytest.raises(DownloadCancelledError):
        rings.build_water_surface_ring_polylines([point(0)], abort_event=event)


def test_export_retains_ring_lines_when_sparse_dots_miss_water(monkeypatch):
    from types import SimpleNamespace

    from zstarview.cli import export_image_layers as layers

    marker = rings.WaterOverlayPolyline("ring", "lake", ())
    captured = {}
    monkeypatch.setattr(
        layers,
        "host",
        lambda: SimpleNamespace(
            _fetch_water_overlay_dots_layer=lambda **kwargs: None,
        ),
    )
    monkeypatch.setattr(
        layers, "_load_or_fetch_water_overlay_footprints", lambda **kwargs: ()
    )
    monkeypatch.setattr(
        layers, "build_water_overlay_polylines", lambda *args, **kwargs: ()
    )
    monkeypatch.setattr(layers, "load_coastline_overlay_polylines", lambda **kwargs: ())

    def generate(*args, **kwargs):
        captured.update(kwargs)
        return (marker,)

    monkeypatch.setattr(layers, "build_water_surface_ring_polylines", generate)

    def sampler(lat, lon):
        return 150.0

    result = layers._fetch_water_overlay_layer(
        viewer_data=SimpleNamespace(
            lat_deg=0,
            lon_deg=0,
            observer_height_m=2,
            ground_elevation_m=100,
            view_center=(0, 0),
            content_fov_deg=90,
        ),
        deadline=None,
        target_ground_sampler=sampler,
    )
    assert result == {"dots": None, "polylines": [marker]}


def test_sea_sampler_keeps_each_radius_for_circumferential_lines(monkeypatch):
    import numpy as np

    from zstarview import water_mask_interface as masks

    distances = np.array([125.0, 500.0, 1000.0])
    monkeypatch.setattr(
        masks, "build_geometric_distance_samples", lambda *args: distances
    )
    monkeypatch.setattr(
        masks,
        "_sample_water_mask_for_lonlat_points_with_stats",
        lambda points, **kwargs: ([True] * len(points), 0),
    )
    points, _ = masks._sample_water_surface_interface_ray_points_for_root_with_stats(
        center_lat_deg=34.6825,
        center_lon_deg=135.1866666666,
        observer_height_m=113.6,
        radius_km=1.0,
        tile_root=masks.DEFAULT_WATER_TILES_ROOT_125M,
        azimuth_step_deg=90.0,
    )
    assert len(points) == 12
    for point in points:
        assert point.scan_distance_m == distances[point.scan_distance_index]
    lines = rings.build_water_surface_ring_polylines(points, azimuth_step_deg=90)
    assert {line.points[0].scan_distance_m for line in lines} == set(distances)
    for line in lines:
        assert len({point.scan_distance_index for point in line.points}) == 1
        assert (
            max(point.alt_deg for point in line.points)
            - min(point.alt_deg for point in line.points)
            < 0.01
        )


def test_sea_resolution_bands_do_not_repeat_nearby_circles(monkeypatch):
    from zstarview import water_mask_interface as masks

    monkeypatch.setattr(
        masks._ZipWaterMaskRoot, "from_cache", staticmethod(lambda: None)
    )
    monkeypatch.setattr(
        masks,
        "_sample_water_mask_for_lonlat_points_with_stats",
        lambda points, **kwargs: ([True] * len(points), 0),
    )
    points, _ = masks.sample_water_surface_interface_points_with_stats(
        observer_lat_deg=34.6825,
        observer_lon_deg=135.1866666666,
        observer_height_m=113.6,
        max_distance_km=32.0,
        azimuth_step_deg=90,
    )
    samples = [(p.scan_distance_m, p.scan_azimuth_index) for p in points]
    assert len(samples) == len(set(samples))
    assert {p.water_category for p in points} == {"sea-125", "sea-250", "sea-500"}
    for p in points:
        if p.water_category == "sea-125":
            assert p.scan_distance_m <= 2000
        elif p.water_category == "sea-250":
            assert 2000 < p.scan_distance_m <= 6000
        else:
            assert p.scan_distance_m > 6000
