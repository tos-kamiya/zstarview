from __future__ import annotations

import pytest

from zstarview.location_resolver.place_projection import (
    PlaceTargetProjection,
    project_place_target_to_altaz,
    project_place_targets_to_altaz,
)


def test_project_place_target_to_altaz_returns_overhead_for_same_latlon_with_height() -> None:
    projection = project_place_target_to_altaz(
        observer_latitude_deg=35.0,
        observer_longitude_deg=139.0,
        observer_height_m=0.0,
        target_latitude_deg=35.0,
        target_longitude_deg=139.0,
        target_height_m=1000.0,
    )

    assert isinstance(projection, PlaceTargetProjection)
    assert projection.alt_deg == pytest.approx(90.0, abs=1e-9)
    assert projection.az_deg == pytest.approx(0.0, abs=1e-9)
    assert projection.distance_km == pytest.approx(1.0, abs=1e-9)
    assert projection.target_latitude_deg == 35.0
    assert projection.target_longitude_deg == 139.0
    assert projection.target_height_m == 1000.0


def test_project_place_target_to_altaz_points_east_for_equatorial_target() -> None:
    projection = project_place_target_to_altaz(
        observer_latitude_deg=0.0,
        observer_longitude_deg=0.0,
        observer_height_m=0.0,
        target_latitude_deg=0.0,
        target_longitude_deg=1.0,
        target_height_m=0.0,
    )

    assert projection.az_deg == pytest.approx(90.0, abs=1e-6)
    assert projection.alt_deg == pytest.approx(-0.5, abs=0.02)
    assert projection.distance_km == pytest.approx(111.3, abs=0.5)


def test_project_place_target_to_altaz_points_north_for_meridional_target() -> None:
    projection = project_place_target_to_altaz(
        observer_latitude_deg=35.0,
        observer_longitude_deg=139.0,
        observer_height_m=0.0,
        target_latitude_deg=36.0,
        target_longitude_deg=139.0,
        target_height_m=0.0,
    )

    assert projection.az_deg == pytest.approx(0.0, abs=1e-6)
    assert projection.alt_deg < 0.0
    assert projection.distance_km == pytest.approx(110.9, abs=1.0)


def test_project_place_targets_to_altaz_matches_scalar_projection() -> None:
    batch = project_place_targets_to_altaz(
        observer_latitude_deg=35.0,
        observer_longitude_deg=139.0,
        observer_height_m=0.0,
        target_latitude_deg=[35.0, 35.0, 36.0],
        target_longitude_deg=[139.0, 140.0, 139.0],
        target_height_m=[1000.0, 0.0, 0.0],
    )

    assert len(batch) == 3
    assert batch[0] == project_place_target_to_altaz(
        observer_latitude_deg=35.0,
        observer_longitude_deg=139.0,
        observer_height_m=0.0,
        target_latitude_deg=35.0,
        target_longitude_deg=139.0,
        target_height_m=1000.0,
    )
    assert batch[1] == project_place_target_to_altaz(
        observer_latitude_deg=35.0,
        observer_longitude_deg=139.0,
        observer_height_m=0.0,
        target_latitude_deg=35.0,
        target_longitude_deg=140.0,
        target_height_m=0.0,
    )
    assert batch[2] == project_place_target_to_altaz(
        observer_latitude_deg=35.0,
        observer_longitude_deg=139.0,
        observer_height_m=0.0,
        target_latitude_deg=36.0,
        target_longitude_deg=139.0,
        target_height_m=0.0,
    )
