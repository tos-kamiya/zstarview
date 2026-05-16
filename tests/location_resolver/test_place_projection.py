from __future__ import annotations

import pytest

from zstarview.location_resolver.place_projection import (
    PlaceTargetProjection,
    project_place_targets_to_altaz,
)


def test_project_place_targets_to_altaz_projects_batch_targets() -> None:
    batch = project_place_targets_to_altaz(
        observer_latitude_deg=35.0,
        observer_longitude_deg=139.0,
        observer_height_m=0.0,
        target_latitude_deg=[35.0, 35.0, 36.0],
        target_longitude_deg=[139.0, 140.0, 139.0],
        target_height_m=[1000.0, 0.0, 0.0],
    )

    assert len(batch) == 3
    first, second, third = batch
    assert isinstance(first, PlaceTargetProjection)
    assert first.alt_deg == pytest.approx(90.0, abs=1e-9)
    assert first.az_deg == pytest.approx(0.0, abs=1e-9)
    assert first.distance_km == pytest.approx(1.0, abs=1e-9)
    assert first.target_latitude_deg == 35.0
    assert first.target_longitude_deg == 139.0
    assert first.target_height_m == 1000.0

    assert second.az_deg == pytest.approx(89.7132, abs=1e-4)
    assert second.alt_deg < 0.0
    assert second.distance_km == pytest.approx(91.3, abs=0.5)

    assert third.az_deg == pytest.approx(0.0, abs=1e-6)
    assert third.alt_deg < 0.0
    assert third.distance_km == pytest.approx(110.9, abs=1.0)
