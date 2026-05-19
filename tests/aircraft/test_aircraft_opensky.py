from zstarview.aircraft.opensky import (
    build_observer_bbox,
    normalize_opensky_state_vectors,
)


def test_build_observer_bbox_uses_one_degree_each_side() -> None:
    bbox = build_observer_bbox(35.47, 133.05)
    assert bbox.min_lat == 34.47
    assert bbox.max_lat == 36.47
    assert bbox.min_lon == 132.05
    assert bbox.max_lon == 134.05
    assert bbox.area_sq_deg == 4.0


def test_build_observer_bbox_clamps_world_edges() -> None:
    bbox = build_observer_bbox(89.5, 179.7)
    assert bbox.min_lat == 88.5
    assert bbox.max_lat == 90.0
    assert bbox.min_lon < 178.7
    assert bbox.max_lon == 180.0
    assert bbox.area_sq_deg <= 25.0


def test_build_observer_bbox_expands_longitude_at_high_latitudes() -> None:
    bbox = build_observer_bbox(70.0, 10.0)
    assert bbox.max_lon - bbox.min_lon > 2.0
    assert bbox.area_sq_deg <= 25.0


def test_normalize_opensky_state_vectors_filters_ground_and_missing_positions() -> None:
    payload = {
        "states": [
            ["abcd01", "JAL123 ", "JP", None, 1710000000, 133.1, 35.5, 10000.0, False, 250.0, 90.0, 0.0],
            ["abcd02", " ", "JP", None, 1710000001, 133.2, None, 9000.0, False, 250.0, 90.0, 0.0],
            ["abcd03", "ANA456", "JP", None, 1710000002, 133.3, 35.7, 0.0, True, 10.0, 90.0, 0.0],
            ["abcd04", "ANA789", "JP", None, 1710000003, 133.4, 35.8, 2000.0, False, 20.0, 90.0, 0.0],
        ]
    }

    snapshots = normalize_opensky_state_vectors(payload)

    assert len(snapshots) == 1
    assert snapshots[0].icao24 == "abcd01"
    assert snapshots[0].callsign == "JAL123"
    assert snapshots[0].latitude == 35.5
    assert snapshots[0].longitude == 133.1
