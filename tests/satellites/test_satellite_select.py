from __future__ import annotations

from zstarview.clouddisc.providers.select import (
    MAX_VISIBLE_CENTRAL_ANGLE_DEG,
    SAT_LON,
    central_angle_deg,
    visible_satellites,
)


def test_goes_visibility_for_europe_examples() -> None:
    # Madrid
    assert central_angle_deg(40.4165, -3.70256, SAT_LON["G19"]) <= MAX_VISIBLE_CENTRAL_ANGLE_DEG
    assert central_angle_deg(40.4165, -3.70256, SAT_LON["G18"]) > MAX_VISIBLE_CENTRAL_ANGLE_DEG

    # London
    assert central_angle_deg(51.5072, -0.1276, SAT_LON["G19"]) <= MAX_VISIBLE_CENTRAL_ANGLE_DEG
    assert central_angle_deg(51.5072, -0.1276, SAT_LON["G18"]) > MAX_VISIBLE_CENTRAL_ANGLE_DEG


def test_visible_satellites_returns_angle_sorted_subset() -> None:
    sats = visible_satellites(35.6764, 139.65, ("G19", "G18", "HIMAWARI"))
    # For Tokyo, HIMAWARI is expected to be the closest visible satellite.
    assert sats
    assert sats[0] == "HIMAWARI"


def test_legacy_g16_alias_normalizes_to_g19() -> None:
    sats = visible_satellites(40.4165, -3.70256, ("G16", "G19", "G18"))

    assert sats[0] == "G19"
