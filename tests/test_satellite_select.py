from __future__ import annotations

from zstarview.clouddisc.providers.select import is_satellite_visible, visible_satellites


def test_goes_visibility_for_europe_examples() -> None:
    # Madrid
    assert is_satellite_visible(40.4165, -3.70256, "G16")
    assert not is_satellite_visible(40.4165, -3.70256, "G18")

    # London
    assert is_satellite_visible(51.5072, -0.1276, "G16")
    assert not is_satellite_visible(51.5072, -0.1276, "G18")


def test_visible_satellites_returns_angle_sorted_subset() -> None:
    sats = visible_satellites(35.6764, 139.65, ("G16", "G18", "HIMAWARI"))
    # For Tokyo, HIMAWARI is expected to be the closest visible satellite.
    assert sats
    assert sats[0] == "HIMAWARI"
