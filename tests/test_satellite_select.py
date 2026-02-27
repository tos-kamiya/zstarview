from __future__ import annotations

from zstarview.clouddisc.providers.select import (
    is_satellite_visible,
    pick_satellite,
    visible_satellites,
)


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


def test_meteosat_visibility_available_in_experimental_mode() -> None:
    # London should see METEOSAT geometry when experimental satellites are included.
    assert is_satellite_visible(
        51.5072,
        -0.1276,
        "METEOSAT",
        include_experimental=True,
    )
    sats = visible_satellites(
        51.5072,
        -0.1276,
        ("G16", "G18", "METEOSAT"),
        include_experimental=True,
    )
    assert "METEOSAT" in sats


def test_pick_satellite_auto_can_include_meteosat_when_requested() -> None:
    sat = pick_satellite(
        51.5072,
        -0.1276,
        ("AUTO",),
        include_experimental=True,
    )
    assert sat == "METEOSAT"


def test_pick_satellite_can_select_meteosat_in_manual_mode() -> None:
    sat = pick_satellite(
        51.5072,
        -0.1276,
        ("METEOSAT", "G16"),
        include_experimental=True,
    )
    assert sat == "METEOSAT"
