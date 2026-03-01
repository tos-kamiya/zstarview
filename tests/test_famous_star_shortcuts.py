import polars as pl

from zstarview.ui.famous_star_shortcuts import (
    DEC_BAND_EQUATOR,
    DEC_BAND_NORTH,
    DEC_BAND_SOUTH,
    build_famous_star_shortcuts,
    classify_declination_band,
)


def test_classify_declination_band_thresholds() -> None:
    assert classify_declination_band(20.0) == DEC_BAND_NORTH
    assert classify_declination_band(19.99) == DEC_BAND_EQUATOR
    assert classify_declination_band(-20.0) == DEC_BAND_SOUTH
    assert classify_declination_band(-19.99) == DEC_BAND_EQUATOR


def test_build_famous_star_shortcuts_filters_groups_and_dedupes() -> None:
    df = pl.DataFrame(
        {
            "Name": ["", "Sirius", "Sirius", "Vega", "Canopus", "Altair"],
            "RAh": [1.0, 6.75, 6.76, 18.6, 6.4, 19.8],
            "Dec": [10.0, -16.7, -16.8, 38.7, -52.7, 8.9],
            "Vmag": [1.0, -1.44, -1.30, 0.03, -0.62, 2.50],
        }
    )

    grouped = build_famous_star_shortcuts(df, max_vmag=2.0)

    assert len(grouped[DEC_BAND_NORTH]) == 1
    assert grouped[DEC_BAND_NORTH][0].name == "Vega"

    assert len(grouped[DEC_BAND_EQUATOR]) == 1
    assert grouped[DEC_BAND_EQUATOR][0].name == "Sirius"
    assert grouped[DEC_BAND_EQUATOR][0].vmag == -1.44

    assert len(grouped[DEC_BAND_SOUTH]) == 1
    assert grouped[DEC_BAND_SOUTH][0].name == "Canopus"
