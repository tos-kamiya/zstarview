import polars as pl

from zstarview.location_resolver import PlaceSearchCandidate
from zstarview.gui.famous_star_shortcuts import (
    DEC_BAND_EQUATOR,
    DEC_BAND_NORTH,
    DEC_BAND_SOUTH,
    build_place_search_jump_targets,
    build_search_jump_targets,
    build_named_star_shortcuts,
    classify_declination_band,
    flatten_named_star_shortcuts,
)
from zstarview.gui.window_inputs import prepare_window_catalogs


def test_classify_declination_band_thresholds() -> None:
    assert classify_declination_band(20.0) == DEC_BAND_NORTH
    assert classify_declination_band(19.99) == DEC_BAND_EQUATOR
    assert classify_declination_band(-20.0) == DEC_BAND_SOUTH
    assert classify_declination_band(-19.99) == DEC_BAND_EQUATOR


def test_build_named_star_shortcuts_filters_groups_and_dedupes() -> None:
    df = pl.DataFrame(
        {
            "Name": ["", "Sirius", "Sirius", "Vega", "Canopus", "Altair"],
            "RAh": [1.0, 6.75, 6.76, 18.6, 6.4, 19.8],
            "Dec": [10.0, -16.7, -16.8, 38.7, -52.7, 8.9],
            "Vmag": [1.0, -1.44, -1.30, 0.03, -0.62, 2.50],
        }
    )

    grouped = build_named_star_shortcuts(df, max_vmag=2.0)

    assert len(grouped[DEC_BAND_NORTH]) == 1
    assert grouped[DEC_BAND_NORTH][0].name == "Vega"

    assert len(grouped[DEC_BAND_EQUATOR]) == 1
    assert grouped[DEC_BAND_EQUATOR][0].name == "Sirius"
    assert grouped[DEC_BAND_EQUATOR][0].vmag == -1.44

    assert len(grouped[DEC_BAND_SOUTH]) == 1
    assert grouped[DEC_BAND_SOUTH][0].name == "Canopus"


def test_build_named_star_shortcuts_includes_satellite_entries() -> None:
    df = pl.DataFrame(
        {
            "Name": ["Sirius"],
            "RAh": [6.75],
            "Dec": [-16.7],
            "Vmag": [-1.44],
        }
    )

    grouped = build_named_star_shortcuts(df, max_vmag=2.0, include_satellites=True)

    names = [s.name for s in grouped[DEC_BAND_EQUATOR]]
    assert {"ISS", "JWST", "Voyager 1", "Voyager 2", "Parker"}.issubset(set(names))


def test_flatten_named_star_shortcuts_sorts_globally() -> None:
    df = pl.DataFrame(
        {
            "Name": ["Canopus", "Vega", "Sirius"],
            "RAh": [6.4, 18.6, 6.75],
            "Dec": [-52.7, 38.7, -16.7],
            "Vmag": [-0.62, 0.03, -1.44],
        }
    )
    grouped = build_named_star_shortcuts(df, max_vmag=2.0)
    flat = flatten_named_star_shortcuts(grouped)
    assert [s.name for s in flat] == ["Sirius", "Canopus", "Vega"]


def test_build_named_star_shortcuts_without_vmag_limit_keeps_faint_named() -> None:
    df = pl.DataFrame(
        {
            "Name": ["BrightStar", "FaintStar"],
            "RAh": [1.0, 2.0],
            "Dec": [30.0, -30.0],
            "Vmag": [1.0, 7.9],
        }
    )
    grouped = build_named_star_shortcuts(df, max_vmag=None)
    flat = flatten_named_star_shortcuts(grouped)
    assert [s.name for s in flat] == ["BrightStar", "FaintStar"]


def test_build_search_jump_targets_includes_asterisms() -> None:
    df = pl.DataFrame(
        {
            "Name": ["Theta", "Seven", "Gamma", "Kappa", "Lambda", "TX", "Iota", "Sirius"],
            "RAh": [23.46, 23.80, 23.29, 23.45, 23.70, 23.77, 23.66, 6.75],
            "Dec": [6.38, -2.76, 3.28, 1.25, 1.78, 3.49, 5.63, -16.7],
            "Vmag": [4.27, 5.49, 3.70, 4.95, 4.49, 4.95, 4.13, -1.44],
            "SourceId": ["HIP115830", "HIP117375", "HIP114971", "HIP115738", "HIP116928", "HIP117245", "HIP116771", "HIP32349"],
        }
    )

    targets = build_search_jump_targets(df, include_satellites=True)

    assert any(t.label == "Sirius" and t.kind == "star" for t in targets)
    assert any(t.label == "ISS" and t.kind == "satellite" for t in targets)
    circlet = next(t for t in targets if t.label == "Circlet of Pisces")
    assert circlet.kind == "asterism"
    assert circlet.subtitle == "Asterism"


def test_prepare_window_catalogs_search_targets_exclude_satellite_shortcuts() -> None:
    df = pl.DataFrame(
        {
            "Name": ["Sirius"],
            "RAh": [6.75],
            "Dec": [-16.7],
            "Vmag": [-1.44],
            "B-V": [0.0],
            "SourceId": ["HIP32349"],
        }
    )

    catalogs = prepare_window_catalogs(df)

    assert all(target.kind != "satellite" for target in catalogs.named_stars_search_all)


def test_build_place_search_jump_targets_uses_importance_and_coordinates() -> None:
    targets = build_place_search_jump_targets(
        [
            PlaceSearchCandidate(
                name="Tokyo Station",
                display_name="Tokyo Station, Chiyoda, Tokyo, Japan",
                latitude_deg=35.681236,
                longitude_deg=139.767125,
                category="railway",
                type_name="station",
                importance=0.9,
            ),
            PlaceSearchCandidate(
                name="Tokyo",
                display_name="Tokyo, Japan",
                latitude_deg=35.67686,
                longitude_deg=139.76389,
                category="boundary",
                type_name="administrative",
                importance=0.7,
            ),
        ]
    )

    assert [target.label for target in targets] == ["Tokyo Station", "Tokyo"]
    assert targets[0].kind == "place"
    assert targets[0].subtitle == "Place / railway / station"
    assert targets[0].latitude_deg == 35.681236
    assert targets[0].longitude_deg == 139.767125
    assert targets[0].object_key == "Tokyo Station, Chiyoda, Tokyo, Japan"
