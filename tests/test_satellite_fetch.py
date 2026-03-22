from __future__ import annotations

from skyfield.api import EarthSatellite

from zstarview.satellites.fetch import (
    CELESTRAK_GROUP_BY_KEY,
    build_celestrak_group_url,
    build_earth_satellites,
    filter_records_for_group,
    normalize_celestrak_omm_payload,
)


def _sample_record() -> dict[str, object]:
    return {
        "OBJECT_NAME": "ISS (ZARYA)",
        "OBJECT_ID": "1998-067A",
        "EPOCH": "2026-03-22T12:00:00.000000",
        "MEAN_MOTION": "15.50000000",
        "ECCENTRICITY": "0.0006703",
        "INCLINATION": "51.6400",
        "RA_OF_ASC_NODE": "257.7243",
        "ARG_OF_PERICENTER": "130.5360",
        "MEAN_ANOMALY": "325.0288",
        "EPHEMERIS_TYPE": "0",
        "CLASSIFICATION_TYPE": "U",
        "NORAD_CAT_ID": "25544",
        "ELEMENT_SET_NO": "999",
        "REV_AT_EPOCH": "12345",
        "BSTAR": "0.00010234",
        "MEAN_MOTION_DOT": "0.00002182",
        "MEAN_MOTION_DDOT": "0.00000000",
    }


def test_celestrak_group_mapping_covers_initial_layers() -> None:
    assert CELESTRAK_GROUP_BY_KEY == {
        "station": "stations",
        "starlink": "starlink",
    }


def test_build_celestrak_group_url_uses_group_and_json_format() -> None:
    url = build_celestrak_group_url("stations")

    assert "GROUP=stations" in url
    assert "FORMAT=json" in url


def test_normalize_celestrak_omm_payload_filters_non_dict_rows() -> None:
    records = normalize_celestrak_omm_payload([_sample_record(), "bad", 123, {"OBJECT_NAME": "X"}])

    assert len(records) == 2
    assert records[0]["OBJECT_NAME"] == "ISS (ZARYA)"
    assert records[1]["OBJECT_NAME"] == "X"


def test_build_earth_satellites_from_omm_records() -> None:
    satellites = build_earth_satellites([_sample_record()])

    assert len(satellites) == 1
    assert isinstance(satellites[0], EarthSatellite)
    assert satellites[0].name == "ISS (ZARYA)"


def test_filter_records_for_station_keeps_only_supported_station_targets() -> None:
    css_tianhe = dict(_sample_record())
    css_tianhe["OBJECT_NAME"] = "CSS (TIANHE)"
    css_tianhe["NORAD_CAT_ID"] = "48274"
    crew_dragon = dict(_sample_record())
    crew_dragon["OBJECT_NAME"] = "CREW DRAGON 12"
    crew_dragon["NORAD_CAT_ID"] = "99999"

    filtered = filter_records_for_group("station", [_sample_record(), css_tianhe, crew_dragon])

    assert [record["OBJECT_NAME"] for record in filtered] == ["ISS (ZARYA)", "CSS (TIANHE)"]
