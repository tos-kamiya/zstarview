from __future__ import annotations

import logging

from skyfield.api import EarthSatellite

from zstarview.satellites.fetch import (
    CELESTRAK_GROUP_BY_KEY,
    WHERETHEISS_API_URL,
    build_celestrak_group_url,
    build_earth_satellites,
    build_wheretheiss_tle_url,
    extract_record_source,
    fetch_iss_records,
    filter_records_for_group,
    normalize_celestrak_omm_payload,
    normalize_wheretheiss_tle_payload,
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


def _sample_wheretheiss_payload() -> dict[str, object]:
    return {
        "requested_timestamp": 1774176000,
        "tle_timestamp": 1774166400,
        "id": "25544",
        "name": "iss",
        "header": "ISS (ZARYA)",
        "line1": "1 25544U 98067A   26081.50000000  .00010234  00000-0  10234-3 0  9991",
        "line2": "2 25544  51.6400 257.7243 0006703 130.5360 325.0288 15.50000000123456",
    }


def test_celestrak_group_mapping_covers_initial_layers() -> None:
    assert CELESTRAK_GROUP_BY_KEY == {
        "iss": "stations",
    }


def test_build_celestrak_group_url_uses_group_and_json_format() -> None:
    url = build_celestrak_group_url("stations")

    assert "GROUP=stations" in url
    assert "FORMAT=json" in url


def test_build_wheretheiss_tle_url_targets_iss_tles_endpoint() -> None:
    url = build_wheretheiss_tle_url()

    assert url == f"{WHERETHEISS_API_URL}/satellites/25544/tles"


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


def test_normalize_wheretheiss_tle_payload_builds_tle_backed_record() -> None:
    records = normalize_wheretheiss_tle_payload(_sample_wheretheiss_payload())

    assert len(records) == 1
    assert records[0]["OBJECT_NAME"] == "ISS (ZARYA)"
    assert extract_record_source(records) == "wheretheiss"


def test_build_earth_satellites_from_wheretheiss_tle_records() -> None:
    satellites = build_earth_satellites(normalize_wheretheiss_tle_payload(_sample_wheretheiss_payload()))

    assert len(satellites) == 1
    assert isinstance(satellites[0], EarthSatellite)
    assert satellites[0].name == "ISS (ZARYA)"


def test_filter_records_for_iss_keeps_only_iss() -> None:
    css_tianhe = dict(_sample_record())
    css_tianhe["OBJECT_NAME"] = "CSS (TIANHE)"
    css_tianhe["NORAD_CAT_ID"] = "48274"
    crew_dragon = dict(_sample_record())
    crew_dragon["OBJECT_NAME"] = "CREW DRAGON 12"
    crew_dragon["NORAD_CAT_ID"] = "99999"

    filtered = filter_records_for_group("iss", [_sample_record(), css_tianhe, crew_dragon])

    assert [record["OBJECT_NAME"] for record in filtered] == ["ISS (ZARYA)"]


def test_fetch_iss_records_uses_celestrak_as_fallback(monkeypatch) -> None:
    monkeypatch.setattr(
        "zstarview.satellites.fetch.fetch_wheretheiss_iss_tle",
        lambda **kwargs: (_ for _ in ()).throw(RuntimeError("primary down")),
    )
    monkeypatch.setattr(
        "zstarview.satellites.fetch.fetch_celestrak_group_by_key",
        lambda group_key, **kwargs: [dict(_sample_record(), _SOURCE="celestrak-fallback")] if group_key == "iss" else [],
    )

    records = fetch_iss_records("iss")

    assert len(records) == 1
    assert extract_record_source(records) == "celestrak-fallback"


def test_fetch_iss_records_logs_primary_timeout_before_fallback(caplog, monkeypatch) -> None:
    monkeypatch.setattr(
        "zstarview.satellites.fetch.fetch_wheretheiss_iss_tle",
        lambda **kwargs: (_ for _ in ()).throw(RuntimeError("timed out")),
    )
    monkeypatch.setattr(
        "zstarview.satellites.fetch.fetch_celestrak_group_by_key",
        lambda group_key, **kwargs: [dict(_sample_record(), _SOURCE="celestrak-fallback")] if group_key == "iss" else [],
    )

    with caplog.at_level(logging.WARNING):
        records = fetch_iss_records("iss")

    assert len(records) == 1
    assert "ISS fetch failed via wheretheiss.at: timed out" in caplog.text


def test_fetch_iss_records_logs_fallback_failure_too(caplog, monkeypatch) -> None:
    monkeypatch.setattr(
        "zstarview.satellites.fetch.fetch_wheretheiss_iss_tle",
        lambda **kwargs: (_ for _ in ()).throw(RuntimeError("primary unavailable")),
    )
    monkeypatch.setattr(
        "zstarview.satellites.fetch.fetch_celestrak_group_by_key",
        lambda group_key, **kwargs: (_ for _ in ()).throw(RuntimeError("fallback timed out")),
    )

    with caplog.at_level(logging.INFO):
        try:
            fetch_iss_records("iss")
        except RuntimeError as exc:
            assert str(exc) == "fallback timed out"
        else:
            raise AssertionError("expected fallback failure")

    assert "ISS fetch failed via wheretheiss.at: primary unavailable" in caplog.text
    assert "ISS fetch failed via CelesTrak fallback: fallback timed out" in caplog.text
