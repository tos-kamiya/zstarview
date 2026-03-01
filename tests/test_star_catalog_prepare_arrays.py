from __future__ import annotations

import astropy.time
import astropy.units as u
import numpy as np
import polars as pl
from astropy.coordinates import AltAz, EarthLocation, SkyCoord

from zstarview.astro import (
    apply_icrs_to_altaz_matrix,
    build_icrs_to_altaz_matrix,
    prepare_star_catalog_arrays,
)


def test_prepare_star_catalog_arrays_returns_numpy_arrays() -> None:
    df = pl.DataFrame(
        {
            "RAh": [0.0, 1.0],
            "Dec": [10.0, 20.0],
            "Vmag": [5.5, 7.1],
            "B-V": [0.2, None],
            "Name": ["a", "b"],
        }
    )

    out = prepare_star_catalog_arrays(df)

    assert isinstance(out["ra_h"], np.ndarray)
    assert isinstance(out["dec"], np.ndarray)
    assert isinstance(out["vmag"], np.ndarray)
    assert isinstance(out["bv"], np.ndarray)
    assert isinstance(out["name"], np.ndarray)
    assert out["name"].tolist() == ["a", "b"]
    assert np.isnan(out["bv"][1])


def test_prepare_star_catalog_arrays_respects_max_vmag() -> None:
    df = pl.DataFrame(
        {
            "RAh": [0.0, 1.0, 2.0],
            "Dec": [10.0, 20.0, 30.0],
            "Vmag": [5.5, 6.0, 6.1],
            "B-V": [0.2, 0.3, 0.4],
            "Name": ["a", "b", "c"],
        }
    )

    out = prepare_star_catalog_arrays(df, max_vmag=6.0)

    assert out["name"].tolist() == ["a", "b"]
    assert out["vmag"].tolist() == [5.5, 6.0]


def test_prepare_star_catalog_arrays_computes_trig_fields() -> None:
    df = pl.DataFrame(
        {
            "RAh": [0.0, 1.0],
            "Dec": [10.0, 20.0],
            "Vmag": [1.0, 2.0],
            "B-V": [0.1, 0.2],
            "Name": ["a", "b"],
        }
    )

    out = prepare_star_catalog_arrays(df)

    expected_ra_rad = np.radians(np.array([0.0, 1.0]) * 15.0)
    expected_dec_rad = np.radians(np.array([10.0, 20.0]))
    assert np.allclose(out["ra_rad"], expected_ra_rad)
    assert np.allclose(out["dec_rad"], expected_dec_rad)
    assert np.allclose(out["sin_dec"], np.sin(expected_dec_rad))
    assert np.allclose(out["cos_dec"], np.cos(expected_dec_rad))


def test_prepare_star_catalog_arrays_computes_unit_vectors() -> None:
    df = pl.DataFrame(
        {
            "RAh": [0.0, 1.0],
            "Dec": [10.0, 20.0],
            "Vmag": [1.0, 2.0],
            "B-V": [0.1, 0.2],
            "Name": ["a", "b"],
        }
    )

    out = prepare_star_catalog_arrays(df)

    ra_rad = np.radians(np.array([0.0, 1.0]) * 15.0)
    dec_rad = np.radians(np.array([10.0, 20.0]))
    expected_vectors = np.column_stack(
        (
            np.cos(dec_rad) * np.cos(ra_rad),
            np.cos(dec_rad) * np.sin(ra_rad),
            np.sin(dec_rad),
        )
    )

    np.testing.assert_allclose(out["unit_vectors"], expected_vectors)


def test_icrs_to_altaz_matrix_matches_skycoord() -> None:
    df = pl.DataFrame(
        {
            "RAh": [2.0, 5.0, 10.5],
            "Dec": [-10.0, 15.0, 30.5],
            "Vmag": [1.2, 2.3, 3.0],
            "B-V": [0.0, 0.1, 0.2],
            "Name": ["s1", "s2", "s3"],
        }
    )

    arrays = prepare_star_catalog_arrays(df)
    time_obj = astropy.time.Time("2026-03-01T12:30:00", scale="utc")
    lat = 51.0
    lon = 0.127
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

    matrix = build_icrs_to_altaz_matrix(time_obj, location)
    alt, az = apply_icrs_to_altaz_matrix(arrays["unit_vectors"], matrix)

    coords = SkyCoord(ra=np.degrees(arrays["ra_rad"]) * u.deg, dec=np.degrees(arrays["dec_rad"]) * u.deg, frame="icrs")
    altaz = coords.transform_to(AltAz(obstime=time_obj, location=location))
    alt_ref = altaz.alt.deg
    az_ref = altaz.az.deg

    np.testing.assert_allclose(alt, alt_ref, atol=1e-2)
    az_diff = (az - az_ref + 180.0) % 360.0 - 180.0
    np.testing.assert_allclose(az_diff, 0.0, atol=1e-2)
