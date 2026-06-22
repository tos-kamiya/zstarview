# -*- coding: utf-8 -*-
"""Tests for alt/az coordinate projection helpers."""

from __future__ import annotations

import numpy as np

from zstarview.clouddisc.altaz_projection import (
    altaz_to_bin_indices,
    geodetic_to_altaz,
    geodetic_to_altaz_array,
)


def test_geodetic_to_altaz_zenith():
    # A point directly above the observer maps to (alt=90, az=any).
    alt, az = geodetic_to_altaz(35.0, 135.0, 6374.0, 35.0, 135.0)
    assert abs(alt - 90.0) < 1e-6


def test_geodetic_to_altaz_north():
    # A point slightly north of the observer is near the northern horizon.
    alt, az = geodetic_to_altaz(36.0, 135.0, 6374.0, 35.0, 135.0)
    assert 0.0 < alt < 5.0
    assert abs(az - 0.0) < 1.0 or abs(az - 360.0) < 1.0


def test_geodetic_to_altaz_east():
    # A point slightly east of the observer is near the eastern horizon.
    alt, az = geodetic_to_altaz(35.0, 136.0, 6374.0, 35.0, 135.0)
    assert 0.0 < alt < 5.0
    assert abs(az - 90.0) < 1.0


def test_geodetic_to_altaz_array_matches_scalar():
    lat_arr = np.array([[35.0, 35.1], [35.0, 35.2]])
    lon_arr = np.array([[135.0, 135.0], [136.0, 135.0]])
    alt_arr, az_arr = geodetic_to_altaz_array(lat_arr, lon_arr, 6374.0, 35.0, 135.0)

    for i in range(lat_arr.shape[0]):
        for j in range(lat_arr.shape[1]):
            alt_s, az_s = geodetic_to_altaz(lat_arr[i, j], lon_arr[i, j], 6374.0, 35.0, 135.0)
            assert abs(alt_arr[i, j] - alt_s) < 1e-5
            # Azimuth is ill-defined near the zenith; allow a looser tolerance there.
            az_tol = 0.5 if alt_s > 89.9 else 1e-5
            assert abs(az_arr[i, j] - az_s) < az_tol


def test_altaz_to_bin_indices_wraps_azimuth():
    alt = np.array([0.0, 45.0, 90.0])
    az = np.array([0.0, 180.0, 359.0])
    alt_idx, az_idx = altaz_to_bin_indices(alt, az, alt_bins=90, az_bins=720)
    assert alt_idx.tolist() == [0, 45, 89]
    assert az_idx[0] == 0
    assert 350 < az_idx[2] <= 719


def test_altaz_to_bin_indices_clips_altitude():
    alt = np.array([-10.0, 100.0])
    az = np.array([0.0, 0.0])
    alt_idx, az_idx = altaz_to_bin_indices(alt, az, alt_bins=90, az_bins=720)
    assert alt_idx.tolist() == [0, 89]
