from __future__ import annotations

import numpy as np
import polars as pl

from zstarview.astro import prepare_star_catalog_arrays


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
