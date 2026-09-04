import astropy.time
import numpy as np

from zstarview import astro


def _star_catalog() -> astro.StarCatalogArrays:
    return {
        "catalog_index": np.array([0, 1, 2], dtype=np.int32),
        "unit_vectors": np.zeros((3, 3), dtype=float),
        "vmag": np.array([1.0, 2.0, 3.0], dtype=float),
        "bv": np.array([0.1, 0.2, 0.3], dtype=float),
        "size_scale": np.array([1.0, 1.0, 1.0], dtype=float),
        "color_base": np.array([1.0, 1.0, 1.0], dtype=float),
    }


def test_positional_static_star_data_policy_keeps_all_magnitude_selected_stars(
    monkeypatch,
) -> None:
    monkeypatch.setattr(astro, "build_icrs_to_altaz_matrix", lambda *_args: object())
    monkeypatch.setattr(
        astro,
        "apply_icrs_to_altaz_matrix",
        lambda *_args: (
            np.array([10.0, 20.0, 30.0], dtype=float),
            np.array([0.0, 120.0, 240.0], dtype=float),
        ),
    )
    monkeypatch.setattr(
        astro,
        "is_in_fov_vectorized",
        lambda alt, *_args, **_kwargs: np.array([True, False, True]),
    )

    stars, _location = astro.calculate_visible_stars(
        _star_catalog(),
        0.0,
        0.0,
        0.0,
        astropy.time.Time("2026-07-08T00:00:00", scale="utc"),
    )

    assert stars["star_index"].tolist() == [0, 1, 2]


def test_scenic_view_scoped_star_data_policy_keeps_all_stars(monkeypatch) -> None:
    monkeypatch.setattr(astro, "build_icrs_to_altaz_matrix", lambda *_args: object())
    monkeypatch.setattr(
        astro,
        "apply_icrs_to_altaz_matrix",
        lambda *_args: (
            np.array([10.0, 20.0, 30.0], dtype=float),
            np.array([0.0, 120.0, 240.0], dtype=float),
        ),
    )
    monkeypatch.setattr(
        astro,
        "is_in_fov_vectorized",
        lambda alt, *_args, **_kwargs: np.array([True, False, True]),
    )

    stars, _location = astro.calculate_visible_stars(
        _star_catalog(),
        0.0,
        0.0,
        0.0,
        astropy.time.Time("2026-07-08T00:00:00", scale="utc"),
    )

    assert stars["star_index"].tolist() == [0, 1, 2]
