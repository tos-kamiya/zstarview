from __future__ import annotations

from dataclasses import FrozenInstanceError
from datetime import timedelta

import numpy as np
import pytest

from zstarview.gui.sky_worker import SkyComputationRequest
from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.types import ScreenGeometry, ViewerData


def _request_inputs() -> dict[str, object]:
    return {
        "ephemeris": object(),
        "viewer_data": ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
        ),
        "geometry": ScreenGeometry(center=(100, 100), radius=90),
        "star_catalog": object(),
        "dso_catalog": None,
        "star_vmag_limit": 6.0,
        "star_subset_indices": np.array([1, 2, 3], dtype=np.int64),
        "delta_t": timedelta(0),
        "sky_update_interval": 60.0,
        "sky_disc_alpha": 0.15,
        "theme": THEME_STYLES_BY_PRESET["night"],
        "star_catalog_meta": None,
        "image_size": [800, 600],
        "sky_disc_render_scale": 0.5,
        "terrain_horizon_profile_altaz": [(1.0, 2.0)],
        "terrain_horizon_profile_distances_m": [100.0],
        "terrain_secondary_ridges_altaz_layers": [[(3.0, 4.0)]],
        "terrain_secondary_ridges_distances_m_layers": [[200.0]],
        "terrain_sample_distances_m": np.array([100.0]),
        "terrain_sample_terrain_elevation_m": np.array([12.0]),
        "night_light_glow_profile": None,
        "night_light_opacity": 0.1,
        "render_generation": 7,
    }


def test_sky_computation_request_snapshots_mutable_inputs() -> None:
    inputs = _request_inputs()
    subset = inputs["star_subset_indices"]
    assert isinstance(subset, np.ndarray)
    request = SkyComputationRequest.from_inputs(**inputs)

    subset[0] = 99
    inputs["terrain_horizon_profile_altaz"].append((5.0, 6.0))  # type: ignore[union-attr]

    assert request.star_subset_indices is not subset
    assert request.star_subset_indices is not None
    assert request.star_subset_indices.tolist() == [1, 2, 3]
    assert request.image_size == (800, 600)
    assert request.terrain_horizon_profile_altaz == ((1.0, 2.0),)
    assert request.terrain_horizon_profile_distances_m == (100.0,)
    assert request.terrain_sample_distances_m is not None
    assert not request.terrain_sample_distances_m.flags.writeable


def test_sky_computation_request_is_frozen() -> None:
    request = SkyComputationRequest.from_inputs(**_request_inputs())

    with pytest.raises(FrozenInstanceError):
        request.render_generation = 8  # type: ignore[misc]
