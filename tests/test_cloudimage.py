from __future__ import annotations

import argparse
from types import SimpleNamespace

import numpy as np

from zstarview import cloudimage


class _FakeAltAzGrid:
    pass


class _FakeCloudDisc:
    def __init__(self) -> None:
        self.fetch_args: tuple[float, float] | None = None
        self.build_grid_args: dict[str, object] | None = None

    def fetch_source(self, *, lat: float, lon: float):
        self.fetch_args = (lat, lon)
        return SimpleNamespace()

    def build_altaz_grid_from_source(
        self,
        *,
        source,
        lat: float,
        lon: float,
        cloud_shells_km,
    ):
        self.build_grid_args = {
            "lat": lat,
            "lon": lon,
            "cloud_shells_km": cloud_shells_km,
        }
        return _FakeAltAzGrid()


def test_parse_observer_spec_requires_at_prefix() -> None:
    with np.testing.assert_raises_regex(argparse.ArgumentTypeError, "@lat,lon"):
        cloudimage.parse_observer_spec("35.0,139.0")


def test_parse_observer_spec_parses_coordinates() -> None:
    lat, lon = cloudimage.parse_observer_spec("@35.0, 139.0")
    assert lat == 35.0
    assert lon == 139.0


def test_compose_on_black_background_flattens_alpha() -> None:
    cloud = np.array(
        [[[255, 255, 255, 0], [255, 255, 255, 128]]],
        dtype=np.uint8,
    )
    out = cloudimage._compose_on_black_background(cloud)
    assert out.dtype == np.uint8
    assert out.shape == (1, 2, 4)
    assert tuple(out[0, 0]) == (0, 0, 0, 255)
    assert tuple(out[0, 1]) == (128, 128, 128, 255)


def _fake_render_altaz_grid_circles(
    altaz_grid,
    *,
    width,
    height,
    center_alt_deg,
    center_az_deg,
    edge_fov_deg,
    mask_fov_deg,
):
    return np.array(
        [
            [[255, 255, 255, 0], [255, 255, 255, 128]],
            [[255, 255, 255, 255], [255, 255, 255, 64]],
        ],
        dtype=np.uint8,
    )


def test_render_cloud_image_uses_single_fov(monkeypatch) -> None:
    fake = _FakeCloudDisc()
    monkeypatch.setattr(cloudimage, "CloudDisc", lambda _cfg: fake)
    monkeypatch.setattr(
        cloudimage,
        "render_altaz_grid_circles",
        _fake_render_altaz_grid_circles,
    )
    image = cloudimage.render_cloud_image(
        observer_lat=35.0,
        observer_lon=139.0,
        alt=25.0,
        az=180.0,
        fov_deg=60.0,
        radius_px=4,
    )
    assert image.shape == (2, 2, 4)
    assert fake.fetch_args == (35.0, 139.0)
    assert fake.build_grid_args is not None
    assert fake.build_grid_args["lat"] == 35.0
    assert fake.build_grid_args["lon"] == 139.0
