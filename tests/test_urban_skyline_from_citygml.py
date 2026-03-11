from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "src" / "zstarview" / "data" / "urban_skyline_from_citygml.py"
    spec = importlib.util.spec_from_file_location("urban_skyline_from_citygml", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_parse_pos_list_ring_reads_lonlat_pairs() -> None:
    mod = _load_module()

    ring = mod.parse_pos_list_ring(
        "35.0 139.0 10 35.0 139.001 10 35.001 139.001 10 35.0 139.0 10"
    )

    assert ring == (
        (139.0, 35.0),
        (139.001, 35.0),
        (139.001, 35.001),
        (139.0, 35.0),
    )


def test_select_tile_envelopes_filters_by_observer_radius(tmp_path: Path) -> None:
    mod = _load_module()
    near = tmp_path / "near.gml"
    far = tmp_path / "far.gml"
    near.write_text(
        (
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            "<core:CityModel xmlns:core=\"http://www.opengis.net/citygml/2.0\" "
            "xmlns:gml=\"http://www.opengis.net/gml\">\n"
            "  <gml:boundedBy>\n"
            "    <gml:Envelope srsName=\"http://www.opengis.net/def/crs/EPSG/0/6697\" srsDimension=\"3\">\n"
            "      <gml:lowerCorner>35.0000 139.0000 0</gml:lowerCorner>\n"
            "      <gml:upperCorner>35.0100 139.0100 30</gml:upperCorner>\n"
            "    </gml:Envelope>\n"
            "  </gml:boundedBy>\n"
            "</core:CityModel>\n"
        ),
        encoding="utf-8",
    )
    far.write_text(
        (
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            "<core:CityModel xmlns:core=\"http://www.opengis.net/citygml/2.0\" "
            "xmlns:gml=\"http://www.opengis.net/gml\">\n"
            "  <gml:boundedBy>\n"
            "    <gml:Envelope srsName=\"http://www.opengis.net/def/crs/EPSG/0/6697\" srsDimension=\"3\">\n"
            "      <gml:lowerCorner>35.2000 139.2000 0</gml:lowerCorner>\n"
            "      <gml:upperCorner>35.2100 139.2100 30</gml:upperCorner>\n"
            "    </gml:Envelope>\n"
            "  </gml:boundedBy>\n"
            "</core:CityModel>\n"
        ),
        encoding="utf-8",
    )

    selected = mod.select_tile_envelopes(
        tmp_path,
        observer_lat_deg=35.005,
        observer_lon_deg=139.005,
        radius_km=2.0,
    )

    assert [envelope.path.name for envelope in selected] == ["near.gml"]


def test_combine_tile_results_uses_max_altitude_per_azimuth() -> None:
    mod = _load_module()
    tower = mod.select_tower("Tokyo Skytree")
    tile_results = (
        mod.TileSkyline(
            envelope=mod.TileEnvelope(Path("a.gml"), 0.0, 0.0, 0.0, 0.0),
            radius_results=(
                mod.SkylineRadiusResult(
                    radius_km=0.1,
                    result=mod.SkylineResult(
                        tower=tower,
                        samples=(
                            mod.SkylineSample(0.0, -10.0),
                            mod.SkylineSample(1.0, -20.0),
                        ),
                        buildings_considered=1,
                        buildings_contributing=1,
                        peak_altitude_deg=-10.0,
                        peak_azimuth_deg=0.0,
                    ),
                ),
            ),
        ),
        mod.TileSkyline(
            envelope=mod.TileEnvelope(Path("b.gml"), 0.0, 0.0, 0.0, 0.0),
            radius_results=(
                mod.SkylineRadiusResult(
                    radius_km=0.1,
                    result=mod.SkylineResult(
                        tower=tower,
                        samples=(
                            mod.SkylineSample(0.0, -15.0),
                            mod.SkylineSample(1.0, -5.0),
                        ),
                        buildings_considered=2,
                        buildings_contributing=1,
                        peak_altitude_deg=-5.0,
                        peak_azimuth_deg=1.0,
                    ),
                ),
            ),
        ),
    )

    radius_results = mod.combine_tile_results(
        tower,
        tile_results,
        radii_km=(0.1,),
        azimuth_step_deg=1.0,
    )
    result = radius_results[0].result

    samples_by_az = {sample.azimuth_deg: sample.altitude_deg for sample in result.samples[:2]}
    assert math.isclose(samples_by_az[0.0], -10.0, abs_tol=1e-9)
    assert math.isclose(samples_by_az[1.0], -5.0, abs_tol=1e-9)
    assert result.buildings_considered == 3
    assert result.buildings_contributing == 2


def test_compute_tile_skylines_runs_sequentially(monkeypatch) -> None:
    mod = _load_module()
    tower = mod.select_tower("Tokyo Tower")
    envelopes = (
        mod.TileEnvelope(Path("a.gml"), 0.0, 0.0, 0.0, 0.0),
        mod.TileEnvelope(Path("b.gml"), 0.0, 0.0, 0.0, 0.0),
    )
    called = []

    def fake_compute_tile_skyline(envelope, **kwargs):
        called.append((envelope.path.name, kwargs["tower"].name))
        return mod.TileSkyline(
            envelope=envelope,
            radius_results=(
                mod.SkylineRadiusResult(
                    radius_km=1.0,
                    result=mod.SkylineResult(
                        tower=tower,
                        samples=(mod.SkylineSample(0.0, -10.0),),
                        buildings_considered=1,
                        buildings_contributing=1,
                        peak_altitude_deg=-10.0,
                        peak_azimuth_deg=0.0,
                    ),
                ),
            ),
        )

    monkeypatch.setattr(mod, "compute_tile_skyline", fake_compute_tile_skyline)
    monkeypatch.setattr(mod, "print_tile_summary", lambda tile_result: None)

    got = mod.compute_tile_skylines(
        envelopes,
        tower=tower,
        cumulative_radii_km=(1.0,),
        radius_band_width_m=30.0,
        azimuth_step_deg=0.1,
        edge_sample_step_m=5.0,
        workers=1,
    )

    assert [item.envelope.path.name for item in got] == ["a.gml", "b.gml"]
    assert called == [("a.gml", "Tokyo Tower"), ("b.gml", "Tokyo Tower")]


def test_compute_tile_skylines_uses_process_pool(monkeypatch) -> None:
    mod = _load_module()
    tower = mod.select_tower("Tokyo Tower")
    envelopes = (
        mod.TileEnvelope(Path("a.gml"), 0.0, 0.0, 0.0, 0.0),
        mod.TileEnvelope(Path("b.gml"), 0.0, 0.0, 0.0, 0.0),
    )
    submitted = []

    class FakeFuture:
        def __init__(self, value):
            self._value = value

        def result(self):
            return self._value

    class FakeExecutor:
        def __init__(self, *, max_workers):
            self.max_workers = max_workers

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def submit(self, fn, envelope, **kwargs):
            submitted.append(
                (
                    self.max_workers,
                    envelope.path.name,
                    kwargs["cumulative_radii_km"],
                    kwargs["radius_band_width_m"],
                )
            )
            return FakeFuture(
                mod.TileSkyline(
                    envelope=envelope,
                    radius_results=(
                        mod.SkylineRadiusResult(
                            radius_km=1.0,
                            result=mod.SkylineResult(
                                tower=tower,
                                samples=(mod.SkylineSample(0.0, -10.0),),
                                buildings_considered=1,
                                buildings_contributing=1,
                                peak_altitude_deg=-10.0,
                                peak_azimuth_deg=0.0,
                            ),
                        ),
                    ),
                )
            )

    monkeypatch.setattr(mod, "ProcessPoolExecutor", FakeExecutor)
    monkeypatch.setattr(mod, "print_tile_summary", lambda tile_result: None)

    got = mod.compute_tile_skylines(
        envelopes,
        tower=tower,
        cumulative_radii_km=(1.0, 2.0),
        radius_band_width_m=30.0,
        azimuth_step_deg=0.1,
        edge_sample_step_m=5.0,
        workers=4,
    )

    assert [item.envelope.path.name for item in got] == ["a.gml", "b.gml"]
    assert submitted == [
        (2, "a.gml", (1.0, 2.0), 30.0),
        (2, "b.gml", (1.0, 2.0), 30.0),
    ]
