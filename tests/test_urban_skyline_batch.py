from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "src" / "zstarview" / "data" / "urban_skyline_batch.py"
    spec = importlib.util.spec_from_file_location("urban_skyline_batch", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_discover_citygml_dirs_finds_udx_bldg(tmp_path: Path) -> None:
    mod = _load_module()
    bldg_dir = tmp_path / "13100_tokyo" / "udx" / "bldg"
    bldg_dir.mkdir(parents=True)

    got = mod.discover_citygml_dirs(roots=(tmp_path,), explicit_dirs=())

    assert got == (bldg_dir.resolve(),)


def test_select_candidate_coverages_filters_by_distance() -> None:
    mod = _load_module()
    tower = next(tower for tower in mod.load_tower_viewpoints() if tower.name == "Tokyo Skytree")
    coverages = (
        mod.CityGmlCoverage(
            bldg_dir=Path("/tmp/tokyo"),
            min_lat_deg=35.6,
            min_lon_deg=139.7,
            max_lat_deg=35.8,
            max_lon_deg=139.9,
        ),
        mod.CityGmlCoverage(
            bldg_dir=Path("/tmp/osaka"),
            min_lat_deg=34.6,
            min_lon_deg=135.4,
            max_lat_deg=34.8,
            max_lon_deg=135.6,
        ),
    )

    got = mod.select_candidate_coverages(tower, coverages, radius_km=30.0)

    assert [coverage.bldg_dir.name for coverage in got] == ["tokyo"]


def test_select_towers_returns_only_covered_japan_towers(monkeypatch) -> None:
    mod = _load_module()
    tower = next(tower for tower in mod.load_tower_viewpoints() if tower.name == "Tokyo Tower")
    outside = next(t for t in mod.load_tower_viewpoints() if t.name == "CN Tower")
    monkeypatch.setattr(mod, "load_tower_viewpoints", lambda: (tower, outside))

    coverages = (
        mod.CityGmlCoverage(
            bldg_dir=Path("/tmp/tokyo"),
            min_lat_deg=35.5,
            min_lon_deg=139.6,
            max_lat_deg=35.8,
            max_lon_deg=139.9,
        ),
    )

    got = mod.select_towers(
        tower_queries=(),
        all_covered_towers=True,
        coverages=coverages,
        radius_km=30.0,
    )

    assert [item.name for item in got] == ["Tokyo Tower"]


def test_compute_tower_skyline_skips_coverages_without_matching_tiles(monkeypatch) -> None:
    mod = _load_module()
    tower = next(tower for tower in mod.load_tower_viewpoints() if tower.name == "Tokyo Skytree")
    coverages = (
        mod.CityGmlCoverage(
            bldg_dir=Path("/tmp/near"),
            min_lat_deg=35.6,
            min_lon_deg=139.7,
            max_lat_deg=35.8,
            max_lon_deg=139.9,
        ),
        mod.CityGmlCoverage(
            bldg_dir=Path("/tmp/far"),
            min_lat_deg=35.0,
            min_lon_deg=138.0,
            max_lat_deg=35.1,
            max_lon_deg=138.1,
        ),
    )

    def fake_select_tile_envelopes(citygml_dir, **kwargs):
        if Path(citygml_dir).name == "near":
            return (
                mod.TileEnvelope(
                    path=Path("/tmp/near/tile.gml"),
                    min_lat_deg=35.7,
                    min_lon_deg=139.8,
                    max_lat_deg=35.71,
                    max_lon_deg=139.81,
                ),
            )
        raise ValueError("No CityGML building tiles found within 19.6 km.")

    monkeypatch.setattr(mod, "select_tile_envelopes", fake_select_tile_envelopes)
    monkeypatch.setattr(
        mod,
        "compute_tile_skylines",
        lambda envelopes, **kwargs: [
            type(
                "TileResult",
                (),
                {
                    "radius_results": (),
                },
            )()
        ]
        if envelopes
        else [],
    )
    monkeypatch.setattr(mod, "combine_tile_results", lambda *args, **kwargs: ("ok",))

    got = mod.compute_tower_skyline(
        tower,
        coverages=coverages,
        cumulative_radii_km=(0.1, 0.15),
        radius_km=19.6,
        radius_band_width_m=90.0,
        azimuth_step_deg=0.1,
        edge_sample_step_m=5.0,
        workers=1,
        min_building_height_m=50.0,
    )

    assert got == ("ok",)
