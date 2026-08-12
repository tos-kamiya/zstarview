from __future__ import annotations

import sys
from datetime import datetime, timezone
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

from zstarview.meteors.projection import project_meteor_observations_to_celestial
from zstarview.meteors.types import GmnLoadResult, MeteorObservation


def _load_probe_module():
    module_path = (
        Path(__file__).resolve().parents[1]
        / "dev-samples"
        / "gmn_meteor_trail_probe.py"
    )
    spec = spec_from_file_location("gmn_meteor_trail_probe", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("failed to load gmn_meteor_trail_probe.py")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _overhead_observation() -> MeteorObservation:
    return MeteorObservation(
        trajectory_id="overhead<&>",
        beginning_utc=datetime(2026, 8, 10, 12, tzinfo=timezone.utc),
        begin_lat_deg=35.0,
        begin_lon_deg=139.0,
        begin_height_km=105.0,
        end_lat_deg=35.0,
        end_lon_deg=139.0,
        end_height_km=85.0,
        duration_s=0.5,
        peak_abs_magnitude=-1.0,
        shower_code="PER",
    )


def test_parse_utc_requires_explicit_timezone() -> None:
    probe = _load_probe_module()

    with pytest.raises(ValueError, match="UTC offset"):
        probe.parse_utc("2026-08-10T12:00:00")
    assert probe.parse_utc("2026-08-10T21:00:00+09:00") == datetime(
        2026, 8, 10, 12, tzinfo=timezone.utc
    )


def test_probe_reports_small_icrs_roundtrip_error() -> None:
    probe = _load_probe_module()
    observation = _overhead_observation()
    trail = project_meteor_observations_to_celestial(
        (observation,),
        observer_lat=35.0,
        observer_lon=139.0,
        observer_height_m=0.0,
    )[0]

    record = probe.build_trail_diagnostic(
        observation,
        trail,
        observer_lat=35.0,
        observer_lon=139.0,
        observer_height_m=0.0,
    )

    assert record["roundtrip_error_deg"]["maximum"] < 1e-5
    assert not record["horizon_clipped"]


def test_probe_counts_and_svg_use_provisional_red() -> None:
    probe = _load_probe_module()
    observation = _overhead_observation()
    display_time = datetime(2026, 8, 10, 13, tzinfo=timezone.utc)
    loaded = GmnLoadResult(
        observations=(observation,),
        source_files=("sample.txt",),
        unavailable_files=(),
    )

    run = probe.build_probe_run(
        loaded,
        display_time_utc=display_time,
        observer_lat=35.0,
        observer_lon=139.0,
        observer_height_m=0.0,
    )
    svg = probe.build_svg(run, size=400)

    assert run.loaded_count == 1
    assert run.candidate_count == 1
    assert run.visible_count == 1
    assert '#f04848' in svg
    assert 'r="2.2"' in svg
    assert "overhead&lt;&amp;&gt;" in svg


def test_polar_svg_places_cardinal_horizon_points() -> None:
    probe = _load_probe_module()

    assert probe.polar_svg_xy(0.0, 0.0, center=100.0, radius=80.0) == pytest.approx(
        (100.0, 20.0)
    )
    assert probe.polar_svg_xy(0.0, 90.0, center=100.0, radius=80.0) == pytest.approx(
        (180.0, 100.0)
    )
