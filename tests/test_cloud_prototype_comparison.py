from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace

import numpy as np
import pytest
from PIL import Image

from scripts import prototype_cloud_cell_shading as prototype


def test_voxel_variants_select_only_the_requested_lights(monkeypatch) -> None:
    calls = []

    def shade(*args, **kwargs):
        calls.append(kwargs)
        return args[7], {}

    monkeypatch.setattr(prototype, "shade_native_voxels", shade)
    sampler = object()
    for mode in prototype.NIGHT_LIGHTING_MODES:
        prototype._render_voxel_variant(
            object(), 35, 135, np.array([45]), np.array([120]), -20, 220,
            np.zeros((1, 3)), mode=mode, moon_altaz_phase=(45, 120, 180),
            ground_sampler=sampler, opacity=0.85, height_layer_transform=True,
            cloud_amount_subtract=0.09,
        )
    assert [call["moon_altaz_phase"] is not None for call in calls] == [False, True, False, True]
    assert [call["ground_light_sampler"] is sampler for call in calls] == [False, False, True, True]
    assert all(call["cloud_amount_subtract"] == 0.09 for call in calls)


def test_ground_sample_cache_requires_identical_geometry(monkeypatch) -> None:
    calls = []

    def sample(*args):
        calls.append(args)
        return np.ones(args[2], dtype=np.float32)

    monkeypatch.setattr(prototype, "sample_cloud_ground_light", sample)
    sampler = prototype._CachedGroundLightSampler()
    args = (35.0, 135.0, (1, 1, 1), np.zeros(2), np.zeros((1, 2)), np.eye(2)[None])
    first = sampler(*args)
    assert sampler(*args) is first
    shifted = args[:4] + (np.ones((1, 2)), args[5])
    assert sampler(*shifted) is not first
    assert len(calls) == 2


@pytest.mark.parametrize("cloud_only", [True, False])
@pytest.mark.parametrize("compare", [True, False])
@pytest.mark.parametrize("moon_override", ["none", "all", "phase"])
def test_comparison_uses_one_scene_and_saves_variants_with_metadata(
    monkeypatch, tmp_path, cloud_only, compare, moon_override
) -> None:
    when = datetime(2026, 10, 10, 12, tzinfo=timezone.utc)
    source = SimpleNamespace(
        time_utc=when - timedelta(minutes=15), satellite="Test", product="B13"
    )
    source_calls = []
    sky_calls = []
    ground_calls = []
    shade_calls = []
    base_colors = []
    monkeypatch.setattr(prototype, "resolve_launch_location", lambda _: SimpleNamespace(
        lat=35.0, lon=135.0, display_name="Test location"
    ))
    actual_moon = (45.0, 120.0, 180.0) if moon_override == "none" else (5.0, 60.0, 5.0)
    monkeypatch.setattr(prototype, "_celestial_lighting", lambda *args: (
        (-20.0, 220.0), actual_moon
    ))
    monkeypatch.setattr(prototype, "CloudDisc", lambda config: object())

    def fetch(*args, **kwargs):
        source_calls.append(args)
        return source

    def sky(**kwargs):
        sky_calls.append(kwargs)
        Image.new("RGBA", kwargs["size"], (10, 20, 30, 255)).save(kwargs["output"])

    def ground(*args):
        ground_calls.append(args)
        return np.ones(args[2], dtype=np.float32)

    def shade(source_arg, lat, lon, alt, az, sun_alt, sun_az, base, **kwargs):
        assert source_arg is source
        assert (sun_alt, sun_az) == (-20.0, 220.0)
        shade_calls.append(kwargs)
        base_colors.append(base.copy())
        ground_max = 0.0
        if kwargs["ground_light_sampler"] is not None:
            sampled = kwargs["ground_light_sampler"](
                lat, lon, (1, 1, 1), np.zeros(2), np.zeros((1, 2)), np.eye(2)[None]
            )
            ground_max = float(sampled.max())
        moon = kwargs["moon_altaz_phase"] is not None
        return base + 0.05 + 0.10 * moon + 0.20 * ground_max, {
            "coverage_ratio": 1.0, "moon_altaz_phase": kwargs["moon_altaz_phase"],
            "ground_light_max": ground_max,
        }

    monkeypatch.setattr(prototype, "run_cloud_source_worker_process", fetch)
    monkeypatch.setattr(prototype, "_base_sky_image", sky)
    monkeypatch.setattr(prototype, "sample_cloud_ground_light", ground)
    monkeypatch.setattr(prototype, "shade_native_voxels", shade)
    output = tmp_path / "compare.png"
    args = [
        "@35,135", "--datetime", when.isoformat(),
        "--size", "64x64", "--output", str(output),
    ]
    if compare:
        args.append("--compare-night-lighting")
    expected_moon = actual_moon
    if moon_override == "all":
        args.extend(["--moon-alt", "45", "--moon-az", "360", "--moon-phase", "180"])
        expected_moon = (45.0, 0.0, 180.0)
    elif moon_override == "phase":
        args.extend(["--moon-phase", "180"])
        expected_moon = (5.0, 60.0, 180.0)
    if cloud_only:
        args.append("--cloud-only")
    assert prototype.main(args) == 0
    assert len(source_calls) == 1
    assert len(sky_calls) == (0 if cloud_only else 1)
    assert len(ground_calls) == 1
    assert len(shade_calls) == (4 if compare else 1)
    for base in base_colors:
        np.testing.assert_array_equal(base, base_colors[0])
    with Image.open(output) as sheet:
        assert sheet.size == ((152, 216) if compare else (64, 64))
    comparison = json.loads(output.with_suffix(".json").read_text())
    assert comparison["night_lighting_mode"] == ("comparison" if compare else "combined")
    if compare:
        assert [panel["mode"] for panel in comparison["panels"]] == list(prototype.NIGHT_LIGHTING_MODES)
    brightness = []
    modes = prototype.NIGHT_LIGHTING_MODES if compare else ("combined",)
    for mode in modes:
        path = tmp_path / f"compare-{mode}.png" if compare else output
        with Image.open(path) as panel:
            brightness.append(np.asarray(panel)[32, 32, 0])
        metadata = json.loads(path.with_suffix(".json").read_text())
        assert metadata["night_lighting_mode"] == mode
        assert metadata["render_time_utc"] == "2026-10-10T12:00:00Z"
        assert metadata["cloud_observation_time_utc"] == "2026-10-10T11:45:00Z"
        assert metadata["moon"]["solar_elongation_deg"] == 180.0
        assert metadata["moon"]["alt_deg"] == expected_moon[0]
        assert metadata["moon"]["az_deg"] == expected_moon[1]
        assert metadata["moon"]["override_applied"] is (moon_override != "none")
        assert metadata["moon"]["actual"] == {
            "alt_deg": actual_moon[0], "az_deg": actual_moon[1],
            "solar_elongation_deg": actual_moon[2],
        }
        assert metadata["model"]["cloud_amount_subtract"] == 0.09
    for call in shade_calls:
        if call["moon_altaz_phase"] is not None:
            assert call["moon_altaz_phase"] == expected_moon
    if compare:
        assert brightness[0] < brightness[1] < brightness[2] < brightness[3]


@pytest.mark.parametrize("option", ["--compare-night-lighting", "--moon-phase=180"])
def test_comparison_rejects_legacy_cloud_model_before_fetching(tmp_path, option) -> None:
    with pytest.raises(SystemExit) as exc:
        prototype.main([
            "Tokyo", "--cloud-model", "altaz", option,
            "--output", str(tmp_path / "compare.png"),
        ])
    assert exc.value.code == 2


@pytest.mark.parametrize("option", [
    "--moon-alt=91", "--moon-alt=-91", "--moon-alt=nan", "--moon-alt=inf",
    "--moon-az=-1", "--moon-az=361", "--moon-az=nan", "--moon-az=-inf",
    "--moon-phase=-1", "--moon-phase=181", "--moon-phase=nan", "--moon-phase=inf",
])
def test_invalid_moon_override_is_rejected_before_acquisition(
    monkeypatch, tmp_path, capsys, option
) -> None:
    def forbidden(*args):
        raise AssertionError("invalid input must not start celestial or cloud acquisition")
    monkeypatch.setattr(prototype, "_celestial_lighting", forbidden)
    monkeypatch.setattr(prototype, "resolve_launch_location", forbidden)
    with pytest.raises(SystemExit) as exc:
        prototype.main(["Tokyo", option, "--output", str(tmp_path / "moon.png")])
    assert exc.value.code == 2
    assert "must be finite and between" in capsys.readouterr().err
