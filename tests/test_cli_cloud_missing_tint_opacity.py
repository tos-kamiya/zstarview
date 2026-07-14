from __future__ import annotations

import math

import pytest

from zstarview.cli.args import parse_args
from zstarview.paths import CLOUD_MISSING_TINT_RGBA


def test_parse_args_cloud_missing_tint_opacity_default(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    args = parse_args()
    expected = float(CLOUD_MISSING_TINT_RGBA[3]) / 255.0
    assert math.isclose(float(args.cloud_missing_tint_opacity), expected, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_cloud_missing_tint_opacity_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--cloud-missing-tint-opacity", "0.4"])
    args = parse_args()
    assert math.isclose(float(args.cloud_missing_tint_opacity), 0.4, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_cloud_stripe_default(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    args = parse_args()
    assert args.cloud_stripe == ("width", 50, 0.85)


def test_parse_args_cloud_opacity_default(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    args = parse_args()
    assert math.isclose(float(args.cloud_opacity), 0.10, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_cloud_opacity_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--cloud-opacity", "0.1"])
    args = parse_args()
    assert math.isclose(float(args.cloud_opacity), 0.1, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_cloud_stripe_accepts_mode_only(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--cloud-stripe", "alpha"])
    args = parse_args()
    assert args.cloud_stripe == ("alpha", 50, 0.25)


def test_parse_args_cloud_stripe_accepts_halftone_mode_only(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--cloud-stripe", "halftone"])
    args = parse_args()
    assert args.cloud_stripe == ("halftone", 30, 1.7)


def test_parse_args_cloud_stripe_accepts_mode_and_count(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--cloud-stripe", "width,24"])
    args = parse_args()
    assert args.cloud_stripe == ("width", 24, 0.85)


def test_parse_args_cloud_stripe_accepts_mode_count_and_width(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--cloud-stripe", "alpha,12,0.35"])
    args = parse_args()
    assert args.cloud_stripe == ("alpha", 12, 0.35)


def test_parse_args_overlay_visibility_defaults(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    args = parse_args()
    assert args.show_dso_initial is None
    assert args.show_asterisms_initial is None
    assert args.show_guidelines_initial is None
    assert args.observation_info == "bottom"


def test_parse_args_overlay_visibility_override(monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        ["zstarview", "--show-dso-initial", "false", "--show-asterisms-initial", "true"],
    )
    args = parse_args()
    assert args.show_dso_initial is False
    assert args.show_asterisms_initial is True


def test_parse_args_show_guidelines_initial_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-guidelines-initial", "false"])
    args = parse_args()
    assert args.show_guidelines_initial is False


def test_parse_args_observation_info_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--observation-info", "off"])
    args = parse_args()
    assert args.observation_info == "off"


def test_parse_args_observer_height_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--observer-height-m", "123.4"])
    args = parse_args()
    assert math.isclose(float(args.observer_height_m), 123.4, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_use_building_top(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--use-building-top"])
    args = parse_args()
    assert args.use_building_top is True


def test_parse_args_urban_outline_opacity_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--urban-outline-opacity", "0.6"])
    args = parse_args()
    assert math.isclose(float(args.urban_outline_opacity), 0.6, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_urban_outline_opacity_default(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    args = parse_args()
    assert math.isclose(float(args.urban_outline_opacity), 0.2, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_satellite_opacity_default(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    args = parse_args()
    assert math.isclose(float(args.satellite_opacity), 0.7, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_satellite_opacity_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--satellite-opacity", "0.3"])
    args = parse_args()
    assert math.isclose(float(args.satellite_opacity), 0.3, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_tropical_cyclone_opacity_default(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    args = parse_args()
    assert math.isclose(float(args.tropical_cyclone_opacity), 0.7, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_tropical_cyclone_opacity_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--tropical-cyclone-opacity", "0.25"])
    args = parse_args()
    assert math.isclose(float(args.tropical_cyclone_opacity), 0.25, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_urban_outline_radius_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--urban-outline-radius-km", "1.8"])
    args = parse_args()
    assert math.isclose(float(args.urban_outline_radius_km), 1.8, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_place_short_option(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "-p", "Tokyo Tower"])
    args = parse_args()
    assert args.place == "Tokyo Tower"


def test_parse_args_urban_outline_radius_short_option(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "-r", "2.2"])
    args = parse_args()
    assert math.isclose(float(args.urban_outline_radius_km), 2.2, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_urban_outline_skyscraper_radius_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--urban-outline-skyscraper-radius-km", "48"])
    args = parse_args()
    assert math.isclose(float(args.urban_outline_skyscraper_radius_km), 48.0, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_urban_outline_skyscraper_radius_zero_disables_lookup(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--urban-outline-skyscraper-radius-km", "0"])
    args = parse_args()
    assert math.isclose(float(args.urban_outline_skyscraper_radius_km), 0.0, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_rejects_skyscraper_radius_smaller_than_base_radius(monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        ["zstarview", "--urban-outline-radius-km", "12", "--urban-outline-skyscraper-radius-km", "8"],
    )
    with pytest.raises(SystemExit):
        parse_args()


def test_parse_args_urban_outline_min_height_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--urban-outline-min-building-height-m", "40"])
    args = parse_args()
    assert math.isclose(float(args.urban_outline_min_height_m), 40.0, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_urban_outline_min_height_short_option(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "-b", "25"])
    args = parse_args()
    assert math.isclose(float(args.urban_outline_min_height_m), 25.0, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_urban_outline_feature_type_default(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    args = parse_args()
    assert args.urban_outline_feature_type == "both"


def test_parse_args_urban_outline_feature_type_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--urban-outline-feature-type", "building"])
    args = parse_args()
    assert args.urban_outline_feature_type == "building"
