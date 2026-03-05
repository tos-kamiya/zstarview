from __future__ import annotations

import math

from zstarview.paths import CLOUD_MISSING_TINT_RGBA
from zstarview.zstarview import parse_args


def test_parse_args_cloud_missing_tint_opacity_default(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    args = parse_args()
    expected = float(CLOUD_MISSING_TINT_RGBA[3]) / 255.0
    assert math.isclose(float(args.cloud_missing_tint_opacity), expected, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_cloud_missing_tint_opacity_override(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--cloud-missing-tint-opacity", "0.4"])
    args = parse_args()
    assert math.isclose(float(args.cloud_missing_tint_opacity), 0.4, rel_tol=0.0, abs_tol=1e-9)


def test_parse_args_overlay_visibility_defaults(monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview"])
    args = parse_args()
    assert args.show_dso_initial is None
    assert args.show_asterisms_initial is None


def test_parse_args_overlay_visibility_override(monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        ["zstarview", "--show-dso-initial", "false", "--show-asterisms-initial", "true"],
    )
    args = parse_args()
    assert args.show_dso_initial is False
    assert args.show_asterisms_initial is True
