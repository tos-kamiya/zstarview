from __future__ import annotations

from datetime import datetime, timezone
from types import SimpleNamespace

from zstarview.gui.window import DEFAULT_CLOUD_ALT_MIN_DEG, SkyWindow


def _dummy_window(cloud_state):
    return SimpleNamespace(
        cloud_state=cloud_state,
        _predicted_cloud_satellite=lambda: "G19",
    )


def test_cloud_status_line_shows_downloading() -> None:
    state = SimpleNamespace(
        current_satellite="HIMAWARI",
        banner_text="Clouds: downloading…",
        meta=None,
        coverage_ratio=None,
    )
    got = SkyWindow._cloud_status_line(_dummy_window(state))
    assert got == "Clouds [HIMAWARI]: downloading"


def test_cloud_status_line_shows_partial_when_coverage_incomplete() -> None:
    meta = SimpleNamespace(
        satellite="HIMAWARI",
        product="ISatSS-B13",
        time_utc=datetime(2026, 3, 5, 1, 20, tzinfo=timezone.utc),
    )
    state = SimpleNamespace(
        current_satellite=None,
        banner_text=None,
        meta=meta,
        coverage_ratio=0.72,
    )
    got = SkyWindow._cloud_status_line(_dummy_window(state))
    assert got == "Clouds [HIMAWARI]: partial 72% 01:20Z"


def test_cloud_status_line_shows_idle_without_meta_or_banner() -> None:
    state = SimpleNamespace(
        current_satellite=None,
        banner_text=None,
        meta=None,
        coverage_ratio=None,
    )
    got = SkyWindow._cloud_status_line(_dummy_window(state))
    assert got == "Clouds [G19]: idle"


def test_default_cloud_horizon_cutoff_is_one_degree() -> None:
    assert DEFAULT_CLOUD_ALT_MIN_DEG == 1.0
