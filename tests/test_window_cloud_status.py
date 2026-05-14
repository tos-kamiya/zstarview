from __future__ import annotations

from datetime import datetime, timezone
from types import SimpleNamespace

from zstarview.gui.window import DEFAULT_CLOUD_ALT_MIN_DEG, SkyWindow


def _dummy_window(cloud_state):
    return SimpleNamespace(
        cloud_state=cloud_state,
        cloud_disc_alpha=0.2,
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
    assert got == "☁ HIMAWARI downloading"


def test_cloud_status_line_shows_download_failed() -> None:
    state = SimpleNamespace(
        current_satellite="HIMAWARI",
        banner_text="Clouds: Clouds fetch timed out",
        meta=None,
        coverage_ratio=None,
    )
    got = SkyWindow._cloud_status_line(_dummy_window(state))
    assert got == "☁ HIMAWARI failed"


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
    assert got == "☁ HIMAWARI 72% 01:20Z"


def test_cloud_status_line_shows_idle_without_meta_or_banner() -> None:
    state = SimpleNamespace(
        current_satellite=None,
        banner_text=None,
        meta=None,
        coverage_ratio=None,
    )
    got = SkyWindow._cloud_status_line(_dummy_window(state))
    assert got == "☁ G19 idle"


def test_satellite_status_line_shows_epoch_with_icon() -> None:
    state = SimpleNamespace(
        banner_text=None,
        element_epoch_utc=datetime(2026, 3, 5, 1, 20, tzinfo=timezone.utc),
    )
    dummy = SimpleNamespace(satellite_state=state, satellite_opacity=0.2)

    got = SkyWindow._satellite_status_line(dummy)

    assert got == "🛰 01:20Z"


def test_aircraft_status_line_shows_last_success_with_icon() -> None:
    state = SimpleNamespace(
        banner_text=None,
        last_success_utc=datetime(2026, 3, 5, 1, 20, tzinfo=timezone.utc),
    )
    dummy = SimpleNamespace(aircraft_state=state, aircraft_opacity=0.2)

    got = SkyWindow._aircraft_status_line(dummy)

    assert got == "✈ 01:20Z"


def test_terrain_and_urban_status_lines_show_icons() -> None:
    terrain_state = SimpleNamespace(banner_text="Terrain horizon: loading DEM...", current_source="DEM")
    urban_state = SimpleNamespace(banner_text="Urban outline: downloading...", current_source="overture")
    dummy = SimpleNamespace(
        terrain_horizon_state=terrain_state,
        terrain_horizon_opacity=0.2,
        urban_outline_state=urban_state,
        urban_outline_opacity=0.2,
    )

    assert SkyWindow._terrain_horizon_status_line(dummy) == "▲ loading DEM..."
    assert SkyWindow._urban_outline_status_line(dummy) == "🂓 downloading..."


def test_water_status_line_shows_only_count_when_enabled() -> None:
    dummy = SimpleNamespace(
        water_overlay_state=SimpleNamespace(
            banner_text=None,
            points=[object(), object(), object()],
            current_source="Water: overpass",
            current_mode="dem",
        ),
        water_overlay_opacity=0.2,
    )

    assert SkyWindow._water_overlay_status_line(dummy) == "W 3"


def test_hidden_water_status_line_shows_placeholder_icon() -> None:
    dummy = SimpleNamespace(
        water_overlay_state=SimpleNamespace(
            banner_text=None,
            points=None,
            current_source=None,
            current_mode=None,
        ),
        water_overlay_opacity=0.0,
    )

    assert SkyWindow._water_overlay_status_line(dummy) == "W ---"


def test_hidden_status_lines_show_placeholder_icons() -> None:
    dummy = SimpleNamespace(
        cloud_state=SimpleNamespace(banner_text="Clouds: downloading…"),
        cloud_disc_alpha=0.0,
        _predicted_cloud_satellite=lambda: "G19",
        satellite_state=SimpleNamespace(banner_text="Satellites: partial"),
        satellite_opacity=0.0,
        aircraft_state=SimpleNamespace(banner_text="Aircraft: unavailable"),
        aircraft_opacity=0.0,
        terrain_horizon_state=SimpleNamespace(banner_text="Terrain horizon: unavailable", current_source=None),
        terrain_horizon_opacity=0.0,
        urban_outline_state=SimpleNamespace(banner_text="Urban outline: unavailable", current_source=None),
        urban_outline_opacity=0.0,
    )

    assert SkyWindow._cloud_status_line(dummy) == "☁ ---"
    assert SkyWindow._satellite_status_line(dummy) == "🛰 ---"
    assert SkyWindow._aircraft_status_line(dummy) == "✈ ---"
    assert SkyWindow._terrain_horizon_status_line(dummy) == "▲ ---"
    assert SkyWindow._urban_outline_status_line(dummy) == "🂓 ---"


def test_default_cloud_horizon_cutoff_is_one_degree() -> None:
    assert DEFAULT_CLOUD_ALT_MIN_DEG == 1.0
