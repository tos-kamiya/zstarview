from __future__ import annotations

from datetime import datetime, timedelta, timezone
from types import SimpleNamespace

from zstarview.gui.aircraft_state import AircraftState
from zstarview.gui.water_overlay_state import WaterOverlayState
from zstarview.gui.window import DEFAULT_CLOUD_ALT_MIN_DEG, SkyWindow
from zstarview.gui.window_updates import SkyWindowUpdatesMixin
from zstarview.tropical_cyclones.cache import TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS


def _dummy_window(cloud_state):
    return SimpleNamespace(
        cloud_state=cloud_state,
        cloud_disc_alpha=0.2,
        _geo_satellite_mode_active=lambda: False,
        _predicted_cloud_satellite=lambda: "G19",
    )


def test_cloud_status_line_shows_downloading() -> None:
    state = SimpleNamespace(
        current_satellite="G19",
        banner_text="Clouds: downloading…",
        meta=None,
        coverage_ratio=None,
    )
    got = SkyWindow._cloud_status_line(_dummy_window(state))
    assert got == "☁ GOES G19 downloading"


def test_cloud_status_line_shows_download_failed() -> None:
    state = SimpleNamespace(
        current_satellite="G18",
        banner_text="Clouds: Clouds fetch timed out",
        meta=None,
        coverage_ratio=None,
    )
    got = SkyWindow._cloud_status_line(_dummy_window(state))
    assert got == "☁ GOES G18 failed"


def test_cloud_status_line_shows_partial_when_coverage_incomplete() -> None:
    meta = SimpleNamespace(
        satellite="G19",
        product="ISatSS-B13",
        time_utc=datetime(2026, 3, 5, 1, 20, tzinfo=timezone.utc),
    )
    state = SimpleNamespace(
        current_satellite=None,
        banner_text=None,
        meta=meta,
        coverage_ratio=0.72,
        source_completeness_ratio=1.0,
    )
    got = SkyWindow._cloud_status_line(_dummy_window(state))
    assert got == "☁ GOES G19 72% 01:20Z"


def test_cloud_status_line_shows_question_mark_for_partial_source() -> None:
    meta = SimpleNamespace(
        satellite="HIMAWARI",
        product="ISatSS-B13",
        time_utc=datetime(2026, 3, 5, 1, 20, tzinfo=timezone.utc),
    )
    state = SimpleNamespace(
        current_satellite=None,
        banner_text=None,
        meta=meta,
        coverage_ratio=1.0,
        source_expected_count=88,
        source_available_count=82,
        source_completeness_ratio=82.0 / 88.0,
    )
    got = SkyWindow._cloud_status_line(_dummy_window(state))
    assert got == "☁ HIMAWARI ? 01:20Z"


def test_cloud_status_line_shows_idle_without_meta_or_banner() -> None:
    state = SimpleNamespace(
        current_satellite=None,
        banner_text=None,
        meta=None,
        coverage_ratio=None,
    )
    got = SkyWindow._cloud_status_line(_dummy_window(state))
    assert got == "☁ GOES G19 idle"


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


def test_aircraft_status_line_compacts_stale_cache_banner() -> None:
    state = SimpleNamespace(
        banner_text="Aircraft: cache",
        last_success_utc=datetime(2026, 3, 5, 1, 20, tzinfo=timezone.utc),
    )
    dummy = SimpleNamespace(aircraft_state=state, aircraft_opacity=0.2)

    assert SkyWindow._aircraft_status_line(dummy) == "✈ cache"


def test_aircraft_rate_limited_skip_does_not_store_other_process_time() -> None:
    refreshed_at = datetime(2026, 3, 5, 1, 20, tzinfo=timezone.utc)
    calls: list[str] = []
    dummy = SimpleNamespace(
        aircraft_state=AircraftState(),
        aircraft_opacity=0.2,
        _schedule_next_aircraft_refresh=lambda *args, **kwargs: calls.append(
            "schedule"
        ),
        reproject_aircraft_overlay=lambda: calls.append("reproject"),
        request_client_update=lambda: calls.append("update"),
        _queue_aircraft_debug_snapshot=lambda payload: calls.append("debug"),
    )

    SkyWindow._on_aircraft_ready(
        dummy,
        {
            "snapshots": [],
            "bbox": None,
            "refreshed_at_utc": refreshed_at,
            "banner": "Aircraft: deferred",
            "source": "rate-limited-skip",
        },
    )

    assert dummy.aircraft_state.last_success_utc is None
    assert SkyWindow._aircraft_status_line(dummy) == "✈ deferred"


def test_terrain_and_urban_status_lines_show_icons() -> None:
    terrain_state = SimpleNamespace(
        banner_text="Terrain horizon: loading DEM...", current_source="Dem: cache"
    )
    urban_state = SimpleNamespace(
        banner_text=None,
        outlines=[object(), object(), object()],
        base_outline_count=2,
        skyscraper_outline_count=3,
        current_source="Urban: cache",
    )
    dummy = SimpleNamespace(
        terrain_horizon_state=terrain_state,
        terrain_horizon_opacity=0.2,
        urban_outline_state=urban_state,
        urban_outline_opacity=0.2,
    )

    assert SkyWindow._terrain_horizon_status_line(dummy) == "△ loading DEM..."
    assert SkyWindow._urban_outline_status_line(dummy) == "🂓 Overture Maps 2+3"


def test_terrain_status_line_compacts_stale_cache_source() -> None:
    dummy = SimpleNamespace(
        terrain_horizon_state=SimpleNamespace(
            banner_text=None,
            current_source="Dem: cache-stale",
        ),
        terrain_horizon_opacity=0.2,
    )

    assert SkyWindow._terrain_horizon_status_line(dummy) == "△ cache"


def test_urban_status_line_falls_back_to_merged_count_when_split_counts_missing() -> (
    None
):
    urban_state = SimpleNamespace(
        banner_text=None,
        outlines=[object(), object(), object()],
        base_outline_count=None,
        skyscraper_outline_count=None,
        current_source="Urban: cache",
    )
    dummy = SimpleNamespace(
        urban_outline_state=urban_state,
        urban_outline_opacity=0.2,
    )

    assert SkyWindow._urban_outline_status_line(dummy) == "🂓 Overture Maps 3"


def test_urban_status_line_shows_plateau_source() -> None:
    urban_state = SimpleNamespace(
        banner_text=None,
        outlines=[object(), object()],
        base_outline_count=2,
        skyscraper_outline_count=None,
        current_source="Urban: PLATEAU",
    )
    dummy = SimpleNamespace(
        urban_outline_state=urban_state,
        urban_outline_opacity=0.2,
    )

    assert SkyWindow._urban_outline_status_line(dummy) == "🂓 PLATEAU 2"


def test_tropical_cyclone_status_line_shows_no_entry_for_empty_collection() -> None:
    state = SimpleNamespace(
        banner_text=None,
        snapshots=(),
        snapshot_collection=SimpleNamespace(),
    )
    dummy = SimpleNamespace(
        show_tropical_cyclone_overlay=True,
        tropical_cyclone_state=state,
    )

    got = SkyWindow._tropical_cyclone_status_line(dummy)

    assert got == "TC none"


def test_tropical_cyclone_status_line_keeps_idle_before_first_result() -> None:
    state = SimpleNamespace(
        banner_text=None,
        snapshots=(),
        snapshot_collection=None,
    )
    dummy = SimpleNamespace(
        show_tropical_cyclone_overlay=True,
        tropical_cyclone_state=state,
    )

    got = SkyWindow._tropical_cyclone_status_line(dummy)

    assert got == "TC idle"


def test_tropical_cyclone_failure_clears_overlay_and_schedules_retry() -> None:
    now = datetime.now(timezone.utc)
    state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=object(),
        banner_text=None,
        cached_at_utc=now,
        last_checked_utc=now,
        next_check_utc=now,
        next_refresh_utc=now,
        projection_next_refresh_utc=now,
        source_url="https://example.invalid",
        current_source="test",
        set_error_banner=lambda text: (
            setattr(state, "snapshots", ()),
            setattr(state, "snapshot_collection", None),
            setattr(state, "cached_at_utc", None),
            setattr(state, "last_checked_utc", None),
            setattr(state, "next_check_utc", None),
            setattr(state, "next_refresh_utc", None),
            setattr(state, "projection_next_refresh_utc", None),
            setattr(state, "banner_text", text),
        ),
    )
    dummy = SimpleNamespace(tropical_cyclone_state=state, request_client_update=lambda: None)

    before = datetime.now(timezone.utc)
    SkyWindow._on_tropical_cyclone_failed(dummy, {"banner": "Typhoon: unavailable"})
    after = datetime.now(timezone.utc)

    assert state.snapshots == ()
    assert state.banner_text == "Typhoon: unavailable"
    assert state.next_check_utc is not None
    assert before + timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS) <= state.next_check_utc <= after + timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS)


def test_water_status_line_shows_only_count_when_enabled() -> None:
    dummy = SimpleNamespace(
        water_overlay_state=WaterOverlayState(
            sea_level_dots=[object(), object()],
            inland_dots=[object()],
            banner_text=None,
            current_source="Water: overpass",
            current_mode="sea",
        ),
        water_overlay_opacity=0.2,
    )

    assert SkyWindow._water_overlay_status_line(dummy) == "W 2+1"


def test_water_status_line_uses_question_mark_for_pending_counts() -> None:
    dummy = SimpleNamespace(
        water_overlay_state=WaterOverlayState(
            sea_level_dots=[object(), object()],
            inland_dots=None,
            banner_text=None,
            current_source="Water: overpass",
            current_mode="sea",
        ),
        water_overlay_opacity=0.2,
    )

    assert SkyWindow._water_overlay_status_line(dummy) == "W 2+?"


def test_water_status_line_identifies_cache_fallback() -> None:
    dummy = SimpleNamespace(
        water_overlay_state=WaterOverlayState(
            sea_level_dots=[object(), object()],
            inland_dots=[object()],
            banner_text=None,
            current_source="Water: cache-stale",
            current_mode="sea",
        ),
        water_overlay_opacity=0.2,
    )

    assert SkyWindow._water_overlay_status_line(dummy) == "W cache"


def test_water_status_line_compacts_unavailable_banner() -> None:
    dummy = SimpleNamespace(
        water_overlay_state=WaterOverlayState(
            banner_text="Water: unavailable",
            current_source=None,
        ),
        water_overlay_opacity=0.2,
    )

    assert SkyWindow._water_overlay_status_line(dummy) == "W unavailable"


def test_water_overlay_started_does_not_override_visible_points_banner() -> None:
    calls: list[str] = []
    dummy = SimpleNamespace(
        water_overlay_state=WaterOverlayState(
            dots=[object()],
            banner_text=None,
            current_source="Water: Overpass",
            current_mode="sea",
        ),
        request_client_update=lambda: calls.append("update"),
    )

    SkyWindowUpdatesMixin._on_water_overlay_started(
        dummy, {"banner": "Water: loading..."}
    )

    assert dummy.water_overlay_state.banner_text is None
    assert calls == ["update"]


def test_water_overlay_started_sets_banner_before_points_exist() -> None:
    calls: list[str] = []
    dummy = SimpleNamespace(
        water_overlay_state=WaterOverlayState(),
        request_client_update=lambda: calls.append("update"),
    )

    SkyWindowUpdatesMixin._on_water_overlay_started(
        dummy, {"banner": "Water: loading..."}
    )

    assert dummy.water_overlay_state.banner_text == "Water: loading..."
    assert calls == ["update"]


def test_hidden_water_status_line_shows_placeholder_icon() -> None:
    dummy = SimpleNamespace(
        water_overlay_state=SimpleNamespace(
            banner_text=None,
            dots=None,
            current_source=None,
            current_mode=None,
        ),
        water_overlay_opacity=0.0,
    )

    assert SkyWindow._water_overlay_status_line(dummy) == "W ---"


def test_road_night_lights_status_line_reports_cache_and_count() -> None:
    dummy = SimpleNamespace(
        road_night_lights_opacity=0.12,
        road_night_lights_status="cache 17",
    )

    assert SkyWindow._road_night_lights_status_line(dummy) == "R cache 17"


def test_road_night_lights_status_line_is_hidden_when_disabled() -> None:
    dummy = SimpleNamespace(road_night_lights_opacity=0.0)

    assert SkyWindow._road_night_lights_status_line(dummy) == "R ---"


def test_hidden_status_lines_show_placeholder_icons() -> None:
    dummy = SimpleNamespace(
        cloud_state=SimpleNamespace(banner_text="Clouds: downloading…"),
        cloud_disc_alpha=0.0,
        _predicted_cloud_satellite=lambda: "G19",
        satellite_state=SimpleNamespace(banner_text="Satellites: partial"),
        satellite_opacity=0.0,
        aircraft_state=SimpleNamespace(banner_text="Aircraft: unavailable"),
        aircraft_opacity=0.0,
        terrain_horizon_state=SimpleNamespace(
            banner_text="Terrain horizon: unavailable", current_source=None
        ),
        terrain_horizon_opacity=0.0,
        urban_outline_state=SimpleNamespace(
            banner_text="Urban outline: unavailable", current_source=None
        ),
        urban_outline_opacity=0.0,
    )

    assert SkyWindow._cloud_status_line(dummy) == "☁ ---"
    assert SkyWindow._satellite_status_line(dummy) == "🛰 ---"
    assert SkyWindow._aircraft_status_line(dummy) == "✈ ---"
    assert SkyWindow._terrain_horizon_status_line(dummy) == "△ ---"
    assert SkyWindow._urban_outline_status_line(dummy) == "🂓 ---"


def test_default_cloud_horizon_cutoff_is_one_degree() -> None:
    assert DEFAULT_CLOUD_ALT_MIN_DEG == 1.0
