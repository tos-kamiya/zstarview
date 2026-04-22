from __future__ import annotations

from datetime import datetime, timezone

import pytest

from zstarview.search.jpl import resolve_jpl_target_altaz
from zstarview.search.jpl import resolve_jpl_target_state_vector
from zstarview.search.models import SearchJumpTarget
from zstarview.search.query import parse_search_query, search_target_matches_query
from zstarview.search.resolver import compute_search_target_altaz, resolve_search_targets


def test_parse_search_query_supports_label_selector() -> None:
    spec = parse_search_query("label=Ceres")

    assert spec.selector == "label"
    assert spec.value == "Ceres"
    assert spec.normalized == "ceres"


def test_parse_search_query_supports_id_selector() -> None:
    spec = parse_search_query("id=2000001")

    assert spec.selector == "id"
    assert spec.value == "2000001"
    assert spec.normalized == "2000001"


def test_search_target_matches_query_uses_label_and_id() -> None:
    target = SearchJumpTarget(
        label="Ceres",
        kind="jpl_small_body",
        sort_key=(0.0, "ceres"),
        object_key="2000001",
        command="DES=2000001;",
    )

    assert search_target_matches_query(target, parse_search_query("Ceres"))
    assert search_target_matches_query(target, parse_search_query("label=Ceres"))
    assert search_target_matches_query(target, parse_search_query("id=2000001"))


def test_resolve_search_targets_prefers_local_matches() -> None:
    local_targets = [
        SearchJumpTarget(
            label="Sirius",
            kind="star",
            sort_key=(0.0, "sirius"),
            subtitle="Vmag -1.44",
        ),
        SearchJumpTarget(
            label="Ceres",
            kind="jpl_small_body",
            sort_key=(0.0, "ceres"),
            object_key="2000001",
            command="DES=2000001;",
        ),
    ]
    resolution = resolve_search_targets(
        "Ceres",
        local_targets,
        jpl_search_callback=lambda _query: [],
    )

    assert len(resolution.candidates) == 1
    assert resolution.selected_target is not None
    assert resolution.selected_target.label == "Ceres"


def test_resolve_search_targets_prefers_known_satellites() -> None:
    local_targets = []
    satellite_calls: list[str] = []
    jpl_calls: list[str] = []

    def fake_satellite_search(query: str):
        satellite_calls.append(query)
        return [
            SearchJumpTarget(
                label="JWST",
                kind="satellite",
                sort_key=(0.0, "jwst"),
                object_key="JWST",
                alt_deg=10.0,
                az_deg=20.0,
            )
        ]

    def fake_jpl_search(query: str):
        jpl_calls.append(query)
        return []

    resolution = resolve_search_targets(
        "ISS",
        local_targets,
        satellite_search_callback=fake_satellite_search,
        jpl_search_callback=fake_jpl_search,
    )

    assert satellite_calls == ["ISS"]
    assert jpl_calls == []
    assert len(resolution.candidates) == 1
    assert resolution.selected_target is not None
    assert resolution.selected_target.label == "JWST"


def test_resolve_search_targets_skips_solar_system_bodies() -> None:
    local_targets = []
    satellite_calls: list[str] = []
    jpl_calls: list[str] = []

    def fake_satellite_search(query: str):
        satellite_calls.append(query)
        return []

    def fake_jpl_search(query: str):
        jpl_calls.append(query)
        return []

    resolution = resolve_search_targets(
        "Mars",
        local_targets,
        satellite_search_callback=fake_satellite_search,
        jpl_search_callback=fake_jpl_search,
    )

    assert resolution.candidates == ()
    assert satellite_calls == []
    assert jpl_calls == []


def test_resolve_search_targets_propagates_missing_satellite_position() -> None:
    def fake_satellite_search(_query: str):
        raise RuntimeError("Satellite position unavailable for JWST")

    with pytest.raises(RuntimeError, match="Satellite position unavailable"):
        resolve_search_targets(
            "ISS",
            [],
            satellite_search_callback=fake_satellite_search,
            jpl_search_callback=lambda _query: [],
        )


def test_compute_search_target_altaz_uses_satellite_altaz() -> None:
    target = SearchJumpTarget(
        label="JWST",
        kind="satellite",
        sort_key=(0.0, "jwst"),
        alt_deg=12.5,
        az_deg=278.0,
    )

    assert compute_search_target_altaz(
        target,
        observer_lat=35.0,
        observer_lon=135.0,
        observer_height_m=50.0,
    ) == (12.5, 278.0)


def test_compute_search_target_altaz_uses_satellite_resolver_when_missing_altaz() -> None:
    target = SearchJumpTarget(
        label="JWST",
        kind="satellite",
        sort_key=(0.0, "jwst"),
        object_key="JWST",
    )

    def fake_resolver(resolved_target: SearchJumpTarget) -> tuple[float, float] | None:
        assert resolved_target.label == "JWST"
        return 13.5, 279.0

    assert compute_search_target_altaz(
        target,
        observer_lat=35.0,
        observer_lon=135.0,
        observer_height_m=50.0,
        satellite_altaz_resolver=fake_resolver,
    ) == (13.5, 279.0)


def test_compute_search_target_altaz_uses_jpl_state_vector_when_present() -> None:
    target = SearchJumpTarget(
        label="Ceres",
        kind="jpl_small_body",
        sort_key=(0.0, "ceres"),
        object_key="2000001",
        command="DES=2000001;",
        horizons_epoch_utc=datetime(2026, 4, 21, 21, 15, 10, tzinfo=timezone.utc),
        horizons_position_km=(1.0, 2.0, 3.0),
        horizons_velocity_km_s=(0.1, 0.2, 0.3),
    )

    def fake_project(
        resolved_target: SearchJumpTarget,
        **kwargs,
    ) -> tuple[float, float] | None:
        assert resolved_target.label == "Ceres"
        assert kwargs["observer_lat"] == 35.0
        assert kwargs["observer_lon"] == 135.0
        assert kwargs["observer_height_m"] == 50.0
        return 21.5, 181.0

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            "zstarview.search.resolver.project_jpl_target_altaz_from_state_vector",
            fake_project,
        )
        assert compute_search_target_altaz(
            target,
            observer_lat=35.0,
            observer_lon=135.0,
            observer_height_m=50.0,
        ) == (21.5, 181.0)


def test_resolve_jpl_target_altaz_logs_target_time_and_observer(caplog) -> None:
    target = SearchJumpTarget(
        label="Voyager 1 (spacecraft)",
        kind="jpl_body",
        sort_key=(0.0, "voyager 1 (spacecraft)"),
        command="-31",
        target_time_utc=datetime(2026, 4, 21, 21, 15, 10, tzinfo=timezone.utc),
        jpl_group="mb",
    )

    def fake_observer_fetch(command: str, **kwargs):
        assert command == "-31"
        assert kwargs["target_time_utc"] == datetime(2026, 4, 21, 21, 15, 10, tzinfo=timezone.utc)
        assert kwargs["observer_lat"] == 35.28
        assert kwargs["observer_lon"] == 133.03
        assert kwargs["observer_height_m"] == 0.0
        return [["2026-Apr-21 21:15:10", "*", "m", "233.1", "42.1"]]

    with caplog.at_level("INFO"):
        altaz = resolve_jpl_target_altaz(
            target,
            observer_lat=35.28,
            observer_lon=133.03,
            observer_height_m=0.0,
            observer_fetch=fake_observer_fetch,
        )

    assert altaz == (42.1, 233.1)
    assert "Resolving JPL target alt/az: label=Voyager 1 (spacecraft)" in caplog.text
    assert "target_time_utc=2026-04-21T21:15:10+00:00" in caplog.text
    assert "observer=(lat=35.28 lon=133.03 height_m=0.0)" in caplog.text
    assert "Resolved JPL target alt/az: label=Voyager 1 (spacecraft) command=-31 alt=42.1 az=233.1" in caplog.text


def test_resolve_jpl_target_state_vector_logs_command_and_time(caplog) -> None:
    target = SearchJumpTarget(
        label="Voyager 1 (spacecraft)",
        kind="jpl_body",
        sort_key=(0.0, "voyager 1 (spacecraft)"),
        command="-31",
        target_time_utc=datetime(2026, 4, 21, 21, 15, 10, tzinfo=timezone.utc),
        jpl_group="mb",
    )

    def fake_vector_fetch(command: str, **kwargs):
        assert command == "-31"
        assert kwargs["target_time_utc"] == datetime(2026, 4, 21, 21, 15, 10, tzinfo=timezone.utc)
        return [["2026-Apr-21 21:15:10", "1.0", "2.0", "3.0", "0.1", "0.2", "0.3"]]

    with caplog.at_level("INFO"):
        state_vector = resolve_jpl_target_state_vector(
            target,
            vector_fetch=fake_vector_fetch,
        )

    assert state_vector == (
        datetime(2026, 4, 21, 21, 15, 10, tzinfo=timezone.utc),
        (1.0, 2.0, 3.0),
        (0.1, 0.2, 0.3),
    )
    assert "Resolving JPL target state vector: label=Voyager 1 (spacecraft)" in caplog.text
    assert "target_time_utc=2026-04-21T21:15:10+00:00" in caplog.text
    assert "Resolved JPL target state vector: label=Voyager 1 (spacecraft) command=-31 x=1.000 y=2.000 z=3.000" in caplog.text


def test_extract_horizons_state_vector_skips_julian_date_prefix() -> None:
    from zstarview.search.jpl import extract_horizons_state_vector

    rows = [
        [
            "2460792.500000000",
            "2026-Apr-21 00:00:00.0000",
            "1.0",
            "2.0",
            "3.0",
            "0.1",
            "0.2",
            "0.3",
            "0.4",
            "0.5",
            "0.6",
        ]
    ]

    assert extract_horizons_state_vector(rows) == (
        (1.0, 2.0, 3.0),
        (0.1, 0.2, 0.3),
    )
