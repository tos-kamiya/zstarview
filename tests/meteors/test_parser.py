from __future__ import annotations

from datetime import date, datetime, timezone

from zstarview.meteors.parser import (
    parse_gmn_daily_index,
    parse_gmn_trajectory_summary,
)


def test_parse_daily_index_extracts_only_solrange_files() -> None:
    text = """
    <a href="traj_summary_20260809_solrange_137.0-138.0.txt">day</a>
    <a href="traj_summary_latest_daily.txt">latest</a>
    <a href="../">parent</a>
    <a href="traj_summary_20260810_solrange_138.0-139.0.txt">day</a>
    """

    files = parse_gmn_daily_index(text)

    assert [item.filename for item in files] == [
        "traj_summary_20260809_solrange_137.0-138.0.txt",
        "traj_summary_20260810_solrange_138.0-139.0.txt",
    ]
    assert files[0].nominal_date == date(2026, 8, 9)


def test_parse_trajectory_summary_uses_documented_columns(sample_summary_text) -> None:
    observations = parse_gmn_trajectory_summary(sample_summary_text)

    assert len(observations) == 1
    observation = observations[0]
    assert observation.trajectory_id == "20260810120000_TEST1"
    assert observation.beginning_utc == datetime(2026, 8, 10, 12, tzinfo=timezone.utc)
    assert observation.shower_code == "PER"
    assert observation.begin_lat_deg == 35.0
    assert observation.begin_lon_deg == 139.0
    assert observation.begin_height_km == 105.0
    assert observation.end_height_km == 85.0
    assert observation.duration_s == 0.8
    assert observation.peak_abs_magnitude == -2.5
    assert observation.initial_speed_km_s == 59.0


def test_parse_trajectory_summary_skips_malformed_rows(summary_row_factory) -> None:
    valid = summary_row_factory(trajectory_id="valid")
    invalid_coordinates = summary_row_factory(trajectory_id="bad", begin_lat=95.0)
    short = "too;short"

    observations = parse_gmn_trajectory_summary(
        "\n".join((short, invalid_coordinates, valid))
    )

    assert [item.trajectory_id for item in observations] == ["valid"]


def test_parse_trajectory_summary_normalizes_sporadic_code(summary_row_factory) -> None:
    observations = parse_gmn_trajectory_summary(
        summary_row_factory(shower_code="...")
    )

    assert observations[0].shower_code is None
