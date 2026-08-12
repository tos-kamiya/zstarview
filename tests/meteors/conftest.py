from __future__ import annotations

import pytest


def build_summary_row(
    *,
    trajectory_id: str = "20260810120000_TEST1",
    beginning_utc: str = "2026-08-10 12:00:00.000000",
    shower_code: str = "PER",
    begin_lat: float = 35.0,
    begin_lon: float = 139.0,
    begin_height_km: float = 105.0,
    end_lat: float = 35.1,
    end_lon: float = 139.1,
    end_height_km: float = 85.0,
    duration_s: float = 0.8,
    peak_abs_magnitude: float = -2.5,
) -> str:
    fields = [""] * 86
    fields[0] = trajectory_id
    fields[1] = "2461263.0"
    fields[2] = beginning_utc
    fields[3] = "7"
    fields[4] = shower_code
    fields[59] = "59.0"
    fields[63] = str(begin_lat)
    fields[65] = str(begin_lon)
    fields[67] = str(begin_height_km)
    fields[69] = str(end_lat)
    fields[71] = str(end_lon)
    fields[73] = str(end_height_km)
    fields[75] = str(duration_s)
    fields[76] = str(peak_abs_magnitude)
    return ";".join(fields)


@pytest.fixture
def summary_row_factory():
    return build_summary_row


@pytest.fixture
def sample_summary_text(summary_row_factory) -> str:
    return "\n".join(
        (
            "# Summary generated on 2026-08-12 05:38:52+00:00 UTC",
            "# synthetic fixture using the documented GMN column positions",
            summary_row_factory(),
        )
    )
