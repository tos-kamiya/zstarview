from __future__ import annotations

import datetime as dt

from zstarview.clouddisc.types import RenderKey, SourceKey, round_down_utc_to_slot


def test_round_down_utc_to_slot_rounds_to_10min() -> None:
    t = dt.datetime(2026, 3, 4, 12, 37, 59, tzinfo=dt.timezone.utc)
    got = round_down_utc_to_slot(t, slot_minutes=10)
    assert got == dt.datetime(2026, 3, 4, 12, 30, 0, tzinfo=dt.timezone.utc)


def test_round_down_utc_to_slot_accepts_naive_as_utc() -> None:
    t = dt.datetime(2026, 3, 4, 12, 37, 59)
    got = round_down_utc_to_slot(t, slot_minutes=10)
    assert got == dt.datetime(2026, 3, 4, 12, 30, 0, tzinfo=dt.timezone.utc)


def test_source_key_normalizes_timeslot_and_is_hashable() -> None:
    k1 = SourceKey(satellite="G19", timeslot_utc=dt.datetime(2026, 3, 4, 12, 39, tzinfo=dt.timezone.utc))
    k2 = SourceKey(satellite="G19", timeslot_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc))
    assert k1 == k2
    assert len({k1, k2}) == 1


def test_render_key_normalizes_az_alt_radius() -> None:
    source = SourceKey(satellite="HIMAWARI", timeslot_utc=dt.datetime(2026, 3, 4, 1, 2, tzinfo=dt.timezone.utc))
    rk = RenderKey(source=source, alt_deg=999.0, az_deg=-45.0, radius_px=0)
    assert rk.alt_deg == 89.999
    assert rk.az_deg == 315.0
    assert rk.radius_px == 1
