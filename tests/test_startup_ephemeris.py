from __future__ import annotations

import pytest

from zstarview.astro import _starfield_load
from zstarview.launch_location_time import LaunchSetupError
from zstarview.gui.viewer import _verify_ephemeris_for_launch


def test_de442s_uses_naif_planets_url() -> None:
    from zstarview.paths import EPHEMERIS_URL

    assert _starfield_load.build_url("de442s.bsp") == EPHEMERIS_URL


def test_startup_verify_ephemeris_passes_when_loader_succeeds(monkeypatch) -> None:
    calls: list[str] = []

    def fake_loader(filename: str) -> object:
        calls.append(filename)
        return object()

    monkeypatch.setattr("zstarview.gui.viewer._starfield_load", fake_loader)
    monkeypatch.setattr("zstarview.gui.viewer.EPHEMERIS_FILENAME", "de442s.bsp")

    _verify_ephemeris_for_launch()

    assert calls == ["de442s.bsp"]


def test_startup_verify_ephemeris_aborts_on_oserror(monkeypatch) -> None:
    def fake_loader(_filename: str) -> object:
        raise OSError("network blocked")

    monkeypatch.setattr("zstarview.gui.viewer._starfield_load", fake_loader)

    with pytest.raises(LaunchSetupError):
        _verify_ephemeris_for_launch()
