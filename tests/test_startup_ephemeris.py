from __future__ import annotations

import pytest

from zstarview.astro import _starfield_load
from zstarview.startup import StartupAbortError, _startup_verify_ephemeris


def test_de442s_uses_naif_planets_url() -> None:
    from zstarview.paths import EPHEMERIS_URL

    assert _starfield_load.build_url("de442s.bsp") == EPHEMERIS_URL


def test_startup_verify_ephemeris_passes_when_loader_succeeds(monkeypatch) -> None:
    calls: list[str] = []

    def fake_loader(filename: str) -> object:
        calls.append(filename)
        return object()

    monkeypatch.setattr("zstarview.startup._starfield_load", fake_loader)
    monkeypatch.setattr("zstarview.startup.EPHEMERIS_FILENAME", "de442s.bsp")

    _startup_verify_ephemeris()

    assert calls == ["de442s.bsp"]


def test_startup_verify_ephemeris_aborts_on_oserror(monkeypatch) -> None:
    def fake_loader(_filename: str) -> object:
        raise OSError("network blocked")

    monkeypatch.setattr("zstarview.startup._starfield_load", fake_loader)

    with pytest.raises(StartupAbortError):
        _startup_verify_ephemeris()
