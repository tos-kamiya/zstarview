from __future__ import annotations

import io
import sys

import pytest

from zstarview.astro import _starfield_load, load_ephemeris
from zstarview.launch_time import LaunchSetupError
from zstarview.gui.viewer import _verify_ephemeris_for_launch


def test_de442s_uses_naif_planets_url() -> None:
    from zstarview.paths import EPHEMERIS_URL

    assert _starfield_load.build_url("de442s.bsp") == EPHEMERIS_URL


def test_startup_verify_ephemeris_passes_when_loader_succeeds(monkeypatch) -> None:
    calls: list[str] = []

    def fake_load_ephemeris() -> object:
        calls.append("de442s.bsp")
        return object()

    monkeypatch.setattr("zstarview.gui.viewer.load_ephemeris", fake_load_ephemeris)

    _verify_ephemeris_for_launch()

    assert calls == ["de442s.bsp"]


def test_startup_verify_ephemeris_aborts_on_oserror(monkeypatch) -> None:
    def fake_load_ephemeris() -> object:
        raise OSError("network blocked")

    monkeypatch.setattr("zstarview.gui.viewer.load_ephemeris", fake_load_ephemeris)

    with pytest.raises(LaunchSetupError):
        _verify_ephemeris_for_launch()


def test_load_ephemeris_provides_dummy_standard_streams_when_missing(monkeypatch) -> None:
    captured: list[tuple[object | None, object | None, str]] = []

    def fake_loader(filename: str) -> object:
        captured.append((sys.stdout, sys.stderr, filename))
        return object()

    monkeypatch.setattr("zstarview.astro._starfield_load", fake_loader)
    monkeypatch.setattr("sys.stdout", None)
    monkeypatch.setattr("sys.stderr", None)

    load_ephemeris()

    assert len(captured) == 1
    stdout_obj, stderr_obj, filename = captured[0]
    assert isinstance(stdout_obj, io.StringIO)
    assert isinstance(stderr_obj, io.StringIO)
    assert filename == "de442s.bsp"
    assert sys.stdout is None
    assert sys.stderr is None


def test_load_ephemeris_caches_loaded_kernel(monkeypatch) -> None:
    calls: list[tuple[object | None, object | None, str]] = []
    sentinel = object()

    def fake_loader(filename: str) -> object:
        calls.append((sys.stdout, sys.stderr, filename))
        return sentinel

    monkeypatch.setattr("zstarview.astro._starfield_load", fake_loader)
    monkeypatch.setattr("zstarview.astro._cached_ephemeris", None)

    first = load_ephemeris()
    second = load_ephemeris()

    assert first is sentinel
    assert second is sentinel
    assert len(calls) == 1
