from __future__ import annotations

import pytest

import zstarview.gui.viewer as viewer_module
from zstarview.launch_time import LaunchSetupError


def _missing(_bin: str) -> str:
    raise FileNotFoundError("missing")


def test_missing_overturemaps_stops_launch(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    monkeypatch.setattr(viewer_module, "resolve_overturemaps_executable_path", _missing)

    with pytest.raises(LaunchSetupError):
        viewer_module._verify_overturemaps_for_launch(0.2)

    assert "overturemaps CLI not found" in caplog.text
    assert "--urban-outline-opacity 0" in caplog.text


def test_missing_overturemaps_is_ignored_when_urban_outlines_are_off(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(viewer_module, "resolve_overturemaps_executable_path", _missing)

    viewer_module._verify_overturemaps_for_launch(0.0)


def test_available_overturemaps_allows_launch(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        viewer_module,
        "resolve_overturemaps_executable_path",
        lambda _bin: "/usr/bin/overturemaps",
    )

    viewer_module._verify_overturemaps_for_launch(0.2)
