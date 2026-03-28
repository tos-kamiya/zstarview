from __future__ import annotations

from types import SimpleNamespace

import zstarview.gui.viewer as viewer
from zstarview.location_resolver import LocationResolveError


class _DummyApp:
    def __init__(self) -> None:
        self.quit_on_last = []

    def setQuitOnLastWindowClosed(self, value: bool) -> None:
        self.quit_on_last.append(bool(value))


class _DummySplash:
    def __init__(self) -> None:
        self.closed = False

    def close(self) -> None:
        self.closed = True


class _DummyRootLogger:
    def __init__(self) -> None:
        self.removed = []

    def removeHandler(self, handler: object) -> None:
        self.removed.append(handler)


def test_main_keeps_splash_visible_for_location_error_then_cleans_up(monkeypatch) -> None:
    args = SimpleNamespace(
        city="Singapole",
        place=None,
        place_countrycode=None,
        place_lang="en",
        theme="night",
        clear_long_lived_cache=False,
    )
    app = _DummyApp()
    splash = _DummySplash()
    root_logger = _DummyRootLogger()
    splash_handler = object()
    sleep_calls: list[float] = []
    splash_context: list[str] = []

    monkeypatch.setattr(viewer, "parse_args", lambda: args)
    monkeypatch.setattr(viewer, "_handle_dataset_query_cli", lambda _args: None)
    monkeypatch.setattr(viewer, "setup_root_logger", lambda: root_logger)
    monkeypatch.setattr("zstarview.splash.setup_app", lambda _app_name: app)
    monkeypatch.setattr(
        "zstarview.splash.setup_splash_and_attach_logger",
        lambda _app, _app_name, _root_logger, _theme: (
            splash,
            splash_handler,
            lambda line: splash_context.append(line),
        ),
    )
    monkeypatch.setattr("zstarview.gui.window.SkyWindow", object)
    monkeypatch.setattr(
        viewer,
        "resolve_launch_location",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(LocationResolveError()),
    )
    monkeypatch.setattr(viewer.time, "sleep", lambda seconds: sleep_calls.append(float(seconds)))

    viewer.main()

    assert app.quit_on_last == [False]
    assert sleep_calls == [3.0]
    assert splash.closed is True
    assert root_logger.removed == [splash_handler]
    assert splash_context == []
