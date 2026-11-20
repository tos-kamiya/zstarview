#!/usr/bin/env python3
"""Run the GUI with counters on the star renderer and the GUI image caches.

This is an observation helper for real interaction (idle, hover, drag, zoom,
sky updates).  It wraps, without changing behavior:

- ``SkyWindowRenderCacheMixin._render_cached_image``: hit/miss per cache slot
  (frame, present frame, star surface, fast frame, ...).
- ``render.stars._draw_stars_render``: calls classified as low-level image
  cache ``hit``, ``regen`` (cache entry replaced), or ``uncached`` (no cache
  or explicit ``screen_positions``), with wall time per class.

Counters are printed every ``--trace-interval`` seconds and at exit.  All
other arguments are passed to the normal ``zstarview`` command line.

Example:
  uv run -p .venv/bin/python dev-samples/trace_star_render_cache.py \
      --trace-interval 5 [zstarview options] [city]
  QT_QPA_PLATFORM=offscreen timeout 40 uv run -p .venv/bin/python \
      dev-samples/trace_star_render_cache.py --trace-exit-after 30
"""

from __future__ import annotations

import atexit
import collections
import sys
import threading
import time

from zstarview.gui import window_render_cache
from zstarview.render import stars as render_stars

_LOCK = threading.Lock()
_GUI_CACHE: collections.Counter[str] = collections.Counter()
_STAR_CALLS: collections.Counter[str] = collections.Counter()
_STAR_TIME_NS: collections.Counter[str] = collections.Counter()
_STARTED = time.monotonic()


def _install() -> None:
    mixin = window_render_cache.SkyWindowRenderCacheMixin
    original_cached = mixin._render_cached_image

    def traced_cached(self, *, image_size, frame_key, render_fn, cache_key_attr, cache_image_attr):  # type: ignore[no-untyped-def]
        before = self.__dict__.get(cache_image_attr)
        image = original_cached(
            self,
            image_size=image_size,
            frame_key=frame_key,
            render_fn=render_fn,
            cache_key_attr=cache_key_attr,
            cache_image_attr=cache_image_attr,
        )
        outcome = "hit" if image is before else "miss"
        with _LOCK:
            _GUI_CACHE[f"{cache_key_attr}:{outcome}"] += 1
        return image

    mixin._render_cached_image = traced_cached  # type: ignore[method-assign]

    original_render = render_stars._draw_stars_render

    def traced_render(painter, geometry, celestial_data, viewer_data, star_base_radius, **kwargs):  # type: ignore[no-untyped-def]
        cache = kwargs.get("render_cache")
        before = None if cache is None else cache.entry
        start = time.perf_counter_ns()
        original_render(painter, geometry, celestial_data, viewer_data, star_base_radius, **kwargs)
        elapsed = time.perf_counter_ns() - start
        if cache is None or kwargs.get("screen_positions") is not None:
            kind = "uncached"
        elif cache.entry is not before:
            kind = "regen"
        else:
            kind = "hit"
        variant = []
        if kwargs.get("fast_mode"):
            variant.append("fast")
        if kwargs.get("light_background_outline"):
            variant.append("light")
        if kwargs.get("draw_vmag_min_exclusive") is not None:
            variant.append("faint")
        label = f"{kind}[{','.join(variant) or 'normal'}]"
        with _LOCK:
            _STAR_CALLS[label] += 1
            _STAR_TIME_NS[label] += elapsed

    render_stars._draw_stars_render = traced_render  # type: ignore[assignment]


_LAST: dict[str, collections.Counter[str]] = {
    "gui": collections.Counter(),
    "calls": collections.Counter(),
    "time": collections.Counter(),
}


def _report(final: bool = False) -> None:
    """Print counters accumulated since the previous report (or totals at exit)."""
    with _LOCK:
        gui = collections.Counter(_GUI_CACHE)
        calls = collections.Counter(_STAR_CALLS)
        times = collections.Counter(_STAR_TIME_NS)
    if final:
        d_gui, d_calls, d_time = gui, calls, times
        title = "total"
    else:
        d_gui = gui - _LAST["gui"]
        d_calls = calls - _LAST["calls"]
        d_time = times - _LAST["time"]
        _LAST.update(gui=gui, calls=calls, time=times)
        title = "delta"
    stamp = time.strftime("%H:%M:%S")
    print(
        f"[star-trace] {stamp} t={time.monotonic() - _STARTED:.1f}s {title} gui-cache={dict(sorted(d_gui.items()))}",
        file=sys.stderr,
        flush=True,
    )
    for label, count in sorted(d_calls.items()):
        mean_ms = d_time[label] / count / 1e6
        print(
            f"[star-trace]   {label}: calls={count} mean={mean_ms:.2f}ms total={d_time[label] / 1e6:.0f}ms",
            file=sys.stderr,
            flush=True,
        )


def _pop_float_option(argv: list[str], name: str) -> float | None:
    if name in argv:
        index = argv.index(name)
        value = float(argv[index + 1])
        del argv[index : index + 2]
        return value
    return None


def main() -> None:
    argv = sys.argv[1:]
    interval = _pop_float_option(argv, "--trace-interval")
    exit_after = _pop_float_option(argv, "--trace-exit-after")
    sys.argv = [sys.argv[0], *argv]
    _install()
    atexit.register(_report, True)

    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QApplication

    def _arm_timers() -> None:
        app = QApplication.instance()
        if app is None:
            QTimer.singleShot(200, _arm_timers)
            return
        if interval:
            timer = QTimer(app)
            timer.timeout.connect(_report)
            timer.start(int(interval * 1000))
        if exit_after:
            QTimer.singleShot(int(exit_after * 1000), app.quit)

    # QTimer.singleShot needs an application; poll from a thread-free hook by
    # patching QApplication construction time via a deferred import hook.
    original_init = QApplication.__init__

    def patched_init(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        original_init(self, *args, **kwargs)
        QTimer.singleShot(0, _arm_timers)

    QApplication.__init__ = patched_init  # type: ignore[method-assign]

    from zstarview.gui.viewer import main as viewer_main

    viewer_main()


if __name__ == "__main__":
    main()
