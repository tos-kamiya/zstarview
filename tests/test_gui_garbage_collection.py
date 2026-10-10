from __future__ import annotations

import gc
import subprocess
import sys
import threading
import weakref
from concurrent.futures import ThreadPoolExecutor

import pytest

from zstarview.gui import garbage_collection, viewer
from zstarview.gui.garbage_collection import GuiGarbageCollector


def test_worker_cycles_are_collected_on_owner_thread() -> None:
    owner = threading.get_ident()
    finalized = []
    refs = []

    class Cycle:
        def __init__(self):
            self.cycle = self

        def __del__(self):
            finalized.append(threading.get_ident())

    def allocate():
        for _ in range(1000):
            value = Cycle()
            refs.append(weakref.ref(value))

    was_enabled = gc.isenabled()
    gc.enable()
    try:
        with GuiGarbageCollector(threading.Lock()) as collector:
            with ThreadPoolExecutor(max_workers=1) as pool:
                pool.submit(allocate).result()
            assert not finalized
            collector._last_full_collection -= 31.0
            collector.collect_if_due()
            assert len(finalized) == 1000
            assert set(finalized) == {owner}
            assert all(ref() is None for ref in refs)
        assert gc.isenabled()
    finally:
        if not was_enabled:
            gc.disable()


def test_collection_skips_busy_native_work_without_losing_due_full_gc(monkeypatch) -> None:
    calls = []
    lock = threading.Lock()
    monkeypatch.setattr(gc, "collect", lambda generation: calls.append(generation))
    with GuiGarbageCollector(lock) as collector:
        collector._last_full_collection -= 31
        with lock:
            collector.collect_if_due()
            assert not calls
        collector.collect_if_due()
        assert calls == [2]


def test_young_and_full_collection_cadences(monkeypatch) -> None:
    calls = []
    lock = threading.Lock()
    clock = [0.0]
    monkeypatch.setattr(garbage_collection.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(gc, "get_threshold", lambda: (700, 10, 10))
    count = [0]
    monkeypatch.setattr(gc, "get_count", lambda: (count[0], 0, 0))

    def collect(generation):
        assert lock.locked()
        calls.append(generation)

    monkeypatch.setattr(gc, "collect", collect)
    with GuiGarbageCollector(lock) as collector:
        collector.collect_if_due()
        assert not calls
        count[0] = 701
        collector.collect_if_due()
        assert calls == [0]
        count[0] = 0
        clock[0] = 31
        collector.collect_if_due()
        assert calls == [0, 2]


@pytest.mark.parametrize("initial_enabled", [False, True])
def test_original_gc_state_is_restored_after_failure(monkeypatch, initial_enabled) -> None:
    original = gc.isenabled()
    (gc.enable if initial_enabled else gc.disable)()
    calls = []
    monkeypatch.setattr(gc, "collect", lambda generation: calls.append(generation))
    try:
        with pytest.raises(ValueError, match="startup failed"):
            with GuiGarbageCollector(threading.Lock()) as collector:
                assert not gc.isenabled()
                collector._last_full_collection -= 31
                collector.collect_if_due()
                # Even if external code toggles GC, restore the entry state.
                gc.enable()
                raise ValueError("startup failed")
        assert gc.isenabled() is initial_enabled
        assert calls == ([2] if initial_enabled else [])
    finally:
        (gc.enable if original else gc.disable)()


def test_collection_is_rejected_on_a_worker_thread() -> None:
    with GuiGarbageCollector(threading.Lock()) as collector:
        with ThreadPoolExecutor(max_workers=1) as pool:
            future = pool.submit(collector.collect_if_due)
            with pytest.raises(RuntimeError, match="owner thread"):
                future.result()


def test_viewer_drains_workers_before_restoring_gc(monkeypatch) -> None:
    events = []

    class Services:
        native_work_lock = threading.Lock()

        def shutdown(self, *, wait):
            assert wait
            assert not gc.isenabled()
            events.append("shutdown")

    def fail(**kwargs):
        assert not gc.isenabled()
        events.append("launch")
        raise ValueError("launch failed")

    original = gc.isenabled()
    monkeypatch.setattr(viewer, "ApplicationServices", Services)
    monkeypatch.setattr(viewer, "_main", fail)
    with pytest.raises(ValueError, match="launch failed"):
        viewer.main()
    assert events == ["launch", "shutdown"]
    assert gc.isenabled() is original


def test_real_qt_timer_collects_qobject_cycles_only_on_gui_thread() -> None:
    # Isolate real binding lifetimes from pytest's own QApplication/test state.
    code = '''
import gc, json, threading, time, weakref
from concurrent.futures import ThreadPoolExecutor
from PySide6.QtCore import QCoreApplication, QObject, QRectF, QTimer
from zstarview.gui.garbage_collection import GuiGarbageCollector
from zstarview.render.text import _rect_overlap_count
app = QCoreApplication([])
owner = threading.get_ident()
finalized, collections, refs = [], [], []
class QtCycle(QObject):
    def __init__(self):
        super().__init__()
        self.cycle = self
    def __del__(self):
        finalized.append(threading.get_ident())
def record(phase, info):
    if phase == 'start':
        collections.append(threading.get_ident())
def json_work():
    payload = json.dumps([{'samples': list(range(100))} for _ in range(50)])
    end = time.monotonic() + 2
    while time.monotonic() < end:
        json.loads(payload)
        value = []
        value.append(value)
def rect_work():
    rects = [QRectF(i, i, 20, 20) for i in range(30)]
    for rect in rects:
        _rect_overlap_count(rect, rects)
with GuiGarbageCollector(threading.Lock()) as collector:
    collector.start_timer(app)
    collector._last_full_collection -= 31
    gc.callbacks.append(record)
    for _ in range(100):
        value = QtCycle()
        refs.append(weakref.ref(value))
    del value
    timer = QTimer()
    timer.setInterval(5)
    timer.timeout.connect(rect_work)
    timer.start()
    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            first, second = pool.submit(json_work), pool.submit(json_work)
            QTimer.singleShot(2200, app.quit)
            app.exec()
            first.result()
            second.result()
        assert len(finalized) == 100, finalized
        assert set(finalized) == {owner}, finalized
        assert all(ref() is None for ref in refs)
        assert collections and set(collections) == {owner}, collections
    finally:
        timer.stop()
        gc.callbacks.remove(record)
assert gc.isenabled()
print('Qt timer/JSON/rectangle stress passed')
'''
    completed = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, timeout=20
    )
    assert completed.returncode == 0, completed.stderr
