from __future__ import annotations

import threading
import time

from zstarview.gui.worker_pool import shutdown_gui_worker_pool, submit_gui_work


def test_shared_gui_worker_pool_serializes_tasks() -> None:
    first_started = threading.Event()
    second_started = threading.Event()
    release_first = threading.Event()

    def first_task() -> None:
        first_started.set()
        release_first.wait(timeout=1.0)

    def second_task() -> None:
        second_started.set()

    first_future = submit_gui_work(first_task)
    assert first_started.wait(timeout=1.0)

    second_future = submit_gui_work(second_task)
    time.sleep(0.05)
    assert not second_started.is_set()

    release_first.set()
    assert first_future.result(timeout=1.0) is None
    assert second_started.wait(timeout=1.0)
    assert second_future.result(timeout=1.0) is None

    shutdown_gui_worker_pool(wait=True)

    third_started = threading.Event()

    def third_task() -> None:
        third_started.set()

    third_future = submit_gui_work(third_task)
    assert third_started.wait(timeout=1.0)
    assert third_future.result(timeout=1.0) is None
