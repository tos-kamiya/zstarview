from __future__ import annotations

import threading

import pytest

from zstarview.gui.application_services import ApplicationServices


def test_shared_gui_worker_pool_allows_two_concurrent_tasks() -> None:
    services = ApplicationServices()
    first_started = threading.Event()
    second_started = threading.Event()
    release_first = threading.Event()

    def first_task() -> None:
        first_started.set()
        release_first.wait(timeout=1.0)

    def second_task() -> None:
        second_started.set()

    first_future = services.submit(first_task)
    assert first_started.wait(timeout=1.0)

    second_future = services.submit(second_task)
    assert second_started.wait(timeout=1.0)

    release_first.set()
    assert first_future.result(timeout=1.0) is None
    assert second_future.result(timeout=1.0) is None

    services.shutdown(wait=True)

    with pytest.raises(RuntimeError, match="application services are shut down"):
        services.submit(lambda: None)
