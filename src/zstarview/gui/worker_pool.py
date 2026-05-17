from __future__ import annotations

from concurrent.futures import Future, ThreadPoolExecutor, wait as wait_for_futures
import threading
from typing import Iterable

_GUI_WORKER_POOL_LOCK = threading.Lock()
_GUI_WORKER_POOL: ThreadPoolExecutor | None = None


def _get_gui_worker_pool() -> ThreadPoolExecutor:
    global _GUI_WORKER_POOL
    with _GUI_WORKER_POOL_LOCK:
        if _GUI_WORKER_POOL is None:
            # Shared GUI worker pool for all long-running background work.
            #
            # The pool is intentionally size 1 so sky/cloud/terrain/search/bootstrap
            # work cannot overlap in native-heavy code paths.
            _GUI_WORKER_POOL = ThreadPoolExecutor(
                max_workers=1,
                thread_name_prefix="zstarview-gui",
            )
        return _GUI_WORKER_POOL


def submit_gui_work(target, /, *args, **kwargs) -> Future:
    return _get_gui_worker_pool().submit(target, *args, **kwargs)


def wait_for_gui_futures(futures: Iterable[Future], timeout_s: float | None) -> None:
    futures = tuple(futures)
    if not futures:
        return
    wait_for_futures(futures, timeout=timeout_s)


def shutdown_gui_worker_pool(*, wait: bool = True) -> None:
    global _GUI_WORKER_POOL
    with _GUI_WORKER_POOL_LOCK:
        pool = _GUI_WORKER_POOL
        _GUI_WORKER_POOL = None
    if pool is not None:
        pool.shutdown(wait=wait, cancel_futures=False)
