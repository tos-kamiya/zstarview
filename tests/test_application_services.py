from __future__ import annotations

import pytest

from zstarview.gui.application_services import ApplicationServices, EphemerisProvider


def test_ephemeris_provider_loads_once_per_instance() -> None:
    calls: list[object] = []
    value = object()
    provider = EphemerisProvider(loader=lambda: calls.append(value) or value)

    assert provider.load() is value
    assert provider.load() is value
    assert calls == [value]


def test_application_services_own_independent_native_locks() -> None:
    first = ApplicationServices(max_workers=1)
    second = ApplicationServices(max_workers=1)
    try:
        assert first.native_work_lock is not second.native_work_lock
    finally:
        first.shutdown()
        second.shutdown()


def test_application_services_reject_work_after_shutdown() -> None:
    services = ApplicationServices(max_workers=1)
    services.shutdown()

    with pytest.raises(RuntimeError, match="application services are shut down"):
        services.submit(lambda: None)
