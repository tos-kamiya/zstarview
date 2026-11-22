from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _block_ip_location_lookup(monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep tests offline: the default location now tries IP lookup first."""

    def _offline() -> dict[str, object]:
        raise RuntimeError("IP location lookup is disabled in tests")

    monkeypatch.setattr(
        "zstarview.location_resolver.resolve.fetch_location_by_ip", _offline
    )
