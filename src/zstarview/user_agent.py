from __future__ import annotations

from .__about__ import __version__

APP_USER_AGENT = f"zstarview/{__version__}"


def build_user_agent(service_suffix: str | None = None) -> str:
    """Build a zstarview User-Agent with an optional service suffix."""
    if service_suffix is None:
        return APP_USER_AGENT
    suffix = service_suffix.strip()
    if not suffix:
        return APP_USER_AGENT
    return f"{APP_USER_AGENT} (+{suffix})"
