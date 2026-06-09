from __future__ import annotations

from zstarview.__about__ import __version__
from zstarview.user_agent import APP_USER_AGENT, build_user_agent


def test_build_user_agent_uses_app_version() -> None:
    assert APP_USER_AGENT == f"zstarview/{__version__}"


def test_build_user_agent_appends_service_suffix() -> None:
    assert build_user_agent("water-overlay") == f"zstarview/{__version__} (+water-overlay)"
