"""Offline mode: start without any network access.

Offline mode is a process-wide switch carried by the ``ZSTARVIEW_OFFLINE``
environment variable, so worker subprocesses started by the viewer inherit it.

When it is on:

- A socket guard rejects every non-loopback connection with ``OfflineError``.
  ``OfflineError`` is an ``OSError``, so existing "use the cached copy when the
  download fails" paths keep working without changes.
- Astropy is told not to download IERS or leap-second tables.
- Layers built from short-lived data (aircraft, satellites, forecast
  precipitation, satellite clouds, tropical cyclones, meteor trails) are
  turned off, and the Moon/Sun hover images are not fetched.

Long-lived data (terrain, buildings, night lights, road lights, water) is used
only when it is already in the cache.
"""

from __future__ import annotations

import argparse
import ipaddress
import logging
import os
import socket
from typing import Any

logger = logging.getLogger(__name__)

OFFLINE_ENV_VAR = "ZSTARVIEW_OFFLINE"

# Layers whose data goes stale within hours; offline mode does not draw them.
SHORT_LIVED_OPACITY_ARGS = (
    "aircraft_opacity",
    "satellite_opacity",
    "precipitation_opacity",
    "cloud_opacity",
    "tropical_cyclone_opacity",
    "meteor_trails_opacity",
)


class OfflineError(OSError):
    """Raised when code tries to use the network in offline mode."""

    def __init__(self, target: object = None) -> None:
        detail = f" ({target})" if target is not None else ""
        super().__init__(f"Network access is disabled in offline mode{detail}")


class OfflineCacheMissError(OfflineError):
    """Raised when offline mode needs cached data that is not present."""

    def __init__(self, what: str) -> None:
        OSError.__init__(self, f"No cached {what} is available offline")


def is_offline_failure(exc: BaseException | None) -> bool:
    """Return True when ``exc`` was caused by offline mode blocking the network.

    The exception chain (``__cause__``, ``__context__``, and ``URLError.reason``)
    is searched, because fetch helpers often wrap the original error.
    """
    if exc is None or not is_offline():
        return False
    stack: list[BaseException] = [exc]
    seen: set[int] = set()
    while stack:
        current = stack.pop()
        if id(current) in seen:
            continue
        seen.add(id(current))
        if isinstance(current, OfflineError):
            return True
        for linked in (
            current.__cause__,
            current.__context__,
            getattr(current, "reason", None),
        ):
            if isinstance(linked, BaseException):
                stack.append(linked)
    return False


def failure_log_level(exc: BaseException | None, default: int = logging.WARNING) -> int:
    """Log offline-caused failures at INFO; others at ``default``."""
    return logging.INFO if is_offline_failure(exc) else default


def is_offline() -> bool:
    """Return True when offline mode is on for this process."""
    return os.environ.get(OFFLINE_ENV_VAR, "").strip() == "1"


def ensure_online(target: object = None) -> None:
    """Raise OfflineError when offline mode is on."""
    if is_offline():
        raise OfflineError(target)


_LOOPBACK_NAMES = frozenset({"localhost", "localhost.localdomain", "ip6-localhost"})


def _is_loopback_host(host: object) -> bool:
    if isinstance(host, bytes):
        try:
            host = host.decode("ascii")
        except UnicodeDecodeError:
            return False
    if not isinstance(host, str):
        return False
    name = host.strip().strip("[]").lower()
    if not name or name in _LOOPBACK_NAMES:
        return True
    try:
        return ipaddress.ip_address(name.split("%", 1)[0]).is_loopback
    except ValueError:
        return False


def _address_host(address: Any) -> object:
    if isinstance(address, tuple) and address:
        return address[0]
    return None


_original_socket_functions: tuple[Any, Any, Any] | None = None


def install_network_guard() -> None:
    """Make every non-loopback socket connection fail with OfflineError."""
    global _original_socket_functions
    if _original_socket_functions is not None:
        return
    original_getaddrinfo = socket.getaddrinfo
    original_connect = socket.socket.connect
    original_connect_ex = socket.socket.connect_ex
    _original_socket_functions = (
        original_getaddrinfo,
        original_connect,
        original_connect_ex,
    )

    def guarded_getaddrinfo(host: Any, *args: Any, **kwargs: Any) -> Any:
        if host is not None and not _is_loopback_host(host):
            raise OfflineError(host)
        return original_getaddrinfo(host, *args, **kwargs)

    def _check(sock: socket.socket, address: Any) -> None:
        if sock.family not in (socket.AF_INET, socket.AF_INET6):
            return
        host = _address_host(address)
        if not _is_loopback_host(host):
            raise OfflineError(host)

    def guarded_connect(self: socket.socket, address: Any) -> None:
        _check(self, address)
        original_connect(self, address)

    def guarded_connect_ex(self: socket.socket, address: Any) -> int:
        _check(self, address)
        return original_connect_ex(self, address)

    socket.getaddrinfo = guarded_getaddrinfo  # type: ignore[assignment]
    socket.socket.connect = guarded_connect  # type: ignore[method-assign, assignment]
    socket.socket.connect_ex = guarded_connect_ex  # type: ignore[method-assign, assignment]


def uninstall_network_guard() -> None:
    """Restore the original socket functions (used by tests)."""
    global _original_socket_functions
    if _original_socket_functions is None:
        return
    getaddrinfo, connect, connect_ex = _original_socket_functions
    socket.getaddrinfo = getaddrinfo  # type: ignore[assignment]
    socket.socket.connect = connect  # type: ignore[method-assign]
    socket.socket.connect_ex = connect_ex  # type: ignore[method-assign]
    _original_socket_functions = None


def _disable_astropy_downloads() -> None:
    try:
        from astropy.utils import iers
        from astropy.utils.data import conf as data_conf
    except Exception:  # pragma: no cover - astropy is a hard dependency
        return
    iers.conf.auto_download = False
    # Without downloads, the bundled IERS-A predictions age past the default
    # 30-day limit and astropy raises. The UT1 error from stale predictions is
    # far below what the viewer can show, so accept them.
    iers.conf.auto_max_age = None
    iers.conf.iers_degraded_accuracy = "warn"
    data_conf.allow_internet = False


def install_from_env() -> None:
    """Apply the offline guards when the environment requests offline mode."""
    if not is_offline():
        return
    install_network_guard()
    _disable_astropy_downloads()


def enable_offline_mode() -> None:
    """Turn offline mode on for this process and its child processes."""
    os.environ[OFFLINE_ENV_VAR] = "1"
    install_from_env()
    logger.info("Offline mode: network access is disabled")


def apply_offline_overrides(args: argparse.Namespace) -> list[str]:
    """Turn off layers that need fresh network data. Return the disabled names."""
    disabled: list[str] = []
    for name in SHORT_LIVED_OPACITY_ARGS:
        value = getattr(args, name, None)
        if value is None:
            continue
        try:
            enabled = float(value) > 0.0
        except (TypeError, ValueError):
            enabled = False
        if enabled:
            disabled.append(name)
        setattr(args, name, 0.0)
    if getattr(args, "moon_style", None) == "image":
        args.moon_style = "sphere"
        disabled.append("moon_style=image")
    if disabled:
        logger.info(
            "Offline mode: turned off layers that need fresh data: %s",
            ", ".join(disabled),
        )
    return disabled


def add_offline_argument(parser: argparse._ActionsContainer) -> None:
    parser.add_argument(
        "--offline",
        action="store_true",
        default=False,
        help=(
            "Start without network access. Cached long-lived data (terrain, "
            "buildings, night lights, road lights, water) is used when present; "
            "aircraft, satellites, forecast precipitation, satellite clouds, "
            "tropical cyclones, meteor trails, and Moon/Sun hover images are "
            "not used."
        ),
    )
