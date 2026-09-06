"""Runtime diagnostics for reproducible crash and worker reports."""

from __future__ import annotations

import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import sysconfig
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from .__about__ import __version__

_DEPENDENCY_DISTRIBUTIONS = (
    "PySide6",
    "numpy",
    "astropy",
    "skyfield",
    "numba",
    "Pillow",
    "polars",
    "xarray",
    "h5py",
    "pyproj",
    "rasterio",
)


def _git_revision() -> str:
    configured = os.getenv("ZSTARVIEW_APP_REVISION", "").strip()
    if configured:
        return configured
    package_root = Path(__file__).resolve().parents[2]
    try:
        completed = subprocess.run(
            ("git", "rev-parse", "--short", "HEAD"),
            cwd=package_root,
            capture_output=True,
            check=True,
            text=True,
            timeout=0.5,
        )
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    revision = completed.stdout.strip()
    return revision or "unknown"


def _gil_state() -> str:
    probe = getattr(sys, "_is_gil_enabled", None)
    if callable(probe):
        try:
            return "enabled" if bool(probe()) else "disabled"
        except Exception:
            return "unknown"
    configured = sysconfig.get_config_var("Py_GIL_DISABLED")
    if configured in {0, "0", False, None}:
        return "enabled"
    if configured in {1, "1", True}:
        return "disabled"
    return "unknown"


def _dependency_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for distribution in _DEPENDENCY_DISTRIBUTIONS:
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[distribution] = "not-installed"
        except Exception:
            versions[distribution] = "unknown"
    return versions


def collect_runtime_diagnostics(
    *,
    session_id: str | None = None,
    worker_epoch: int | None = None,
    app_revision: str | None = None,
) -> dict[str, Any]:
    """Collect stable, side-effect-free process and dependency information."""
    diagnostics: dict[str, Any] = {
        "app_version": __version__,
        "app_revision": app_revision or _git_revision(),
        "pid": os.getpid(),
        "parent_pid": os.getppid(),
        "session_id": session_id or "unknown",
        "worker_epoch": worker_epoch if worker_epoch is not None else "unknown",
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "python_executable": os.path.abspath(sys.executable),
        "gil_state": _gil_state(),
        "os": platform.system(),
        "os_release": platform.release(),
        "machine": platform.machine(),
        "platform": platform.platform(aliased=True),
        "cwd": os.getcwd(),
        "dependencies": _dependency_versions(),
    }
    return json.loads(json.dumps(diagnostics, ensure_ascii=True, sort_keys=True))


def format_runtime_diagnostics(diagnostics: Mapping[str, Any]) -> str:
    """Format diagnostics as ASCII JSON suitable for logs and bug reports."""
    return json.dumps(diagnostics, ensure_ascii=True, sort_keys=True, indent=2)
