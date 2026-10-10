from __future__ import annotations

import subprocess
import sys
from types import SimpleNamespace

import pytest

from zstarview.clouddisc.workers import cloud_voxel_worker


def test_parent_transport_import_does_not_load_cloud_rendering_modules() -> None:
    completed = subprocess.run(
        [sys.executable, "-c", (
            "import sys; "
            "from zstarview.clouddisc.workers import cloud_voxel_worker; "
            "assert 'zstarview.render.cloud_voxels' not in sys.modules; "
            "assert 'zstarview.render.cloud_ground_light' not in sys.modules"
        )],
        capture_output=True, text=True, timeout=30,
    )
    assert completed.returncode == 0, completed.stderr


def test_render_child_crash_is_reported_without_in_process_fallback(monkeypatch) -> None:
    monkeypatch.setattr(cloud_voxel_worker.subprocess, "run", lambda *args, **kwargs:
                        SimpleNamespace(returncode=-11, stderr="child received SIGSEGV"))
    with pytest.raises(RuntimeError, match="without a result manifest.*SIGSEGV"):
        cloud_voxel_worker.render_cloud_voxels_in_subprocess(
            SimpleNamespace(source_key="crash-test"),
            lat=0, lon=0, view_center=(45, 0), edge_fov_deg=90, content_fov_deg=110,
            radius_px=1, sun_alt_deg=-20, sun_az_deg=0, sunlight_mix=0,
            night_color_rgb=(0.36, 0.385, 0.44), request_id=1,
        )
