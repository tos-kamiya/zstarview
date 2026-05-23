from __future__ import annotations

from pathlib import Path

from zstarview.gui import launch_profile


def test_gui_launch_profile_roundtrip(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(launch_profile, "user_config_dir", lambda *_args: str(tmp_path))

    profile = {
        "city": "Tokyo",
        "window_geometry": "restore",
        "theme": "day",
        "cloud_stripe": "width,50,0.85",
    }

    launch_profile.save_gui_launch_profile(profile)
    loaded = launch_profile.load_gui_launch_profile()
    assert loaded["city"] == "Tokyo"
    assert loaded["window_geometry"] == "restore"
    assert loaded["theme"] == "day"
    assert loaded["cloud_stripe"] == "width,50,0.85"

    launch_profile.reset_gui_launch_profile()
    reset = launch_profile.load_gui_launch_profile()
    assert reset["theme"] == "night"
    assert reset["window_geometry"] is None

    profile_file = Path(tmp_path) / launch_profile.GUI_LAUNCH_PROFILE_FILENAME
    assert profile_file.exists()
