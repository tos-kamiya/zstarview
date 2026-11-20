from __future__ import annotations

from pathlib import Path

from zstarview.gui import launch_profile


def test_gui_launch_profile_roundtrip(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(launch_profile, "user_config_dir", lambda *_args: str(tmp_path))

    profile = {
        "city": "Tokyo",
        "enlarge_moon": True,
        "window_geometry": "restore",
        "theme": "day",
        "cloud_stripe": "cutout,8",
        "tropical_cyclone_opacity": 0.27,
    }

    launch_profile.save_gui_launch_profile(profile)
    loaded = launch_profile.load_gui_launch_profile()
    assert loaded["city"] == "Tokyo"
    assert loaded["window_geometry"] == "restore"
    assert loaded["theme"] == "day"
    assert loaded["cloud_stripe"] == "cutout,8"
    assert loaded["tropical_cyclone_opacity"] == 0.27
    assert "enlarge_moon" not in loaded

    structured_profile = {
        "city": {
            "resolver": "nominatim",
            "query": "Matsue Station",
            "result": {
                "name": "Matsue Station, Matsue, Shimane, Japan",
                "lat": 35.4641778,
                "lon": 133.0628539,
            },
        }
    }
    launch_profile.save_gui_launch_profile(structured_profile)
    structured_loaded = launch_profile.load_gui_launch_profile()
    assert structured_loaded["city"]["resolver"] == "nominatim"
    assert structured_loaded["city"]["query"] == "Matsue Station"
    assert structured_loaded["city"]["result"]["name"] == "Matsue Station, Matsue, Shimane, Japan"

    launch_profile.reset_gui_launch_profile()
    reset = launch_profile.load_gui_launch_profile()
    assert reset["theme"] == "night"
    assert reset["window_geometry"] is None

    profile_file = Path(tmp_path) / launch_profile.GUI_LAUNCH_PROFILE_FILENAME
    assert profile_file.exists()
    assert "enlarge_moon" not in profile_file.read_text(encoding="utf-8")


def test_gui_launch_profile_drops_removed_shell_cloud_settings(
    tmp_path, monkeypatch
) -> None:
    monkeypatch.setattr(launch_profile, "user_config_dir", lambda *_args: str(tmp_path))
    path = Path(tmp_path) / launch_profile.GUI_LAUNCH_PROFILE_FILENAME
    path.write_text(
        '{"city": "Tokyo", "cloud_stripe": "halftone2,30,1.7", '
        '"cloud_missing_tint_opacity": 0.2, "cloud_mode": "shell"}',
        encoding="utf-8",
    )
    assert launch_profile.load_gui_launch_profile() == {"city": "Tokyo"}

    path.write_text('{"cloud_stripe": ["width", 50, 0.85]}', encoding="utf-8")
    assert launch_profile.load_gui_launch_profile() == {}

    path.write_text('{"cloud_stripe": ["cutout", 24, 1.0]}', encoding="utf-8")
    assert launch_profile.load_gui_launch_profile() == {
        "cloud_stripe": ["cutout", 24, 1.0]
    }
