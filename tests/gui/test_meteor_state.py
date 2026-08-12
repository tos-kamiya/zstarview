from datetime import datetime, timezone
from types import SimpleNamespace

from zstarview.gui.meteor_state import MeteorState
from zstarview.gui.window_inputs import prepare_window_user_options
from zstarview.gui.window_updates import SkyWindowUpdatesMixin
from zstarview.meteors.types import MeteorWindowResult


def test_meteor_state_result_clears_banner() -> None:
    state = MeteorState(banner_text="GMN meteors: loading...")
    now = datetime(2026, 8, 12, tzinfo=timezone.utc)
    result = MeteorWindowResult(
        trails=(),
        window_start_utc=now,
        window_end_utc=now,
        source_files=(),
        unavailable_files=(),
    )
    state.set_result(result)
    assert state.result is result
    assert state.banner_text is None


def test_meteor_cli_zero_disables_gui_reenable() -> None:
    options = prepare_window_user_options(
        sky_disc_alpha=0.16,
        sky_disc_style="smooth",
        sky_disc_altaz_rings="dimalt",
        sky_disc_altaz_rings_hover="altaz",
        cloud_disc_alpha=0.2,
        satellite_opacity=0.7,
        aircraft_opacity=0.0,
        meteor_trails_opacity=0.0,
        terrain_horizon_opacity=0.003,
        earth_guide_opacity=0.028,
        urban_outline_opacity=0.2,
        ground_tint_opacity=0.025,
        enlarge_moon=False,
        bright_bodies_mode="outline",
        star_base_radius=4.0,
        vmag_limit=7.0,
        visual_preset="night",
        star_visibility_boost=1.0,
        visibility_boost=1.0,
        show_dso_initial=None,
        show_asterisms_initial=None,
        show_guidelines_initial=None,
        observation_info_mode=None,
        sky_disc_gui_allowed=True,
        cloud_gui_allowed=True,
        satellite_gui_allowed=True,
        aircraft_gui_allowed=False,
        meteor_trails_gui_allowed=False,
        tropical_cyclone_gui_allowed=True,
        terrain_horizon_gui_allowed=True,
        earth_guide_gui_allowed=True,
    )

    assert options.meteor_trails_opacity == 0.0
    assert options.meteor_trails_gui_allowed is False


def test_meteor_status_uses_zero_for_empty_result() -> None:
    now = datetime(2026, 8, 12, tzinfo=timezone.utc)
    state = MeteorState(
        result=MeteorWindowResult(
            trails=(),
            window_start_utc=now,
            window_end_utc=now,
            source_files=("sample.txt",),
            unavailable_files=(),
        )
    )
    dummy = SimpleNamespace(meteor_opacity=0.5, meteor_state=state)

    assert SkyWindowUpdatesMixin._meteor_status_line(dummy) == "M 0"  # type: ignore[arg-type]
