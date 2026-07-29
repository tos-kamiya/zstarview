from zstarview.render.zstarview_pipeline import TIME_OF_DAY_MARKER_SKY_ALT_DEG


def test_time_of_day_marker_samples_sky_at_horizon() -> None:
    assert TIME_OF_DAY_MARKER_SKY_ALT_DEG == 0.0
