from zstarview.gui.display_mode import (
    DISPLAY_MODE_INVERTED_CITY,
    DISPLAY_MODE_NORMAL,
    DISPLAY_MODE_SIMPLE_LABELS,
    DISPLAY_MODE_SIMPLE_NO_LABELS,
    display_mode_cycle,
    next_display_mode,
)


def test_display_mode_cycle_includes_inverted_city_when_outline_is_available() -> None:
    assert display_mode_cycle(urban_outline_available=True) == (
        DISPLAY_MODE_NORMAL,
        DISPLAY_MODE_INVERTED_CITY,
        DISPLAY_MODE_SIMPLE_NO_LABELS,
        DISPLAY_MODE_SIMPLE_LABELS,
    )


def test_display_mode_cycle_skips_inverted_city_without_outline() -> None:
    assert next_display_mode(
        DISPLAY_MODE_NORMAL, urban_outline_available=False
    ) == DISPLAY_MODE_SIMPLE_NO_LABELS


def test_display_mode_cycle_wraps_to_normal() -> None:
    assert next_display_mode(
        DISPLAY_MODE_SIMPLE_LABELS, urban_outline_available=True
    ) == DISPLAY_MODE_NORMAL
