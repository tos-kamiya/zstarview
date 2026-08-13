"""Display modes cycled by the viewer's Space shortcut."""

from __future__ import annotations

DISPLAY_MODE_NORMAL = "normal"
DISPLAY_MODE_INVERTED_CITY = "inverted-city"
DISPLAY_MODE_SIMPLE_NO_LABELS = "simple-no-labels"
DISPLAY_MODE_SIMPLE_LABELS = "simple-labels"


def display_mode_cycle(*, urban_outline_available: bool) -> tuple[str, ...]:
    if urban_outline_available:
        return (
            DISPLAY_MODE_NORMAL,
            DISPLAY_MODE_INVERTED_CITY,
            DISPLAY_MODE_SIMPLE_NO_LABELS,
            DISPLAY_MODE_SIMPLE_LABELS,
        )
    return (
        DISPLAY_MODE_NORMAL,
        DISPLAY_MODE_SIMPLE_NO_LABELS,
        DISPLAY_MODE_SIMPLE_LABELS,
    )


def next_display_mode(current: str, *, urban_outline_available: bool) -> str:
    modes = display_mode_cycle(urban_outline_available=urban_outline_available)
    try:
        index = modes.index(current)
    except ValueError:
        index = 0
    return modes[(index + 1) % len(modes)]


def display_mode_label(mode: str) -> str:
    return {
        DISPLAY_MODE_NORMAL: "Normal",
        DISPLAY_MODE_INVERTED_CITY: "Inverted City",
        DISPLAY_MODE_SIMPLE_NO_LABELS: "Simplified: no labels",
        DISPLAY_MODE_SIMPLE_LABELS: "Simplified: with labels",
    }.get(mode, "Normal")
