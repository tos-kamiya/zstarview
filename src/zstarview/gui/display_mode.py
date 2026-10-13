"""Display modes cycled by the viewer's Space shortcut."""

from __future__ import annotations

DISPLAY_MODE_NORMAL = "normal"
DISPLAY_MODE_INVERTED_CITY = "inverted-city"
DISPLAY_MODE_SIMPLE_NO_LABELS = "simple-no-labels"
DISPLAY_MODE_SIMPLE_LABELS = "simple-labels"
DISPLAY_MODE_LANDSCAPE = "landscape"

DISPLAY_MODES = (
    DISPLAY_MODE_NORMAL,
    DISPLAY_MODE_INVERTED_CITY,
    DISPLAY_MODE_SIMPLE_NO_LABELS,
    DISPLAY_MODE_SIMPLE_LABELS,
    DISPLAY_MODE_LANDSCAPE,
)
DISPLAY_MODE_ALIASES = {
    "default": DISPLAY_MODE_NORMAL,
    "inverted": DISPLAY_MODE_INVERTED_CITY,
    "simple": DISPLAY_MODE_SIMPLE_NO_LABELS,
    "labels": DISPLAY_MODE_SIMPLE_LABELS,
}


def normalize_display_mode(value: str) -> str:
    """Resolve CLI aliases without accepting unknown modes implicitly."""
    key = value.strip().lower()
    return DISPLAY_MODE_ALIASES.get(key, key)


def display_mode_cycle(*, urban_outline_available: bool) -> tuple[str, ...]:
    if urban_outline_available:
        return DISPLAY_MODES
    return tuple(mode for mode in DISPLAY_MODES if mode != DISPLAY_MODE_INVERTED_CITY)


def next_display_mode(current: str, *, urban_outline_available: bool) -> str:
    modes = display_mode_cycle(urban_outline_available=urban_outline_available)
    try:
        index = modes.index(current)
    except ValueError:
        index = 0
    return modes[(index + 1) % len(modes)]


def previous_display_mode(current: str, *, urban_outline_available: bool) -> str:
    modes = display_mode_cycle(urban_outline_available=urban_outline_available)
    try:
        index = modes.index(current)
    except ValueError:
        index = 0
    return modes[(index - 1) % len(modes)]


def display_mode_label(mode: str) -> str:
    return {
        DISPLAY_MODE_NORMAL: "Normal",
        DISPLAY_MODE_INVERTED_CITY: "Inverted City",
        DISPLAY_MODE_SIMPLE_NO_LABELS: "Simplified: no labels",
        DISPLAY_MODE_SIMPLE_LABELS: "Simplified: with labels",
        DISPLAY_MODE_LANDSCAPE: "Landscape view",
    }.get(mode, "Normal")
