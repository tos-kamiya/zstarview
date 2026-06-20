"""Shared helpers for simplified-view state resolution."""

from __future__ import annotations


def resolve_simplified_view_mode(
    *,
    base_enabled: bool,
    labels_enabled: bool,
) -> str:
    """Return the effective simplified-view mode.

    Modes are:
    - normal
    - nolabels
    - labels
    """
    if not base_enabled:
        return "normal"
    return "labels" if labels_enabled else "nolabels"
