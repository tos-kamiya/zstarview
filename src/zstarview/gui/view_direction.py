from __future__ import annotations

from PySide6.QtCore import Qt

from ..paths import OBSERVER_MAX_ALT_DEG, OBSERVER_MIN_ALT_DEG

DEFAULT_VIEW_DIRECTION_STEP_DEG = 5.0
FINE_VIEW_DIRECTION_STEP_DEG = 1.0


def resolve_view_direction_step(
    modifiers: Qt.KeyboardModifier,
    coarse_step_deg: float = DEFAULT_VIEW_DIRECTION_STEP_DEG,
) -> float:
    if modifiers & Qt.KeyboardModifier.ShiftModifier:
        return FINE_VIEW_DIRECTION_STEP_DEG
    return float(coarse_step_deg)


def clamp_view_center_alt_az(
    alt_deg: float,
    az_deg: float,
) -> tuple[float, float]:
    alt = max(OBSERVER_MIN_ALT_DEG, min(OBSERVER_MAX_ALT_DEG, float(alt_deg)))
    az = float(az_deg) % 360.0
    return alt, az
