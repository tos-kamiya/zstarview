import math
from typing import Optional, Tuple

import numpy as np
from PySide6.QtGui import QColor

from ..paths import TEXT_STYLES_BY_PRESET


def bv_to_rgb_vectorized(bv: np.ndarray) -> np.ndarray:
    """Vectorized conversion of B-V color index to RGB tuples."""
    rgb = np.full((bv.shape[0], 3), [255, 204, 111], dtype=int)
    rgb[bv < 0.0] = [170, 191, 255]
    rgb[(bv >= 0.0) & (bv < 0.3)] = [202, 215, 255]
    rgb[(bv >= 0.3) & (bv < 0.6)] = [248, 247, 255]
    rgb[(bv >= 0.6) & (bv < 1.0)] = [255, 210, 161]
    return rgb


def flare_strength_from_vmag(vmag: float) -> float:
    """Return monotonic flare strength [0, 1] for all magnitudes."""
    vmag_bright = -1.5
    vmag_faint = 6.0
    t = (vmag_faint - float(vmag)) / (vmag_faint - vmag_bright)
    t = max(0.0, min(1.0, t))
    return t**1.35


def compute_flare_profile(vmag: float, core_radius_px: float) -> Tuple[float, float]:
    """Compute (core_scale, flare_outer_px) for a star."""
    strength = flare_strength_from_vmag(vmag)
    flare_outer_px = float(core_radius_px) * (0.65 * strength)
    if flare_outer_px < 1.0:
        return 1.0, 0.0
    core_scale = 1.0 / math.sqrt(1.0 + 0.9 * strength)
    return core_scale, flare_outer_px


def planet_disc_style_from_vmag(vmag: Optional[float]) -> Tuple[float, int]:
    """Return (radius_px, alpha) for a planet disc marker."""
    if vmag is None or not math.isfinite(float(vmag)):
        return 3.0, 200
    vmag_clipped = float(np.clip(float(vmag), -1.5, 6.0))
    strength = flare_strength_from_vmag(vmag_clipped)
    radius_px = 2.4 + 3.2 * strength
    alpha = int(np.clip(round(125 + 130 * strength), 110, 255))
    return radius_px, alpha


def planet_bloom_profile_from_vmag(
    vmag: Optional[float], core_radius_px: float
) -> Tuple[float, int, int]:
    """Return bloom profile as (radius_px, center_alpha, mid_alpha)."""
    if vmag is None or not math.isfinite(float(vmag)):
        return 0.0, 0, 0

    vm = float(np.clip(float(vmag), -6.0, 6.0))
    base = flare_strength_from_vmag(float(np.clip(vm, -1.5, 6.0)))
    extra = max(0.0, min(1.0, (-1.5 - vm) / 4.5))
    strength = 0.60 * base + 0.40 * extra

    if strength < 0.12:
        return 0.0, 0, 0

    r_core = max(1.0, float(core_radius_px))
    bloom_radius = r_core * (1.20 + 1.20 * strength)
    center_alpha = int(np.clip(round(10 + 58 * strength), 8, 90))
    mid_alpha = int(np.clip(round(4 + 34 * strength), 3, 60))
    return bloom_radius, center_alpha, mid_alpha


def planet_marker_color(name: str) -> QColor:
    """Return display color for a planet marker."""
    palette = {
        "mercury": QColor(190, 190, 182),
        "venus": QColor(245, 226, 176),
        "mars": QColor(232, 126, 96),
        "jupiter": QColor(224, 188, 141),
        "saturn": QColor(226, 214, 154),
        "uranus": QColor(157, 224, 218),
        "neptune": QColor(108, 152, 234),
        "pluto": QColor(194, 166, 132),
    }
    return QColor(palette.get(name, QColor(*TEXT_STYLES_BY_PRESET["night"].text)))


def small_body_marker_color(name: str) -> QColor:
    """Return a color for a dwarf planet or other named small body marker."""
    normalized = name.strip().casefold()
    palette = {
        "ceres": QColor(122, 182, 118),
        "eris": QColor(123, 198, 226),
        "haumea": QColor(216, 178, 128),
        "makemake": QColor(230, 188, 112),
        "pluto": planet_marker_color("pluto"),
    }
    if normalized in palette:
        return QColor(palette[normalized])
    return QColor(196, 176, 132)


def body_label_text(name: str) -> str:
    """Return a display label for a solar-system body name."""
    label = name.strip()
    if not label:
        return label
    return label.title()
