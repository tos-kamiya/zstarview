from __future__ import annotations

import argparse

import numpy as np
from PIL import Image

DEFAULT_LOGO_MASK_FRAC = (0.16, 0.10)
DEFAULT_CONTRAST_LOW = 5.0
DEFAULT_CONTRAST_HIGH = 98.0


def _parse_fraction_pair(value: str) -> tuple[float, float]:
    text = (value or "").strip()
    parts = [part.strip() for part in text.split(",")]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("Expected two comma-separated fractions, e.g. '0.16,0.10'.")
    try:
        x = float(parts[0])
        y = float(parts[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Fractions must be numeric.") from exc
    if not (0.0 <= x <= 1.0 and 0.0 <= y <= 1.0):
        raise argparse.ArgumentTypeError("Fractions must be between 0 and 1.")
    return (x, y)


def _apply_logo_mask(
    score: np.ndarray,
    *,
    mask_frac: tuple[float, float],
    mask_rect: tuple[int, int, int, int] | None,
) -> None:
    if mask_rect is not None:
        x0, y0, w, h = mask_rect
        x1 = max(x0, min(score.shape[1], x0 + w))
        y1 = max(y0, min(score.shape[0], y0 + h))
        score[max(0, y0):y1, max(0, x0):x1] = 0.0
        return

    mask_w = max(1, int(round(score.shape[1] * float(mask_frac[0]))))
    mask_h = max(1, int(round(score.shape[0] * float(mask_frac[1]))))
    score[:mask_h, :mask_w] = 0.0


def build_cloud_proxy(
    image: Image.Image,
    *,
    mask_logo: bool = True,
    mask_frac: tuple[float, float] = DEFAULT_LOGO_MASK_FRAC,
    mask_rect: tuple[int, int, int, int] | None = None,
    contrast_low: float = DEFAULT_CONTRAST_LOW,
    contrast_high: float = DEFAULT_CONTRAST_HIGH,
) -> Image.Image:
    rgb = np.asarray(image.convert("RGB"), dtype=np.float32) / 255.0
    r = rgb[..., 0]
    g = rgb[..., 1]
    b = rgb[..., 2]

    max_rgb = np.maximum.reduce([r, g, b])
    min_rgb = np.minimum.reduce([r, g, b])
    spread = max_rgb - min_rgb
    brightness = 0.2126 * r + 0.7152 * g + 0.0722 * b
    whiteness = 1.0 - spread
    saturation = np.where(max_rgb > 1e-6, spread / np.maximum(max_rgb, 1e-6), 0.0)
    blue_dominance = np.clip(b - np.maximum(r, g), 0.0, 1.0)
    green_dominance = np.clip(g - np.maximum(r, b), 0.0, 1.0)

    score = (
        0.55 * brightness
        + 0.45 * whiteness
        - 0.70 * saturation
        - 0.30 * blue_dominance
        - 0.20 * green_dominance
    )
    score = np.clip(score, 0.0, 1.0)

    if mask_logo:
        _apply_logo_mask(score, mask_frac=mask_frac, mask_rect=mask_rect)

    finite = score[np.isfinite(score)]
    if finite.size > 0:
        lo = float(np.percentile(finite, float(contrast_low)))
        hi = float(np.percentile(finite, float(contrast_high)))
        if hi > lo + 1e-6:
            score = np.clip((score - lo) / (hi - lo), 0.0, 1.0)

    return Image.fromarray((score * 255.0 + 0.5).astype(np.uint8), mode="L")
