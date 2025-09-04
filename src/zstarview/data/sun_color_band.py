#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates a color gradient chart for the sunlight color model.

This script renders a horizontal color band that visualizes the output of the
`get_sun_color` function over the full range of sun altitudes from -90° to +90°.
The resulting image is useful for debugging and tuning the sky color algorithm.
"""

import argparse
import math
from pathlib import Path
from typing import Tuple

import numpy as np
from PIL import Image, ImageDraw, ImageFont


def _lerp_color(a: np.ndarray, b: np.ndarray, t: float) -> np.ndarray:
    """Linearly interpolates between two colors represented as NumPy arrays."""
    return a + (b - a) * t


def get_sun_color(sun_alt_deg: float) -> np.ndarray:
    """
    Determines the color of sunlight based on the sun's altitude.

    Args:
        sun_alt_deg: The sun's altitude in degrees (-90 to +90).

    Returns:
        A NumPy array representing the (r, g, b) color in the range [0, 1].
    """
    # Define key colors for different sun altitudes.
    zenith_color = np.array([0.3, 0.55, 0.98])  # Deep blue for when the sun is at zenith.
    horizon_color = np.array([0.95, 0.50, 0.30])  # Orange for when the sun is at the horizon.
    night_color = np.array([0.01, 0.02, 0.05])  # Dark blue for night.

    # --- Calculate daytime color ---
    # Normalize sun altitude from -7° (below horizon) to 90° (zenith) to a [0,1] range.
    t = float(np.clip((sun_alt_deg + 7.0) / 97.0, 0.0, 1.0))
    # Apply a power function to create a more natural, non-linear transition.
    t = math.pow(t, 0.35)
    day_color = _lerp_color(horizon_color, zenith_color, t)

    # --- Fade to night color ---
    # Create a rapid fade to darkness as the sun sets (from +1° to -5°).
    fade = float(np.clip((-sun_alt_deg + 1.0) / 6.0, 0.0, 1.0))
    return _lerp_color(day_color, night_color, fade)


def _text_size(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont) -> Tuple[int, int]:
    """
    Calculates the bounding box size of a text string, handling Pillow version differences.

    Pillow versions 9.2.0 and later use `draw.textbbox`, while older versions
    use `draw.textsize`. This function provides a compatible way to get the text size.

    Args:
        draw: The ImageDraw object.
        text: The text string to measure.
        font: The font to use.

    Returns:
        A tuple (width, height) of the rendered text.
    """
    try:
        # Modern Pillow versions
        bbox = draw.textbbox((0, 0), text, font=font)
        return bbox[2] - bbox[0], bbox[3] - bbox[1]
    except AttributeError:
        # Older Pillow versions
        return draw.textsize(text, font=font)


def render_sun_color_band(
    width_px: int = 1800,
    height_px: int = 160,
    margin_px: int = 40,
    outfile: str = "sun_color_band.png",
    show_ticks: bool = True,
    show_title: bool = True,
    transparent: bool = False,
) -> str:
    """
    Renders a horizontal color band visualizing sunlight color vs. sun altitude.

    The band covers the range from -90 to +90 degrees. It can optionally include
    a title and tick marks for reference.

    Args:
        width_px: The total width of the output image.
        height_px: The total height of the output image.
        margin_px: The margin around the color band.
        outfile: The path to save the output PNG file.
        show_ticks: Whether to draw the altitude tick marks and axis.
        show_title: Whether to draw the title text.
        transparent: If True, the background will be transparent; otherwise, it will be white.

    Returns:
        The path to the saved output file.
    """
    # Create a transparent RGBA canvas to draw on.
    img = Image.new("RGBA", (int(width_px), int(height_px)), (255, 255, 255, 0))
    draw = ImageDraw.Draw(img)

    # Define the rectangle for the color band.
    x0, y0 = margin_px, margin_px
    x1, y1 = width_px - margin_px, height_px - margin_px
    band_w = max(1, x1 - x0)

    # Fill the band by drawing a vertical line for each pixel, corresponding to a sun altitude.
    for i in range(band_w):
        t = i / (band_w - 1)  # Normalize pixel position to [0, 1]
        sun_alt = -90.0 + t * 180.0  # Map to altitude range [-90, 90]
        rgb = np.clip(get_sun_color(sun_alt), 0.0, 1.0)
        color = tuple(int(round(c * 255.0)) for c in rgb)
        draw.line([(x0 + i, y0), (x0 + i, y1)], fill=color, width=1)

    # --- Draw annotations (ticks and title) ---
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 14)
    except IOError:
        font = ImageFont.load_default()

    if show_ticks:
        axis_y = y1 + 5
        if axis_y < height_px - 1:
            draw.line([(x0, axis_y), (x1, axis_y)], fill=(0, 0, 0, 255), width=1)
            for deg in range(-90, 91, 30):
                t = (deg + 90.0) / 180.0
                x = int(round(x0 + t * (band_w - 1)))
                draw.line([(x, axis_y), (x, axis_y + 6)], fill=(0, 0, 0, 255), width=1)
                label = f"{deg}°"
                tw, th = _text_size(draw, label, font)
                if axis_y + 8 + th < height_px:
                    draw.text(
                        (x - tw // 2, axis_y + 8),
                        label,
                        fill=(0, 0, 0, 255),
                        font=font,
                    )

    if show_title:
        title = "Sunlight Color vs. Sun Altitude (−90° to +90°)"
        tw, th = _text_size(draw, title, font)
        draw.text(
            ((width_px - tw) // 2, max(2, (margin_px - th) // 2)),
            title,
            fill=(0, 0, 0, 255),
            font=font,
        )

    # --- Save the final image ---
    if transparent:
        img.save(outfile, format="PNG")
    else:
        # Create a white background and paste the RGBA image onto it.
        bg = Image.new("RGB", (width_px, height_px), (255, 255, 255))
        bg.paste(img, mask=img.split()[-1])  # Use alpha channel as mask
        bg.save(outfile, format="PNG")

    return outfile


def main() -> None:
    """Parses command-line arguments and runs the rendering function."""
    p = argparse.ArgumentParser(description="Render a color band showing sunlight color vs. sun altitude.")
    p.add_argument("--outfile", type=Path, default=Path("sun_color_band.png"))
    p.add_argument("--width", type=int, default=1800)
    p.add_argument("--height", type=int, default=160)
    p.add_argument("--margin", type=int, default=40)
    p.add_argument("--no-ticks", action="store_true", help="Hide ticks and axis labels.")
    p.add_argument("--no-title", action="store_true", help="Hide the title.")
    p.add_argument("--transparent", action="store_true", help="Save with a transparent background.")
    args = p.parse_args()

    render_sun_color_band(
        width_px=args.width,
        height_px=args.height,
        margin_px=args.margin,
        outfile=str(args.outfile),
        show_ticks=not args.no_ticks,
        show_title=not args.no_title,
        transparent=args.transparent,
    )
    print(f"Saved color band to '{args.outfile}'")


if __name__ == "__main__":
    main()
