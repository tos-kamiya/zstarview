#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Render a band chart of runtime sky color vs sun altitude."""

import argparse
from pathlib import Path
from typing import Tuple

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from zstarview.render.draw_sky_disc import sky_color_samples


def _text_size(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont) -> Tuple[int, int]:
    try:
        bbox = draw.textbbox((0, 0), text, font=font)
        return bbox[2] - bbox[0], bbox[3] - bbox[1]
    except AttributeError:
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
    img = Image.new("RGBA", (int(width_px), int(height_px)), (255, 255, 255, 0))
    draw = ImageDraw.Draw(img)

    x0, y0 = margin_px, margin_px
    x1, y1 = width_px - margin_px, height_px - margin_px
    band_w = max(1, x1 - x0)

    for i in range(band_w):
        t = i / max(1, band_w - 1)
        sun_alt = -90.0 + t * 180.0
        view_alt = max(0.0, sun_alt)
        rgb = sky_color_samples(
            np.array([view_alt], dtype=np.float32),
            np.array([0.0], dtype=np.float32),
            (float(sun_alt), 0.0),
            exposure=1.0,
            saturation=1.0,
            alpha=1.0,
            eclipse_factor=1.0,
        )[0]
        color = tuple(int(round(float(c) * 255.0)) for c in rgb)
        draw.line([(x0 + i, y0), (x0 + i, y1)], fill=color, width=1)

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
                label = f"{deg}deg"
                tw, th = _text_size(draw, label, font)
                if axis_y + 8 + th < height_px:
                    draw.text((x - tw // 2, axis_y + 8), label, fill=(0, 0, 0, 255), font=font)

    if show_title:
        title = "Runtime Sky Color Toward Sun vs Sun Altitude"
        tw, th = _text_size(draw, title, font)
        draw.text(((width_px - tw) // 2, max(2, (margin_px - th) // 2)), title, fill=(0, 0, 0, 255), font=font)

    if transparent:
        img.save(outfile, format="PNG")
    else:
        bg = Image.new("RGB", (width_px, height_px), (255, 255, 255))
        bg.paste(img, mask=img.split()[-1])
        bg.save(outfile, format="PNG")
    return outfile


def main() -> None:
    p = argparse.ArgumentParser(description="Render a runtime sky-color band vs sun altitude.")
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
