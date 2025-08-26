#!/usr/bin/env python3
import math
import argparse
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont

def _lerp_color(a: np.ndarray, b: np.ndarray, t: float) -> np.ndarray:
    """Linearly interpolates between two colors."""
    return a + (b - a) * t

def get_sun_color(sun_alt_deg: float) -> np.ndarray:
    """
    Determine the color of sunlight based on the sun's altitude.
    Returns (r, g, b) in the range [0,1].
    """
    zenith_color  = np.array([0.3, 0.55, 0.98])   # zenith (blue)
    horizon_color = np.array([0.95, 0.50, 0.30])   # horizon (orange)
    night_color   = np.array([0.01, 0.02, 0.05])  # night (dark blue)

    # Normalize sun altitude from -7° (sunset) to 90° (zenith) → [0,1], with gamma-ish shaping
    t = float(np.clip((sun_alt_deg + 7.0) / 97.0, 0.0, 1.0))
    t = math.pow(t, 0.35)

    day_color = _lerp_color(horizon_color, zenith_color, t)

    # Rapid darkening near/below horizon
    fade = float(np.clip((-sun_alt_deg + 1.0) / 6.0, 0.0, 1.0))  # between +1° and -8°
    return _lerp_color(day_color, night_color, fade)

def _text_size(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont):
    """Pillowのバージョン差を吸収してテキストサイズを返す"""
    try:
        bbox = draw.textbbox((0, 0), text, font=font)
        return bbox[2] - bbox[0], bbox[3] - bbox[1]
    except Exception:
        return draw.textsize(text, font=font)

def render_sun_color_band(
    width_px: int = 1800,
    height_px: int = 160,
    margin_px: int = 40,
    outfile: str = "sun_color_band.png",
    show_ticks: bool = True,
    show_title: bool = True,
    transparent: bool = False,
):
    """
    Render a horizontal color band for sun_alt_deg in [-90, 90] and save as PNG.
    """
    # 下地はRGBAで作って、必要なら最後に白背景へ合成
    img = Image.new("RGBA", (int(width_px), int(height_px)), (255, 255, 255, 0))
    draw = ImageDraw.Draw(img)

    # 帯の矩形
    x0, y0 = margin_px, margin_px
    x1, y1 = width_px - margin_px, height_px - margin_px
    band_w = max(1, x1 - x0)
    band_h = max(1, y1 - y0)

    # 縦線で塗りつぶし（x を -90..90 に線形対応）
    for i in range(band_w):
        t = i / (band_w - 1)
        sun_alt = -90.0 + t * 180.0
        rgb = np.clip(get_sun_color(sun_alt), 0.0, 1.0)
        color = tuple(int(round(c * 255.0)) for c in rgb)
        draw.line([(x0 + i, y0), (x0 + i, y1)], fill=color, width=1)

    # 目盛り・軸
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 14)
    except Exception:
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
                    draw.text((x - tw // 2, axis_y + 8), label, fill=(0, 0, 0, 255), font=font)

    if show_title:
        title = "Sunlight Color vs. Sun Altitude (−90° to +90°)"
        tw, th = _text_size(draw, title, font)
        draw.text(((width_px - tw) // 2, max(2, (margin_px - th) // 2)),
                  title, fill=(0, 0, 0, 255), font=font)

    # 出力
    if transparent:
        img.save(outfile, format="PNG")
    else:
        bg = Image.new("RGB", (width_px, height_px), (255, 255, 255))
        bg.paste(img, mask=img.split()[-1])
        bg.save(outfile, format="PNG")
    return outfile

def main():
    p = argparse.ArgumentParser(description="Render sunlight color band for sun_alt in [-90, 90].")
    p.add_argument("--outfile", type=Path, default=Path("sun_color_band.png"))
    p.add_argument("--width", type=int, default=1800)
    p.add_argument("--height", type=int, default=160)
    p.add_argument("--margin", type=int, default=40)
    p.add_argument("--no-ticks", action="store_true", help="Hide ticks/axis.")
    p.add_argument("--no-title", action="store_true", help="Hide title.")
    p.add_argument("--transparent", action="store_true", help="Save with transparent background.")
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

if __name__ == "__main__":
    main()
