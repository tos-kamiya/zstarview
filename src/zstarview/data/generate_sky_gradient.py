#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate a sky-gradient chart using the runtime sky-color model."""

import numpy as np
from PIL import Image

from zstarview.render.draw_sky_disc import sky_color_samples

IMG_WIDTH = 300
IMG_HEIGHT = 181
SUN_ALTITUDE_DEG = 60.0
EXPOSURE = 1.0
SATURATION = 1.2
OUTPUT_FILENAME = "sky_gradient_chart.png"


def main() -> None:
    print("Generating sky gradient image...")
    print(f"Settings: Sun Altitude = {SUN_ALTITUDE_DEG} deg")

    img = np.zeros((IMG_HEIGHT, IMG_WIDTH, 3), dtype=np.uint8)
    sun_az = 0.0
    view_azimuths = (sun_az, 90.0, 180.0)
    strip_width = IMG_WIDTH // 3
    view_alt = (90.0 - np.arange(IMG_HEIGHT, dtype=np.float32)).astype(np.float32)

    for i, view_az in enumerate(view_azimuths):
        az = np.full_like(view_alt, float(view_az), dtype=np.float32)
        rgb = sky_color_samples(
            view_alt,
            az,
            (SUN_ALTITUDE_DEG, sun_az),
            exposure=EXPOSURE,
            saturation=SATURATION,
            alpha=1.0,
            eclipse_factor=1.0,
        )
        rgb_u8 = np.clip(np.round(rgb * 255.0), 0, 255).astype(np.uint8)
        x0 = i * strip_width
        x1 = IMG_WIDTH if i == 2 else (i + 1) * strip_width
        img[:, x0:x1, :] = rgb_u8[:, None, :]

    Image.fromarray(img, mode="RGB").save(OUTPUT_FILENAME)
    print(f"Image saved as '{OUTPUT_FILENAME}'")


if __name__ == "__main__":
    main()
