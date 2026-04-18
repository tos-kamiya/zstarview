#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Render a sky-disc PNG using the runtime renderer."""

import argparse

from PySide6.QtCore import QRect, Qt
from PySide6.QtGui import QColor, QImage, QPainter

from zstarview.render.geometry import get_screen_geometry
from zstarview.render.sky_disc import draw_sky_color_disc


def main() -> None:
    ap = argparse.ArgumentParser(description="Render a sky-disc test image using runtime code.")
    ap.add_argument("--width", type=int, default=320)
    ap.add_argument("--height", type=int, default=320)
    ap.add_argument("--sun-alt", type=float, default=20.0)
    ap.add_argument("--sun-az", type=float, default=270.0, help="0=N, 90=E")
    ap.add_argument("--center-alt", type=float, default=90.0)
    ap.add_argument("--center-az", type=float, default=0.0)
    ap.add_argument("--exposure", type=float, default=1.0)
    ap.add_argument("--saturation", type=float, default=1.2)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--outfile", type=str, default="sky_test.png")
    args = ap.parse_args()

    w = max(2, int(args.width))
    h = max(2, int(args.height))
    view_center = (float(args.center_alt), float(args.center_az))
    sun_altaz = (float(args.sun_alt), float(args.sun_az))

    geometry = get_screen_geometry(w, h, view_center[0])
    disc = draw_sky_color_disc(
        geometry,
        view_center=view_center,
        sun_altaz=sun_altaz,
        alpha=float(args.alpha),
        exposure=float(args.exposure),
        saturation=float(args.saturation),
        eclipse_factor=1.0,
        content_fov_deg=90.0,
    )

    canvas = QImage(w, h, QImage.Format.Format_ARGB32_Premultiplied)
    canvas.fill(QColor(0, 0, 0, 255))
    painter = QPainter(canvas)
    painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)

    x = int(geometry.center[0] - geometry.radius)
    y = int(geometry.center[1] - geometry.radius)
    d = int(geometry.radius * 2)
    painter.drawImage(QRect(x, y, d, d), disc)
    painter.setPen(Qt.PenStyle.NoPen)
    painter.end()

    ok = canvas.save(args.outfile, "PNG")
    print(f"Saved: {args.outfile}, Success: {ok}")


if __name__ == "__main__":
    main()
