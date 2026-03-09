#!/usr/bin/env python3
"""Render a comparison image for gray-background star color experiments.

This stays outside the main app and reuses the runtime B-V -> RGB mapping so
you can judge whether boosting saturation helps stars stand out on a gray
background.

Example:
  python scripts/render_star_gray_background_experiment.py
  python scripts/render_star_gray_background_experiment.py --gray-levels 80 96 112
  python scripts/render_star_gray_background_experiment.py --gray-levels 96 --sat-gain 1.7
  python scripts/render_star_gray_background_experiment.py --gray-levels 144 --outline-vmag-max 2.0
"""

from __future__ import annotations

import argparse
import colorsys
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QColor, QFont, QImage, QPainter, QPen

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.render.photometry import bv_to_rgb_vectorized  # noqa: E402


MAG2_TO_MAG1_SIZE_SCALE = 10.0 ** 0.12
DEFAULT_BV_SAMPLES = (-0.1, 0.15, 0.45, 0.8, 1.3)
DEFAULT_VMAG_SAMPLES = (0.0, 1.0, 2.0, 3.0, 4.0, 5.0)


@dataclass(frozen=True)
class PanelSpec:
    title: str
    background: QColor
    sat_gain: float
    boost_mode: str
    outline_enabled: bool


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Render a synthetic star comparison for gray-background experiments."
    )
    parser.add_argument("--outfile", default="star_gray_background_experiment.png")
    parser.add_argument(
        "--gray-levels",
        type=int,
        nargs="+",
        default=[80, 96, 112],
        help="Gray background levels to render as separate rows (default: 80 96 112).",
    )
    parser.add_argument("--star-base-radius", type=float, default=4.0, help="Same meaning as the app option.")
    parser.add_argument("--vmag-brightness-scale", type=float, default=-0.39)
    parser.add_argument("--sat-gain", type=float, default=1.55, help="Saturation multiplier for the boosted gray panel.")
    parser.add_argument(
        "--visibility-boost",
        type=float,
        default=1.0,
        help="Additional brightness multiplier applied after the app-like color falloff.",
    )
    parser.add_argument(
        "--priority-vmag",
        type=float,
        default=2.0,
        help="Center magnitude for the bright-star-priority boost rolloff.",
    )
    parser.add_argument(
        "--priority-width",
        type=float,
        default=2.0,
        help="How many magnitudes it takes the priority boost to fade toward baseline.",
    )
    parser.add_argument(
        "--outline-vmag-max",
        type=float,
        default=2.0,
        help="Only draw outline for stars at or brighter than this Vmag.",
    )
    return parser


def _compute_star_profiles(vmag: np.ndarray, scale: float) -> tuple[np.ndarray, np.ndarray]:
    v_ref = 1.0
    raw_luminance = 10.0 ** (float(scale) * (vmag - v_ref))
    size_scale = np.power(np.clip(raw_luminance, 0.0, 1.0), 0.3)
    color_factor_base = np.clip(np.power(raw_luminance, 0.6), 0.0, 1.0)
    return size_scale, color_factor_base


def _boost_rgb(rgb: Sequence[int], sat_gain: float) -> tuple[int, int, int]:
    r, g, b = (channel / 255.0 for channel in rgb)
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    s = min(1.0, s * sat_gain)
    rr, gg, bb = colorsys.hsv_to_rgb(h, s, v)
    return (
        int(round(rr * 255.0)),
        int(round(gg * 255.0)),
        int(round(bb * 255.0)),
    )


def _star_color(
    base_rgb: Sequence[int],
    vmag: float,
    color_factor_base: float,
    visibility_boost: float,
    sat_gain: float,
    boost_mode: str,
    priority_vmag: float,
    priority_width: float,
) -> QColor:
    rgb = tuple(int(c) for c in base_rgb)
    mode_sat_gain = sat_gain
    if boost_mode == "priority":
        width = max(0.2, float(priority_width))
        reference_vmag = float(priority_vmag)
        weight = 1.0 - np.clip((reference_vmag - float(priority_vmag)) / width, 0.0, 1.0)
        mode_sat_gain = 1.0 + (sat_gain - 1.0) * float(weight)
    if mode_sat_gain != 1.0:
        rgb = _boost_rgb(rgb, sat_gain=mode_sat_gain)
    color_factor = np.clip((0.15 + 0.85 * float(color_factor_base)) * float(visibility_boost), 0.0, 1.0)
    scaled = tuple(int(round(channel * color_factor)) for channel in rgb)
    return QColor(*scaled, 255)


def _outline_color(color: QColor) -> QColor:
    h, s, v, _ = color.getHsv()
    sat = max(18, int(s * 0.42))
    val = max(18, int(v * 0.38))
    outline = QColor()
    outline.setHsv(h if h >= 0 else 0, sat, val, 220)
    return outline


def _draw_star_shape(painter: QPainter, center: QPointF, size_px: float, color: QColor, vmag: float) -> None:
    half = max(0.75, size_px * 0.5)
    rect = QRectF(center.x() - half, center.y() - half, half * 2.0, half * 2.0)
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(color)
    painter.drawEllipse(rect)
    if vmag < 2.0:
        diamond = half * (1.15 + max(0.0, 2.0 - vmag) * 0.14)
        points = [
            QPointF(center.x(), center.y() - diamond),
            QPointF(center.x() + diamond, center.y()),
            QPointF(center.x(), center.y() + diamond),
            QPointF(center.x() - diamond, center.y()),
        ]
        painter.drawPolygon(points)


def _draw_star_outline(
    painter: QPainter,
    center: QPointF,
    size_px: float,
    color: QColor,
    vmag: float,
    outline_vmag_max: float,
) -> None:
    if size_px < 3.0 or float(vmag) > float(outline_vmag_max):
        return
    half = max(0.75, size_px * 0.5)
    outer_half = half + 1.0
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.setPen(QPen(_outline_color(color), 1.0))
    painter.drawEllipse(QRectF(center.x() - outer_half, center.y() - outer_half, outer_half * 2.0, outer_half * 2.0))
    if vmag < 2.0:
        diamond = outer_half * (1.15 + max(0.0, 2.0 - vmag) * 0.14)
        points = [
            QPointF(center.x(), center.y() - diamond),
            QPointF(center.x() + diamond, center.y()),
            QPointF(center.x(), center.y() + diamond),
            QPointF(center.x() - diamond, center.y()),
        ]
        painter.drawPolygon(points)


def _panel_specs(gray_level: int, sat_gain: float) -> list[PanelSpec]:
    gray = max(0, min(255, int(gray_level)))
    return [
        PanelSpec("Black / current", QColor(0, 0, 0), 1.0, "none", False),
        PanelSpec(f"Gray {gray} / current", QColor(gray, gray, gray), 1.0, "none", False),
        PanelSpec(f"Gray {gray} / uniform boost", QColor(gray, gray, gray), sat_gain, "uniform", False),
        PanelSpec(f"Gray {gray} / uniform + outline", QColor(gray, gray, gray), sat_gain, "uniform", True),
        PanelSpec(f"Gray {gray} / bright-star boost", QColor(gray, gray, gray), sat_gain, "priority", False),
        PanelSpec(f"Gray {gray} / boost + outline", QColor(gray, gray, gray), sat_gain, "priority", True),
    ]


def _draw_panel(
    painter: QPainter,
    panel_rect: QRectF,
    panel: PanelSpec,
    bv_values: Sequence[float],
    vmag_values: Sequence[float],
    *,
    star_base_radius: float,
    vmag_brightness_scale: float,
    visibility_boost: float,
    priority_vmag: float,
    priority_width: float,
    outline_vmag_max: float,
) -> None:
    painter.fillRect(panel_rect, panel.background)

    label_color = QColor(232, 236, 242) if panel.background.lightness() < 90 else QColor(22, 24, 28)
    grid_color = QColor(label_color)
    grid_color.setAlpha(48 if panel.background.lightness() < 128 else 36)
    painter.setPen(QPen(grid_color, 1))

    left_pad = 88.0
    top_pad = 72.0
    inner = QRectF(
        panel_rect.left() + left_pad,
        panel_rect.top() + top_pad,
        panel_rect.width() - left_pad - 18.0,
        panel_rect.height() - top_pad - 18.0,
    )
    col_step = inner.width() / max(1, len(bv_values))
    row_step = inner.height() / max(1, len(vmag_values))

    title_font = QFont("Sans Serif", 12)
    title_font.setBold(True)
    painter.setFont(title_font)
    painter.setPen(label_color)
    painter.drawText(
        QRectF(panel_rect.left() + 14.0, panel_rect.top() + 10.0, panel_rect.width() - 28.0, 26.0),
        int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter),
        panel.title,
    )

    note_font = QFont("Sans Serif", 9)
    painter.setFont(note_font)
    painter.drawText(
        QRectF(panel_rect.left() + 14.0, panel_rect.top() + 36.0, panel_rect.width() - 28.0, 18.0),
        int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter),
        (
            f"sat x{panel.sat_gain:.2f}"
            if panel.boost_mode != "priority"
            else f"sat x{panel.sat_gain:.2f}, fixed@Vmag{priority_vmag:.1f}"
        ),
    )

    painter.setPen(QPen(grid_color, 1))
    for col in range(len(bv_values) + 1):
        x = inner.left() + col * col_step
        painter.drawLine(QPointF(x, inner.top()), QPointF(x, inner.bottom()))
    for row in range(len(vmag_values) + 1):
        y = inner.top() + row * row_step
        painter.drawLine(QPointF(inner.left(), y), QPointF(inner.right(), y))

    body_font = QFont("Sans Serif", 10)
    painter.setFont(body_font)
    painter.setPen(label_color)

    for i, bv in enumerate(bv_values):
        text_rect = QRectF(inner.left() + i * col_step, inner.top() - 24.0, col_step, 18.0)
        painter.drawText(text_rect, int(Qt.AlignmentFlag.AlignCenter), f"B-V {bv:+.2f}")

    for i, vmag in enumerate(vmag_values):
        text_rect = QRectF(panel_rect.left() + 12.0, inner.top() + i * row_step, left_pad - 18.0, row_step)
        painter.drawText(text_rect, int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter), f"Vmag {vmag:.1f}")

    star_layer = QImage(
        math.ceil(panel_rect.width()),
        math.ceil(panel_rect.height()),
        QImage.Format.Format_ARGB32_Premultiplied,
    )
    star_layer.fill(QColor(0, 0, 0, 0))
    outline_layer = QImage(
        math.ceil(panel_rect.width()),
        math.ceil(panel_rect.height()),
        QImage.Format.Format_ARGB32_Premultiplied,
    )
    outline_layer.fill(QColor(0, 0, 0, 0))

    layer_painter = QPainter(star_layer)
    layer_painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    layer_painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    outline_painter = QPainter(outline_layer)
    outline_painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    outline_painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)

    bv_arr = np.asarray(bv_values, dtype=float)
    base_rgb = bv_to_rgb_vectorized(bv_arr)

    for row, vmag in enumerate(vmag_values):
        size_scale, color_factor_base = _compute_star_profiles(np.asarray([vmag], dtype=float), vmag_brightness_scale)
        size_px = max(1.0, star_base_radius * MAG2_TO_MAG1_SIZE_SCALE * float(size_scale[0]))
        for col, rgb in enumerate(base_rgb):
            center = QPointF(
                (inner.left() - panel_rect.left()) + (col + 0.5) * col_step,
                (inner.top() - panel_rect.top()) + (row + 0.5) * row_step,
            )
            star_color = _star_color(
                rgb,
                vmag=float(vmag),
                color_factor_base=float(color_factor_base[0]),
                visibility_boost=visibility_boost,
                sat_gain=panel.sat_gain,
                boost_mode=panel.boost_mode,
                priority_vmag=priority_vmag,
                priority_width=priority_width,
            )
            if panel.outline_enabled:
                _draw_star_outline(
                    outline_painter,
                    center,
                    size_px,
                    star_color,
                    vmag,
                    outline_vmag_max,
                )
            _draw_star_shape(layer_painter, center, size_px, star_color, vmag)

    layer_painter.end()
    outline_painter.end()

    if panel.outline_enabled:
        painter.drawImage(panel_rect.topLeft(), outline_layer)

    painter.save()
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
    painter.drawImage(panel_rect.topLeft(), star_layer)
    painter.restore()


def render_experiment(
    outfile: Path,
    *,
    gray_levels: Iterable[int],
    star_base_radius: float,
    vmag_brightness_scale: float,
    sat_gain: float,
    visibility_boost: float,
    priority_vmag: float,
    priority_width: float,
    outline_vmag_max: float,
    bv_values: Iterable[float] = DEFAULT_BV_SAMPLES,
    vmag_values: Iterable[float] = DEFAULT_VMAG_SAMPLES,
) -> None:
    gray_levels = tuple(max(0, min(255, int(level))) for level in gray_levels)
    if not gray_levels:
        raise ValueError("At least one gray level is required.")
    bv_values = tuple(float(v) for v in bv_values)
    vmag_values = tuple(float(v) for v in vmag_values)

    panel_width = 360
    panel_height = 420
    gutter = 16
    row_gutter = 20
    row_count = len(gray_levels)
    panels_per_row = 6
    width = panels_per_row * panel_width + (panels_per_row - 1) * gutter
    height = row_count * panel_height + max(0, row_count - 1) * row_gutter

    canvas = QImage(width, height, QImage.Format.Format_ARGB32_Premultiplied)
    canvas.fill(QColor(245, 245, 245))

    painter = QPainter(canvas)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)

    y = 0.0
    for gray_level in gray_levels:
        panels = _panel_specs(gray_level=gray_level, sat_gain=sat_gain)
        x = 0.0
        for panel in panels:
            rect = QRectF(x, y, panel_width, panel_height)
            _draw_panel(
                painter,
                rect,
                panel,
                bv_values,
                vmag_values,
                star_base_radius=star_base_radius,
                vmag_brightness_scale=vmag_brightness_scale,
                visibility_boost=visibility_boost,
                priority_vmag=priority_vmag,
                priority_width=priority_width,
                outline_vmag_max=outline_vmag_max,
            )
            x += panel_width + gutter
        y += panel_height + row_gutter

    painter.end()
    outfile.parent.mkdir(parents=True, exist_ok=True)
    if not canvas.save(str(outfile), "PNG"):
        raise RuntimeError(f"Failed to save output PNG: {outfile}")


def main() -> int:
    args = _build_parser().parse_args()
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    app = QApplication.instance() or QApplication([])
    render_experiment(
        Path(args.outfile),
        gray_levels=args.gray_levels,
        star_base_radius=max(1.0, float(args.star_base_radius)),
        vmag_brightness_scale=float(args.vmag_brightness_scale),
        sat_gain=max(1.0, float(args.sat_gain)),
        visibility_boost=max(0.1, float(args.visibility_boost)),
        priority_vmag=float(args.priority_vmag),
        priority_width=max(0.2, float(args.priority_width)),
        outline_vmag_max=float(args.outline_vmag_max),
    )
    app.quit()
    print(f"Saved comparison image to {args.outfile}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
