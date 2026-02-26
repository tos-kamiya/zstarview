# -*- coding: utf-8 -*-
"""
Sky and cloud compositing utilities and cache.

This module provides:
- Pure helpers to generate hatch tiles and apply them to cloud alpha.
- A function to composite sky and cloud layers without relying on any global state.
- A small cache class that handles scaling, hatching, compositing, and reuse.
"""
from __future__ import annotations

import math
from typing import Optional, Tuple

import numpy as np
from PySide6.QtCore import QRect, Qt
from PySide6.QtGui import QImage, QPainter

from ..paths import CLOUD_HATCH_DEFAULT, HatchConfig
from ..types import ScreenGeometry
from ..utils.qt import np_rgba_to_qimage, qimage_to_np_rgba


def make_hatch_tile_qimage(W: int, H: int, line_px: int, strength: int) -> QImage:
    """Generate a diagonal hatch tile as a QImage (ARGB32_Premultiplied)."""
    norm = math.sqrt(W * W + H * H)
    period = max(1, int(round(norm * 0.5)))
    band_u = max(1, int(round(line_px)))

    xs = np.arange(W, dtype=np.int32)[None, :]
    ys = np.arange(H, dtype=np.int32)[:, None]
    # Enforce exact 45-degree hatch lines: x - y = const.
    u = xs - ys
    u_mod = np.mod(u, period)
    dist = np.minimum(u_mod, period - u_mod)
    mask = dist <= (band_u / 2)

    arr = np.zeros((H, W, 4), dtype=np.uint8)
    arr[..., 0:3] = 0
    arr[..., 3] = 0
    arr[..., 3][mask] = np.uint8(np.clip(strength, 0, 255))

    qimg = QImage(arr.tobytes(), W, H, QImage.Format_ARGB32_Premultiplied)
    return qimg.copy()


def cloud_with_hatched_alpha(cloud_img: QImage, hatch_cfg: HatchConfig) -> QImage:
    """Apply continuous 45-degree hatch directly to alpha (no tile seams)."""
    out = cloud_img if cloud_img.format() == QImage.Format_RGBA8888 else cloud_img.convertToFormat(QImage.Format_RGBA8888)
    cloud = qimage_to_np_rgba(out).astype(np.uint8)

    h, w = cloud.shape[:2]
    xs = np.arange(w, dtype=np.int32)[None, :]
    ys = np.arange(h, dtype=np.int32)[:, None]

    period = max(2, int(round(math.hypot(hatch_cfg.tile_w_px, hatch_cfg.tile_h_px) * 0.5)))
    band = max(1, int(round(hatch_cfg.line_px)))

    # Exact 45-degree lines: x - y = const, repeated by period.
    u = xs - ys
    u_mod = np.mod(u, period)
    dist = np.minimum(u_mod, period - u_mod)
    line_mask = dist <= (band / 2.0)

    erase = float(np.clip(hatch_cfg.strength, 0, 255)) / 255.0
    keep = 1.0 - erase

    # Match DestinationOut semantics: attenuate both alpha and premultiplied color.
    rgb = cloud[..., :3].astype(np.float32)
    alpha = cloud[..., 3].astype(np.float32)
    rgb[line_mask] = rgb[line_mask] * keep
    alpha[line_mask] = alpha[line_mask] * keep

    cloud[..., :3] = np.clip(np.round(rgb), 0, 255).astype(np.uint8)
    cloud[..., 3] = np.clip(np.round(alpha), 0, 255).astype(np.uint8)
    return np_rgba_to_qimage(cloud)


def compose_cloud_over_sky(
    sky_img: QImage,
    cloud_img_rgba: QImage,
    dest_rect: QRect,
    *,
    cloud_opacity: float = 1.0,
    gray_mix: float = 1.0,
) -> QImage:
    """Composite cloud over sky with optional gray desaturation behind clouds.

    - sky is blended with a grayscale version where cloud alpha is present
    - cloud color is additively applied with ``cloud_opacity``
    - final image is clipped to a circle
    """
    w, h = dest_rect.width(), dest_rect.height()

    if sky_img.width() != w or sky_img.height() != h:
        sky_img = sky_img.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)
    if cloud_img_rgba.width() != w or cloud_img_rgba.height() != h:
        cloud_img_rgba = cloud_img_rgba.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)

    sky_np = qimage_to_np_rgba(sky_img).astype(np.uint8)
    cloud_np = qimage_to_np_rgba(cloud_img_rgba).astype(np.uint8)

    r = sky_np[..., 0].astype(np.uint16)
    g = sky_np[..., 1].astype(np.uint16)
    b = sky_np[..., 2].astype(np.uint16)
    gray_u8 = ((77 * r + 150 * g + 29 * b) >> 8).astype(np.uint8)

    a = (cloud_np[..., 3].astype(np.float32) / 255.0) * float(np.clip(gray_mix, 0.0, 1.0))
    a8 = (np.clip(a, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint16)
    inv_a8 = 255 - a8

    sky_rgb_u16 = sky_np[..., :3].astype(np.uint16)
    gray_rgb_u16 = np.repeat(gray_u8[:, :, None], 3, axis=2).astype(np.uint16)
    base_u16 = (inv_a8[:, :, None] * sky_rgb_u16 + a8[:, :, None] * gray_rgb_u16) // 255

    cop = float(np.clip(cloud_opacity, 0.0, 1.0))
    if cop > 0.0:
        add_u16 = (cloud_np[..., :3].astype(np.uint16) * int(round(cop * 255))) // 255
        out_u16 = base_u16 + add_u16
        np.minimum(out_u16, 255, out=out_u16)
    else:
        out_u16 = base_u16

    out = np.empty((h, w, 4), dtype=np.uint8)
    out[..., :3] = out_u16.astype(np.uint8)
    out[..., 3] = 255

    cy, cx = (h - 1) * 0.5, (w - 1) * 0.5
    rr = min(cx, cy)
    y, x = np.arange(h, dtype=np.float32)[:, None], np.arange(w, dtype=np.float32)[None, :]
    r2 = (x - cx) ** 2 + (y - cy) ** 2
    mask = r2 <= (rr + 0.25) ** 2
    out[..., 3][~mask] = 0
    out[..., :3][~mask] = 0

    return np_rgba_to_qimage(out)


class SkyCompositorCache:
    """Manage compositing and reuse the last composited image via a cache key."""

    def __init__(
        self,
        *,
        hatch_cfg: HatchConfig = CLOUD_HATCH_DEFAULT,
        cloud_opacity_scale: float = 0.7,
        gray_mix: float = 1.0,
    ) -> None:
        self._hatch_cfg = hatch_cfg
        self._cloud_opacity_scale = cloud_opacity_scale
        self._gray_mix = gray_mix
        self._composited_img: Optional[QImage] = None
        self._composite_key: Optional[Tuple] = None

    def invalidate(self) -> None:
        self._composite_key = None
        self._composited_img = None

    def draw(
        self,
        painter: QPainter,
        geometry: ScreenGeometry,
        sky_img: Optional[QImage],
        cloud_img: Optional[QImage],
        *,
        cloud_alpha: float,
    ) -> None:
        """Composite the sky/cloud layers (with cache) and draw into painter."""
        if (sky_img is None) and (cloud_img is None or cloud_alpha <= 0.0):
            return

        x = int(geometry.center[0] - geometry.radius)
        y = int(geometry.center[1] - geometry.radius)
        w = h = int(geometry.radius * 2)

        sky_ck = int(sky_img.cacheKey()) if sky_img else 0
        cloud_ck = int(cloud_img.cacheKey()) if cloud_img else 0
        hatch_key = (
            self._hatch_cfg.tile_w_px,
            self._hatch_cfg.tile_h_px,
            self._hatch_cfg.line_px,
            self._hatch_cfg.strength,
        )
        comp_key = ("comp", sky_ck, cloud_ck, w, h, float(cloud_alpha), hatch_key, self._cloud_opacity_scale, self._gray_mix)

        if self._composite_key != comp_key or self._composited_img is None:
            def _scaled(qimg: Optional[QImage]) -> Optional[QImage]:
                if qimg is None:
                    return None
                if qimg.width() == w and qimg.height() == h:
                    return qimg
                return qimg.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)

            sky_s = _scaled(sky_img)
            cloud_s = _scaled(cloud_img)

            if cloud_s is not None and cloud_alpha > 0.0:
                cloud_s = cloud_with_hatched_alpha(cloud_s, self._hatch_cfg)

            if cloud_s is None or cloud_alpha <= 0.0:
                composited = sky_s if sky_s is not None else QImage(w, h, QImage.Format_ARGB32_Premultiplied)
                if composited is None or composited.isNull():
                    composited = QImage(w, h, QImage.Format_ARGB32_Premultiplied)
                    composited.fill(Qt.transparent)
            else:
                composited = compose_cloud_over_sky(
                    sky_img=sky_s if sky_s is not None else QImage(w, h, QImage.Format_ARGB32_Premultiplied),
                    cloud_img_rgba=cloud_s,
                    dest_rect=QRect(0, 0, w, h),
                    cloud_opacity=cloud_alpha * self._cloud_opacity_scale,
                    gray_mix=self._gray_mix,
                )

            self._composited_img = composited
            self._composite_key = comp_key

        painter.save()
        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        painter.drawImage(QRect(x, y, w, h), self._composited_img)
        painter.restore()
