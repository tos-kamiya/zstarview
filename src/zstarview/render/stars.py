import math
from typing import Optional, Tuple

import hashlib
import logging
import numpy as np

from PySide6.QtCore import QPoint, QPointF
from PySide6.QtGui import (
    QImage,
    QPainter,
)

from ..types import ScreenGeometry, CelestialData, ViewerData, CelestialObject
from ..astro import altaz_to_normalized_xy, resolve_star_names, resolve_star_source_ids, is_in_fov, is_in_fov_vectorized
from ..asterisms import ASTERISM_REQUIRED_SOURCE_IDS
from .geometry import (
    _altaz_to_normalized_xy_vectorized,
    _normalized_to_screen_xy_vectorized,
    normalized_to_screen_xy,
)
from .photometry import bv_to_rgb_vectorized

logger = logging.getLogger(__name__)

_star_render_cache: tuple[tuple, QImage] | None = None
_MAG2_TO_MAG1_SIZE_SCALE = 10.0 ** 0.12
_DIAMOND_OVERLAY_GAIN = 0.85
_OUTLINE_DIAMOND_GAIN = 1.25
_OUTLINE_DIAMOND_SCALE = 0.72
_SINGLE_STAR_GAUSSIAN_STRENGTH = 0.12


def _content_fov_deg_from_viewer(viewer_data: ViewerData) -> float:
    return float(viewer_data.content_fov_deg)


def _array_hash(arr: np.ndarray) -> str:
    if arr.size == 0:
        return "empty"
    return hashlib.md5(arr.tobytes()).hexdigest()


def _apply_weak_gaussian3x3_rgb(arr: np.ndarray, strength: float) -> np.ndarray:
    """Apply a very weak 3x3 Gaussian only to reduce salt-like 1px star noise."""
    s = float(np.clip(strength, 0.0, 1.0))
    if s <= 0.0:
        return arr
    p = np.pad(arr, ((1, 1), (1, 1), (0, 0)), mode="constant")
    g = (
        p[:-2, :-2, :]
        + 2.0 * p[:-2, 1:-1, :]
        + p[:-2, 2:, :]
        + 2.0 * p[1:-1, :-2, :]
        + 4.0 * p[1:-1, 1:-1, :]
        + 2.0 * p[1:-1, 2:, :]
        + p[2:, :-2, :]
        + 2.0 * p[2:, 1:-1, :]
        + p[2:, 2:, :]
    ) / 16.0
    return arr * (1.0 - s) + g * s


def _add_rgb_pixel(arr: np.ndarray, x: int, y: int, color: np.ndarray) -> None:
    if x < 0 or y < 0:
        return
    h, w, _ = arr.shape
    if x >= w or y >= h:
        return
    arr[y, x, :] += color


def _draw_rgb_line(
    arr: np.ndarray,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    color: np.ndarray,
) -> None:
    x0_i = int(round(x0))
    y0_i = int(round(y0))
    x1_i = int(round(x1))
    y1_i = int(round(y1))

    dx = abs(x1_i - x0_i)
    dy = abs(y1_i - y0_i)
    sx = 1 if x0_i < x1_i else -1
    sy = 1 if y0_i < y1_i else -1
    err = dx - dy

    x = x0_i
    y = y0_i
    while True:
        _add_rgb_pixel(arr, x, y, color)
        if x == x1_i and y == y1_i:
            break
        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            x += sx
        if e2 < dx:
            err += dx
            y += sy


def _draw_diamond_outline_rgb(
    arr: np.ndarray,
    cx: float,
    cy: float,
    half_diag: float,
    color: np.ndarray,
) -> None:
    if half_diag <= 0.0:
        return
    top = (cx, cy - half_diag)
    right = (cx + half_diag, cy)
    bottom = (cx, cy + half_diag)
    left = (cx - half_diag, cy)
    _draw_rgb_line(arr, *top, *right, color)
    _draw_rgb_line(arr, *right, *bottom, color)
    _draw_rgb_line(arr, *bottom, *left, color)
    _draw_rgb_line(arr, *left, *top, color)


def _star_cache_key(
    alt: np.ndarray,
    az: np.ndarray,
    size_factor: np.ndarray,
    color_factor_base: np.ndarray,
    celestial_time_value: float,
    view_center: Tuple[float, float],
    geometry: ScreenGeometry,
    star_base_radius: float,
    visibility_boost: float,
    outline_bright_bodies: bool,
    outline_render_scale: float,
    draw_vmag_limit: float | None,
) -> tuple:
    return (
        _array_hash(alt),
        _array_hash(az),
        _array_hash(size_factor),
        _array_hash(color_factor_base),
        float(celestial_time_value),
        view_center[0],
        view_center[1],
        geometry.center,
        geometry.radius,
        float(star_base_radius),
        float(visibility_boost),
        bool(outline_bright_bodies),
        float(outline_render_scale),
        None if draw_vmag_limit is None else float(draw_vmag_limit),
    )


def find_highlighted_object(
    celestial_data: Optional[CelestialData],
    viewer_data: ViewerData,
    mouse_pos: QPoint,
    geometry: ScreenGeometry,
) -> Optional[Tuple[CelestialObject, QPointF]]:
    """
    Find the nearest celestial object to the mouse cursor.

    This function searches for the closest star or planet to the given mouse
    position. It performs a vectorized search for stars and a scalar search for
    planets to optimize performance.

    Args:
        celestial_data: A CelestialData object containing information about celestial objects.
        viewer_data: A ViewerData object containing the viewer's location and time.
        mouse_pos: The QPoint representing the current mouse cursor position.
        geometry: A ScreenGeometry object for coordinate transformations.

    Returns:
        A tuple containing the highlighted celestial object (as a dictionary or
        a CelestialObject) and its screen position (QPointF), or None if no
        object is close enough to the cursor.
    """
    min_dist_sq = 30**2  # squared pixels
    highlighted_object: Optional[Tuple[CelestialObject, QPointF]] = None

    def _has_named_star(name: object) -> bool:
        if name is None:
            return False
        return bool(str(name).strip())

    def _is_asterism_member(source_id: object) -> bool:
        if source_id is None:
            return False
        return str(source_id).strip() in ASTERISM_REQUIRED_SOURCE_IDS

    if not celestial_data:
        return None

    # Handle stars first (vectorized)
    stars = celestial_data.stars
    if stars["alt"].size > 0:
        names = resolve_star_names(stars, celestial_data.star_catalog_meta)
        source_ids = resolve_star_source_ids(stars, celestial_data.star_catalog_meta)
        hover_mask = np.array(
            [_has_named_star(name) or _is_asterism_member(source_id) for name, source_id in zip(names, source_ids)],
            dtype=bool,
        )
        if np.any(hover_mask):
            valid_indices = np.nonzero(hover_mask)[0]
            alt_named = stars["alt"][hover_mask]
            az_named = stars["az"][hover_mask]
            nx, ny = _altaz_to_normalized_xy_vectorized(
                alt_named,
                az_named,
                viewer_data.view_center,
                edge_fov_deg=float(viewer_data.edge_fov_deg),
            )
            x, y = _normalized_to_screen_xy_vectorized(nx, ny, geometry)
            dist_sq = (mouse_pos.x() - x) ** 2 + (mouse_pos.y() - y) ** 2
            closest_idx = np.argmin(dist_sq)
            if dist_sq[closest_idx] < min_dist_sq:
                min_dist_sq = dist_sq[closest_idx]
                original_idx = valid_indices[closest_idx]
                highlighted_star: CelestialObject = {key: val[original_idx] for key, val in stars.items()}
                highlighted_star["name"] = names[original_idx]
                highlighted_star["source_id"] = source_ids[original_idx]
                highlighted_object = (highlighted_star, QPointF(x[closest_idx], y[closest_idx]))

    content_fov_deg = _content_fov_deg_from_viewer(viewer_data)

    # Handle planets (scalar)
    for body in celestial_data.planets:
        # For planets, allow below-horizon display as long as they are in the current FOV.
        if not is_in_fov(body.alt, body.az, viewer_data.view_center, fov_deg=content_fov_deg):
            continue
        nx, ny = altaz_to_normalized_xy(
            body.alt,
            body.az,
            viewer_data.view_center,
            edge_fov_deg=float(viewer_data.edge_fov_deg),
        )
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        dist_sq = (mouse_pos.x() - px) ** 2 + (mouse_pos.y() - py) ** 2
        if dist_sq < min_dist_sq:
            min_dist_sq = dist_sq
            highlighted_object = (body, QPointF(px, py))  # body is already an object/dataclass

    return highlighted_object

def draw_stars(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    star_base_radius: float,
    *,
    visibility_boost: float = 1.0,
    outline_bright_bodies: bool = False,
    outline_render_scale: float = 1.0,
    draw_vmag_limit: Optional[float] = None,
    viewport_size: Tuple[int, int] | None = None,
    content_fov_deg: float | None = None,
) -> None:
    """
    Draw stars using a numpy canvas that paints uniformly sized rectangles.

    The visual area of each star is derived from its magnitude (a 1st-magnitude
    star maps to a 6x6 rectangle), and the color intensity is reduced for stars
    whose calculated area falls below 1 pixel^2. The resulting RGB canvas is
    converted into a `QImage` and blended additively to the painter before
    drawing other objects.

    Args:
        painter: The QPainter object for drawing.
        geometry: The screen geometry for coordinate conversions.
        celestial_data: The data containing star information (positions, magnitudes, etc.).
        viewer_data: The viewer's data, including the view center.
        star_base_radius: A base radius for scaling star sizes.
        visibility_boost: Multiplier for the overall star colors.
        outline_bright_bodies: When true, bright stars are rendered as outline-only diamonds.
        outline_render_scale: Scale factor to keep outline-only stars visually comparable to
            the usual downscaled star buffer.
        draw_vmag_limit: Optional override to skip fainter stars entirely.
        viewport_size: Optional `(width, height)` of the drawing area, used to create the numpy canvas; if omitted, defaults to the painter's clip rect.
    """
    stars = celestial_data.stars
    effective_fov_deg = _content_fov_deg_from_viewer(viewer_data) if content_fov_deg is None else float(content_fov_deg)
    visibility_boost = float(np.clip(visibility_boost, 0.7, 2.0))
    outline_render_scale = max(1.0, float(outline_render_scale))

    if draw_vmag_limit is not None:
        draw_mask = stars["vmag"] <= float(draw_vmag_limit)
        if not np.any(draw_mask):
            return
        alt = stars["alt"][draw_mask]
        az = stars["az"][draw_mask]
        vmag = stars["vmag"][draw_mask]
        bv = stars["bv"][draw_mask]
        size_factor = stars["size_factor"][draw_mask]
        color_factor_base = stars["color_factor_base"][draw_mask]
    else:
        alt = stars["alt"]
        az = stars["az"]
        vmag = stars["vmag"]
        bv = stars["bv"]
        size_factor = stars["size_factor"]
        color_factor_base = stars["color_factor_base"]

    width_px = None
    height_px = None
    if viewport_size:
        width_px = max(1, int(viewport_size[0]))
        height_px = max(1, int(viewport_size[1]))
    else:
        width_px = painter.viewport().width()
        height_px = painter.viewport().height()

    if width_px <= 0 or height_px <= 0:
        return

    # Convert directions to pixel centers.
    nx, ny = _altaz_to_normalized_xy_vectorized(
        alt,
        az,
        viewer_data.view_center,
        edge_fov_deg=float(viewer_data.edge_fov_deg),
    )
    x, y = _normalized_to_screen_xy_vectorized(nx, ny, geometry)

    bv_clamped = np.nan_to_num(bv, nan=0.45)
    rgb_colors = bv_to_rgb_vectorized(bv_clamped).astype(np.float32)

    # `star_base_radius` is defined as the apparent square size of a 2nd-magnitude star.
    max_size = max(12.0, float(max(1.0, star_base_radius)))
    size_float = float(star_base_radius) * _MAG2_TO_MAG1_SIZE_SCALE * size_factor
    if outline_bright_bodies:
        size_float *= outline_render_scale
    size_px = np.clip(np.round(size_float), 1, int(max_size)).astype(int)

    color_factor = np.clip(0.15 + 0.85 * color_factor_base, 0.0, 1.0) * visibility_boost
    color_factor = np.clip(color_factor, 0.0, 1.0)
    drawn_area = np.maximum(1.0, size_px.astype(float) ** 2)
    expected_area = size_float**2
    area_ratio = np.clip(expected_area / drawn_area, 0.0, 1.0)
    star_colors = rgb_colors * color_factor[:, None] * area_ratio[:, None]

    cache_key = _star_cache_key(
        alt,
        az,
        size_factor,
        color_factor_base,
        celestial_data.time.jd,
        viewer_data.view_center,
        geometry,
        star_base_radius,
        visibility_boost,
        outline_bright_bodies,
        outline_render_scale,
        draw_vmag_limit,
    )
    global _star_render_cache
    if _star_render_cache and _star_render_cache[0] == cache_key:
        painter.save()
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
        painter.drawImage(0, 0, _star_render_cache[1])
        painter.restore()
        return

    canvas = np.zeros((height_px, width_px, 3), dtype=np.float32)

    ix = np.round(x).astype(int)
    iy = np.round(y).astype(int)

    half = (size_px // 2).astype(int)
    x0 = ix - half
    y0 = iy - half
    x1 = x0 + size_px
    y1 = y0 + size_px

    x0_clamped = np.clip(x0, 0, width_px)
    y0_clamped = np.clip(y0, 0, height_px)
    x1_clamped = np.clip(x1, 0, width_px)
    y1_clamped = np.clip(y1, 0, height_px)

    outside_content = ~is_in_fov_vectorized(alt, az, viewer_data.view_center, fov_deg=effective_fov_deg)
    bright_outline_mask = outline_bright_bodies & (vmag < 2.0)
    valid_base = (
        (x1_clamped > x0_clamped)
        & (y1_clamped > y0_clamped)
        & (size_px > 0)
        & (~outside_content)
    )
    size_one = size_px == 1
    size_two = size_px == 2
    single_mask = valid_base & size_one & (~bright_outline_mask)
    size2_full_fit = (x0 >= 0) & (y0 >= 0) & (x1 <= width_px) & (y1 <= height_px)
    size2_mask = valid_base & size_two & size2_full_fit & (~bright_outline_mask)
    multi_mask = valid_base & ~(size_one | size_two) & (~bright_outline_mask)
    single_indices = np.nonzero(single_mask)[0]
    if single_indices.size > 0:
        single_layer = np.zeros_like(canvas)
        flat_single = single_layer.reshape(-1, 3)
        x_single = x0_clamped[single_indices]
        y_single = y0_clamped[single_indices]
        flat_idx = y_single * width_px + x_single
        np.add.at(flat_single, flat_idx, star_colors[single_indices])
        canvas += _apply_weak_gaussian3x3_rgb(single_layer, _SINGLE_STAR_GAUSSIAN_STRENGTH)

    size2_indices = np.nonzero(size2_mask)[0]
    if size2_indices.size > 0:
        size2_layer = np.zeros_like(canvas)
        flat_size2 = size2_layer.reshape(-1, 3)
        x_size2 = x0[size2_indices]
        y_size2 = y0[size2_indices]
        base_idx = y_size2 * width_px + x_size2
        colors_size2 = star_colors[size2_indices]
        np.add.at(flat_size2, base_idx, colors_size2)
        np.add.at(flat_size2, base_idx + 1, colors_size2)
        np.add.at(flat_size2, base_idx + width_px, colors_size2)
        np.add.at(flat_size2, base_idx + width_px + 1, colors_size2)
        canvas += size2_layer

    multi_indices = np.nonzero(multi_mask)[0]
    for idx in multi_indices:
        canvas[y0_clamped[idx]:y1_clamped[idx], x0_clamped[idx]:x1_clamped[idx], :] += star_colors[idx]

    # Overlay a rotated square (diamond) for bright stars to emphasize stars above 2nd magnitude.
    bright_indices = np.nonzero(valid_base & (vmag < 2.0))[0]
    if bright_indices.size > 0:
        bright_scale = np.power(10.0, 0.12 * np.clip(2.0 - vmag[bright_indices], 0.0, 4.0))
        for local_i, idx in enumerate(bright_indices):
            cx = float(ix[idx])
            cy = float(iy[idx])
            half_diag = max(
                0.5,
                0.5
                * float(size_px[idx])
                * float(bright_scale[local_i])
                * _OUTLINE_DIAMOND_SCALE,
            )
            if outline_bright_bodies:
                _draw_diamond_outline_rgb(
                    canvas,
                    cx,
                    cy,
                    half_diag,
                    star_colors[idx] * _OUTLINE_DIAMOND_GAIN,
                )
                continue
            xmin = max(0, int(math.floor(cx - half_diag)))
            xmax = min(width_px, int(math.ceil(cx + half_diag + 1.0)))
            ymin = max(0, int(math.floor(cy - half_diag)))
            ymax = min(height_px, int(math.ceil(cy + half_diag + 1.0)))
            if xmin >= xmax or ymin >= ymax:
                continue
            xs = np.arange(xmin, xmax, dtype=np.float32)[None, :]
            ys = np.arange(ymin, ymax, dtype=np.float32)[:, None]
            mask = (np.abs(xs - cx) + np.abs(ys - cy)) <= half_diag
            if not np.any(mask):
                continue
            canvas[ymin:ymax, xmin:xmax, :][mask] += star_colors[idx] * _DIAMOND_OVERLAY_GAIN

    np.clip(canvas, 0.0, 255.0, out=canvas)
    canvas_uint8 = np.ascontiguousarray(canvas.astype(np.uint8))
    if canvas_uint8.size == 0:
        return

    rgba = np.zeros((height_px, width_px, 4), dtype=np.uint8)
    alpha = np.max(canvas_uint8, axis=2)
    nonzero_alpha = alpha > 0
    if np.any(nonzero_alpha):
        rgb_float = canvas_uint8.astype(np.float32) / 255.0
        alpha_float = alpha.astype(np.float32) / 255.0
        rgba_rgb = np.zeros_like(rgb_float, dtype=np.float32)
        rgba_rgb[nonzero_alpha] = rgb_float[nonzero_alpha] / alpha_float[nonzero_alpha, None]
        np.clip(rgba_rgb, 0.0, 1.0, out=rgba_rgb)
        rgba[:, :, :3] = np.round(rgba_rgb * 255.0).astype(np.uint8)
    rgba[:, :, 3] = alpha
    image = QImage(rgba.data, width_px, height_px, width_px * 4, QImage.Format_RGBA8888).copy()
    painter.save()
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
    painter.drawImage(0, 0, image)
    logger.debug("Generated star buffer (%d stars) and cached image", len(size_px))
    _star_render_cache = (cache_key, image)

    painter.restore()
    # Reset composition mode.
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)
