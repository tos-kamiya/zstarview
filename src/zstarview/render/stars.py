import hashlib
import logging
import math

import numpy as np
from PySide6.QtCore import QPoint, QPointF, Qt
from PySide6.QtGui import (
    QColor,
    QImage,
    QPainter,
    QPen,
)

from ..asterisms import ASTERISM_REQUIRED_SOURCE_IDS
from ..astro import (
    altaz_to_normalized_xy,
    is_in_fov,
    is_in_fov_vectorized,
    resolve_star_names,
    resolve_star_source_ids,
)
from ..types import CelestialData, CelestialObject, ScreenGeometry, ViewerData
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
_DIAMOND_OVERLAY_SCALE = 0.72
_OUTLINE_DIAMOND_SCALE = 1.0
_LIGHT_BACKGROUND_BRIGHT_VMAG = 2.0
_LIGHT_BACKGROUND_OUTLINE_RGB = (24, 24, 24)
_LIGHT_BACKGROUND_OUTLINE_ALPHA = 85
_LIGHT_BACKGROUND_DARK_DIAMOND_THICKNESS = 3
_LIGHT_BACKGROUND_COLOR_DIAMOND_THICKNESS = 1
_LIGHT_BACKGROUND_COLOR_DIAMOND_ALPHA = 150
_SCENIC_DARK_UNDERLAY_ALPHA = 85
_SCENIC_DARK_UNDERLAY_WIDTH = 1.0


def _content_fov_deg_from_viewer(viewer_data: ViewerData) -> float:
    return float(viewer_data.content_fov_deg)


def _array_hash(arr: np.ndarray) -> str:
    if arr.size == 0:
        return "empty"
    return hashlib.md5(arr.tobytes()).hexdigest()


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


def _draw_rgba_line(
    arr: np.ndarray,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    color: np.ndarray,
    thickness: int = 1,
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
    height, width, _ = arr.shape
    radius = max(0, (int(thickness) - 1) // 2)
    while True:
        for offset_y in range(-radius, radius + 1):
            for offset_x in range(-radius, radius + 1):
                pixel_x = x + offset_x
                pixel_y = y + offset_y
                if 0 <= pixel_x < width and 0 <= pixel_y < height:
                    arr[pixel_y, pixel_x, :] = color
        if x == x1_i and y == y1_i:
            break
        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            x += sx
        if e2 < dx:
            err += dx
            y += sy


def _draw_diamond_outline_rgba(
    arr: np.ndarray,
    cx: float,
    cy: float,
    half_diag: float,
    color: np.ndarray,
    thickness: int = 1,
) -> None:
    if half_diag <= 0.0:
        return
    top = (cx, cy - half_diag)
    right = (cx + half_diag, cy)
    bottom = (cx, cy + half_diag)
    left = (cx - half_diag, cy)
    _draw_rgba_line(arr, *top, *right, color, thickness=thickness)
    _draw_rgba_line(arr, *right, *bottom, color, thickness=thickness)
    _draw_rgba_line(arr, *bottom, *left, color, thickness=thickness)
    _draw_rgba_line(arr, *left, *top, color, thickness=thickness)


def _draw_diamond_fill_rgba(
    arr: np.ndarray,
    cx: float,
    cy: float,
    half_diag: float,
    color: np.ndarray,
) -> None:
    if half_diag <= 0.0:
        return
    height, width, _ = arr.shape
    xmin = max(0, int(math.floor(cx - half_diag)))
    xmax = min(width, int(math.ceil(cx + half_diag + 1.0)))
    ymin = max(0, int(math.floor(cy - half_diag)))
    ymax = min(height, int(math.ceil(cy + half_diag + 1.0)))
    if xmin >= xmax or ymin >= ymax:
        return
    xs = np.arange(xmin, xmax, dtype=np.float32)[None, :]
    ys = np.arange(ymin, ymax, dtype=np.float32)[:, None]
    mask = (np.abs(xs - cx) + np.abs(ys - cy)) <= half_diag
    if np.any(mask):
        arr[ymin:ymax, xmin:xmax, :][mask] = color


def _draw_rect_outline_rgb(
    arr: np.ndarray,
    x0: int,
    y0: int,
    x1: int,
    y1: int,
    color: np.ndarray,
    thickness: int = 1,
) -> None:
    if x1 <= x0 or y1 <= y0:
        return
    t = max(1, int(thickness))
    x1_i = min(x1, x0 + t)
    y1_i = min(y1, y0 + t)
    x0_i = max(x0, x1 - t)
    y0_i = max(y0, y1 - t)
    arr[y0:y1_i, x0:x1, :] += color
    if y0_i > y0:
        arr[y0_i:y1, x0:x1, :] += color
    if x1_i > x0:
        arr[y0:y1, x0:x1_i, :] += color
    if x0_i > x0:
        arr[y0:y1, x0_i:x1, :] += color


def _draw_rect_outline_batch_rgb(
    arr: np.ndarray,
    x0: np.ndarray,
    y0: np.ndarray,
    size_px: int,
    colors: np.ndarray,
    width_px: int,
    thickness: int = 1,
) -> None:
    """Draw same-sized rectangular outlines in batches for better performance."""
    if x0.size == 0:
        return
    if size_px <= 0:
        return

    flat = arr.reshape(-1, 3)
    t = max(1, int(thickness))
    band = min(size_px, t)
    x_offsets = np.arange(size_px, dtype=int)
    color_rows = np.repeat(colors, band * size_px, axis=0)

    x_band = np.repeat(x0, band)
    top_rows = np.repeat(y0, band) + np.tile(np.arange(band, dtype=int), x0.size)
    top_idx = (top_rows[:, None] * width_px + x_band[:, None] + x_offsets[None, :]).reshape(-1)
    np.add.at(flat, top_idx, color_rows)

    bottom_y0 = y0 + size_px - band
    bottom_rows = np.repeat(bottom_y0, band) + np.tile(np.arange(band, dtype=int), x0.size)
    bottom_idx = (bottom_rows[:, None] * width_px + x_band[:, None] + x_offsets[None, :]).reshape(-1)
    np.add.at(flat, bottom_idx, color_rows)

    if size_px <= band * 2:
        return

    inner_y_offsets = np.arange(band, size_px - band, dtype=int)
    side_colors = np.repeat(colors, inner_y_offsets.size, axis=0)

    y_band = np.repeat(y0, inner_y_offsets.size) + np.tile(inner_y_offsets, x0.size)
    left_idx = (y_band * width_px + np.repeat(x0, inner_y_offsets.size)).reshape(-1)
    right_x = x0 + size_px - 1
    right_idx = (y_band * width_px + np.repeat(right_x, inner_y_offsets.size)).reshape(-1)
    np.add.at(flat, left_idx, side_colors)
    np.add.at(flat, right_idx, side_colors)


def _star_cache_key(
    alt: np.ndarray,
    az: np.ndarray,
    size_factor: np.ndarray,
    color_factor_base: np.ndarray,
    celestial_time_value: float,
    view_center: tuple[float, float],
    geometry: ScreenGeometry,
    star_base_radius: float,
    visibility_boost: float,
    outline_bright_bodies: bool,
    outline_render_scale: float,
    draw_vmag_limit: float | None,
    draw_vmag_min_exclusive: float | None,
    fast_mode: bool,
    light_background_outline: bool = False,
    twinkle_targets: tuple[tuple[int, float], ...] = (),
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
        None if draw_vmag_min_exclusive is None else float(draw_vmag_min_exclusive),
        bool(fast_mode),
        bool(light_background_outline),
        tuple(
            (int(star_index), round(float(alpha), 6))
            for star_index, alpha in twinkle_targets
        ),
    )


def _draw_stars_light_background_rgba(
    painter: QPainter,
    *,
    x: np.ndarray,
    y: np.ndarray,
    alt: np.ndarray,
    az: np.ndarray,
    viewer_data: ViewerData,
    view_center: tuple[float, float],
    effective_fov_deg: float,
    size_float: np.ndarray,
    rgb_colors: np.ndarray,
    vmag: np.ndarray,
    outline_bright_bodies: bool,
    width_px: int,
    height_px: int,
    geometry: ScreenGeometry,
    celestial_time_value: float,
    star_base_radius: float,
    draw_vmag_limit: float | None,
    fast_mode: bool,
) -> None:
    apparent_diameter_px = size_float.astype(float, copy=False)
    visible_mask = (
        (apparent_diameter_px >= 1.0)
        & is_in_fov_vectorized(alt, az, view_center, fov_deg=effective_fov_deg)
    )
    if not np.any(visible_mask):
        return

    x = x[visible_mask]
    y = y[visible_mask]
    size_float = apparent_diameter_px[visible_mask]
    rgb_colors = rgb_colors[visible_mask]
    vmag = vmag[visible_mask]

    body_size = np.maximum(1, np.round(size_float).astype(int))
    outline_size = np.maximum(body_size + 2, np.ceil(size_float + 2.0).astype(int))
    ix = np.round(x).astype(int)
    iy = np.round(y).astype(int)

    cache_key = (
        "light-background",
        _array_hash(x.astype(np.float32, copy=False)),
        _array_hash(y.astype(np.float32, copy=False)),
        _array_hash(size_float.astype(np.float32, copy=False)),
        _array_hash(rgb_colors.astype(np.float32, copy=False)),
        _array_hash(vmag.astype(np.float32, copy=False)),
        bool(outline_bright_bodies),
        float(_LIGHT_BACKGROUND_BRIGHT_VMAG),
        float(_DIAMOND_OVERLAY_SCALE),
        float(_OUTLINE_DIAMOND_SCALE),
        int(_LIGHT_BACKGROUND_DARK_DIAMOND_THICKNESS),
        int(_LIGHT_BACKGROUND_COLOR_DIAMOND_THICKNESS),
        int(_LIGHT_BACKGROUND_COLOR_DIAMOND_ALPHA),
        float(celestial_time_value),
        viewer_data.view_center[0],
        viewer_data.view_center[1],
        geometry.center,
        geometry.radius,
        float(star_base_radius),
        None if draw_vmag_limit is None else float(draw_vmag_limit),
        bool(fast_mode),
        int(width_px),
        int(height_px),
    )
    global _star_render_cache
    if _star_render_cache and _star_render_cache[0] == cache_key:
        painter.save()
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)
        painter.drawImage(0, 0, _star_render_cache[1])
        painter.restore()
        return

    rgba = np.zeros((height_px, width_px, 4), dtype=np.uint8)

    bright_mask = vmag <= _LIGHT_BACKGROUND_BRIGHT_VMAG
    regular_mask = ~bright_mask

    outline_colors = np.zeros((outline_size.size, 4), dtype=np.uint8)
    outline_colors[:, :3] = _LIGHT_BACKGROUND_OUTLINE_RGB
    outline_colors[:, 3] = _LIGHT_BACKGROUND_OUTLINE_ALPHA
    if np.any(regular_mask):
        regular_outline_colors = outline_colors[regular_mask]
        regular_outline_sizes = outline_size[regular_mask]
        regular_ix = ix[regular_mask]
        regular_iy = iy[regular_mask]
        for idx in range(regular_outline_sizes.size):
            half = int(regular_outline_sizes[idx] // 2)
            size = int(regular_outline_sizes[idx])
            raw_x0 = int(regular_ix[idx] - half)
            raw_y0 = int(regular_iy[idx] - half)
            x0 = max(0, raw_x0)
            y0 = max(0, raw_y0)
            x1 = min(width_px, raw_x0 + size)
            y1 = min(height_px, raw_y0 + size)
            if x1 > x0 and y1 > y0:
                rgba[y0:y1, x0:x1, :] = regular_outline_colors[idx]

    body_colors = np.zeros((body_size.size, 4), dtype=np.uint8)
    body_colors[:, :3] = np.clip(np.round(rgb_colors), 0, 255).astype(np.uint8)
    body_colors[:, 3] = 255
    if np.any(regular_mask):
        regular_body_sizes = body_size[regular_mask]
        regular_body_colors = body_colors[regular_mask]
        regular_ix = ix[regular_mask]
        regular_iy = iy[regular_mask]
        for idx in range(regular_body_sizes.size):
            half = int(regular_body_sizes[idx] // 2)
            size = int(regular_body_sizes[idx])
            raw_x0 = int(regular_ix[idx] - half)
            raw_y0 = int(regular_iy[idx] - half)
            x0 = max(0, raw_x0)
            y0 = max(0, raw_y0)
            x1 = min(width_px, raw_x0 + size)
            y1 = min(height_px, raw_y0 + size)
            if x1 > x0 and y1 > y0:
                rgba[y0:y1, x0:x1, :] = regular_body_colors[idx]

    bright_indices = np.nonzero(bright_mask)[0]
    if bright_indices.size > 0:
        bright_scale = np.power(
            10.0,
            0.12 * np.clip(_LIGHT_BACKGROUND_BRIGHT_VMAG - vmag[bright_indices], 0.0, 4.0),
        )
        dark_color = np.array(
            (*_LIGHT_BACKGROUND_OUTLINE_RGB, _LIGHT_BACKGROUND_OUTLINE_ALPHA),
            dtype=np.uint8,
        )
        for local_i, idx in enumerate(bright_indices):
            cx = float(ix[idx])
            cy = float(iy[idx])
            marker_color = body_colors[idx].copy()
            marker_color[3] = _LIGHT_BACKGROUND_COLOR_DIAMOND_ALPHA
            diamond_scale = (
                _OUTLINE_DIAMOND_SCALE
                if outline_bright_bodies
                else _DIAMOND_OVERLAY_SCALE
            )
            half_diag = max(
                0.5,
                0.5 * float(body_size[idx]) * float(bright_scale[local_i]) * diamond_scale,
            )
            if outline_bright_bodies:
                _draw_diamond_outline_rgba(
                    rgba,
                    cx,
                    cy,
                    half_diag,
                    dark_color,
                    thickness=_LIGHT_BACKGROUND_DARK_DIAMOND_THICKNESS,
                )
                _draw_diamond_outline_rgba(
                    rgba,
                    cx,
                    cy,
                    half_diag,
                    marker_color,
                    thickness=_LIGHT_BACKGROUND_COLOR_DIAMOND_THICKNESS,
                )
            else:
                _draw_diamond_fill_rgba(
                    rgba,
                    cx,
                    cy,
                    half_diag + 1.0,
                    dark_color,
                )
                _draw_diamond_fill_rgba(
                    rgba,
                    cx,
                    cy,
                    half_diag,
                    marker_color,
                )

    image = QImage(rgba.data, width_px, height_px, width_px * 4, QImage.Format_RGBA8888).copy()
    painter.save()
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)
    painter.drawImage(0, 0, image)
    painter.restore()
    _star_render_cache = (cache_key, image)


def find_highlighted_object(
    celestial_data: CelestialData | None,
    viewer_data: ViewerData,
    mouse_pos: QPoint,
    geometry: ScreenGeometry,
) -> tuple[CelestialObject, QPointF] | None:
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
    highlighted_object: tuple[CelestialObject, QPointF] | None = None

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


def draw_bright_star_underlay(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    star_base_radius: float,
    *,
    outline_bright_bodies: bool = False,
    outline_render_scale: float = 1.0,
    draw_vmag_limit: float = 4.0,
    viewport_size: tuple[int, int] | None = None,
    content_fov_deg: float | None = None,
) -> None:
    """Draw the local dark backing for bright stars before their colored bodies."""
    stars = celestial_data.stars
    vmag_values = np.asarray(stars["vmag"], dtype=float)
    mask = vmag_values <= float(draw_vmag_limit)
    if not np.any(mask):
        return

    alt = stars["alt"][mask]
    az = stars["az"][mask]
    vmag = stars["vmag"][mask]
    size_factor = stars["size_factor"][mask]
    effective_fov_deg = (
        _content_fov_deg_from_viewer(viewer_data)
        if content_fov_deg is None
        else float(content_fov_deg)
    )
    nx, ny = _altaz_to_normalized_xy_vectorized(
        alt,
        az,
        viewer_data.view_center,
        edge_fov_deg=float(viewer_data.edge_fov_deg),
    )
    x, y = _normalized_to_screen_xy_vectorized(nx, ny, geometry)
    size_float = float(star_base_radius) * _MAG2_TO_MAG1_SIZE_SCALE * size_factor
    if outline_bright_bodies:
        size_float *= max(1.0, float(outline_render_scale))
    visible = (
        (size_float >= 1.0)
        & is_in_fov_vectorized(alt, az, viewer_data.view_center, fov_deg=effective_fov_deg)
    )
    if not np.any(visible):
        return

    x = x[visible]
    y = y[visible]
    vmag = vmag[visible]
    size_float = size_float[visible]
    dark = QColor(0, 0, 0, _SCENIC_DARK_UNDERLAY_ALPHA)
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
    pen = QPen(dark)
    pen.setWidthF(_SCENIC_DARK_UNDERLAY_WIDTH)
    painter.setPen(pen)
    painter.setBrush(dark)
    for cx, cy, size, magnitude in zip(x, y, size_float, vmag):
        ix = int(round(float(cx)))
        iy = int(round(float(cy)))
        body_size = max(1.0, float(round(size)))
        if magnitude <= 2.0:
            bright_scale = 10.0 ** (0.12 * np.clip(2.0 - float(magnitude), 0.0, 4.0))
            half_diag = max(1.0, 0.5 * body_size * bright_scale)
            points = [
                (ix, int(round(iy - half_diag))),
                (int(round(ix + half_diag)), iy),
                (ix, int(round(iy + half_diag))),
                (int(round(ix - half_diag)), iy),
            ]
            polygon = [QPoint(px, py) for px, py in points]
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPolygon(polygon)
        else:
            half = max(1, int(round(body_size / 2.0)))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(
                ix - half - 1,
                iy - half - 1,
                int(body_size) + 1,
                int(body_size) + 1,
            )
    painter.restore()


def collect_visible_named_star_labels(
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    geometry: ScreenGeometry,
    star_base_radius: float,
    *,
    outline_bright_bodies: bool = False,
    outline_render_scale: float = 1.0,
    draw_vmag_limit: float | None = None,
    content_fov_deg: float | None = None,
    viewport_size: tuple[int, int] | None = None,
) -> list[tuple[str, QPointF, tuple[int, int, int]]]:
    """Return screen positions for named stars that are currently drawn."""
    stars = celestial_data.stars
    if stars["alt"].size == 0:
        return []

    if draw_vmag_limit is not None:
        draw_mask = stars["vmag"] <= float(draw_vmag_limit)
        if not np.any(draw_mask):
            return []
        alt = stars["alt"][draw_mask]
        az = stars["az"][draw_mask]
        size_factor = stars["size_factor"][draw_mask]
        bv = stars["bv"][draw_mask]
    else:
        alt = stars["alt"]
        az = stars["az"]
        size_factor = stars["size_factor"]
        bv = stars["bv"]

    names = np.asarray(resolve_star_names(stars, celestial_data.star_catalog_meta), dtype=object)
    if draw_vmag_limit is not None:
        names = names[draw_mask]

    if names.size == 0:
        return []

    effective_fov_deg = (
        _content_fov_deg_from_viewer(viewer_data)
        if content_fov_deg is None
        else float(content_fov_deg)
    )
    nx, ny = _altaz_to_normalized_xy_vectorized(
        alt,
        az,
        viewer_data.view_center,
        edge_fov_deg=float(viewer_data.edge_fov_deg),
    )
    x, y = _normalized_to_screen_xy_vectorized(nx, ny, geometry)
    bv_clamped = np.nan_to_num(bv, nan=0.45)
    rgb_colors = bv_to_rgb_vectorized(bv_clamped).astype(np.float32)

    max_size = max(12.0, float(max(1.0, star_base_radius)))
    size_float = float(star_base_radius) * _MAG2_TO_MAG1_SIZE_SCALE * size_factor
    if outline_bright_bodies:
        size_float *= float(outline_render_scale)
    size_px = np.clip(np.round(size_float), 1, int(max_size)).astype(int)
    star_colors = rgb_colors

    if viewport_size is None:
        width_px = max(1, int(geometry.center[0] * 2))
        height_px = max(1, int(geometry.center[1] * 2))
    else:
        width_px = max(1, int(viewport_size[0]))
        height_px = max(1, int(viewport_size[1]))
    x0 = np.round(x).astype(int) - (size_px // 2).astype(int)
    y0 = np.round(y).astype(int) - (size_px // 2).astype(int)
    x1 = x0 + size_px
    y1 = y0 + size_px
    x0_clamped = np.clip(x0, 0, width_px)
    y0_clamped = np.clip(y0, 0, height_px)
    x1_clamped = np.clip(x1, 0, width_px)
    y1_clamped = np.clip(y1, 0, height_px)
    outside_content = ~is_in_fov_vectorized(
        alt,
        az,
        viewer_data.view_center,
        fov_deg=effective_fov_deg,
    )
    valid_base = (
        (x1_clamped > x0_clamped)
        & (y1_clamped > y0_clamped)
        & (size_px > 0)
        & (~outside_content)
    )
    label_mask = valid_base & np.array(
        [bool(str(name).strip()) for name in names],
        dtype=bool,
    )
    if not np.any(label_mask):
        return []

    positions: list[tuple[str, QPointF, tuple[int, int, int]]] = []
    label_indices = np.nonzero(label_mask)[0]
    for idx in label_indices:
        color = tuple(int(v) for v in np.clip(np.round(star_colors[idx]), 0, 255))
        positions.append(
            (str(names[idx]).strip(), QPointF(float(x[idx]), float(y[idx])), color)
        )
    return positions

def _draw_stars_render(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    star_base_radius: float,
    *,
    visibility_boost: float = 1.0,
    outline_bright_bodies: bool = False,
    outline_render_scale: float = 1.0,
    light_background_outline: bool = False,
    draw_vmag_limit: float | None = None,
    draw_vmag_min_exclusive: float | None = None,
    fast_mode: bool = False,
    viewport_size: tuple[int, int] | None = None,
    content_fov_deg: float | None = None,
    twinkle_targets: tuple[tuple[int, float], ...] = (),
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
    star_time = celestial_data.star_time or celestial_data.time
    effective_fov_deg = _content_fov_deg_from_viewer(viewer_data) if content_fov_deg is None else float(content_fov_deg)
    visibility_boost = float(np.clip(visibility_boost, 0.7, 2.0))
    outline_render_scale = max(1.0, float(outline_render_scale))
    fast_mode = bool(fast_mode)

    if draw_vmag_limit is not None:
        draw_mask = stars["vmag"] <= float(draw_vmag_limit)
        if draw_vmag_min_exclusive is not None:
            draw_mask &= stars["vmag"] > float(draw_vmag_min_exclusive)
        if not np.any(draw_mask):
            return
        source_index = stars["star_index"][draw_mask]
        alt = stars["alt"][draw_mask]
        az = stars["az"][draw_mask]
        vmag = stars["vmag"][draw_mask]
        bv = stars["bv"][draw_mask]
        size_factor = stars["size_factor"][draw_mask]
        color_factor_base = stars["color_factor_base"][draw_mask]
    else:
        source_index = stars["star_index"]
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
    if light_background_outline:
        _draw_stars_light_background_rgba(
            painter,
            x=x,
            y=y,
            alt=alt,
            az=az,
            viewer_data=viewer_data,
            view_center=viewer_data.view_center,
            effective_fov_deg=effective_fov_deg,
            size_float=size_float,
            rgb_colors=rgb_colors,
            vmag=vmag,
            outline_bright_bodies=outline_bright_bodies,
            width_px=width_px,
            height_px=height_px,
            geometry=geometry,
            celestial_time_value=star_time.jd,
            star_base_radius=star_base_radius,
            draw_vmag_limit=draw_vmag_limit,
            fast_mode=fast_mode,
        )
        return
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
        star_time.jd,
        viewer_data.view_center,
        geometry,
        star_base_radius,
        visibility_boost,
        outline_bright_bodies,
        outline_render_scale,
        draw_vmag_limit,
        draw_vmag_min_exclusive,
        fast_mode,
        light_background_outline=light_background_outline,
        twinkle_targets=twinkle_targets,
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
    bright_outline_mask = outline_bright_bodies & (vmag <= 2.0)
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
    if fast_mode:
        size3_to_6_mask = valid_base & (size_px >= 3) & (~bright_outline_mask)
        outline_mask = np.zeros_like(valid_base, dtype=bool)
    else:
        size3_to_6_mask = valid_base & (size_px >= 3) & (size_px <= 6) & (~bright_outline_mask)
        outline_mask = valid_base & (size_px >= 7) & (~bright_outline_mask)
    single_indices = np.nonzero(single_mask)[0]
    if single_indices.size > 0:
        single_layer = np.zeros_like(canvas)
        flat_single = single_layer.reshape(-1, 3)
        x_single = x0_clamped[single_indices]
        y_single = y0_clamped[single_indices]
        flat_idx = y_single * width_px + x_single
        np.add.at(flat_single, flat_idx, star_colors[single_indices])
        canvas += single_layer

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

    size3_to_6_indices = np.nonzero(size3_to_6_mask)[0]
    if size3_to_6_indices.size > 0:
        for idx in size3_to_6_indices:
            canvas[
                y0_clamped[idx]:y1_clamped[idx],
                x0_clamped[idx]:x1_clamped[idx],
                :,
            ] += star_colors[idx]

    outline_indices = np.nonzero(outline_mask)[0]
    if outline_indices.size > 0:
        outline_sizes = size_px[outline_indices]
        outline_x0 = x0[outline_indices]
        outline_y0 = y0[outline_indices]
        outline_x1 = x1[outline_indices]
        outline_y1 = y1[outline_indices]
        outline_colors = star_colors[outline_indices]

        outline_full_fit = (
            (outline_x0 >= 0)
            & (outline_y0 >= 0)
            & (outline_x1 <= width_px)
            & (outline_y1 <= height_px)
        )
        if np.any(outline_full_fit):
            full_fit_indices = np.nonzero(outline_full_fit)[0]
            for size_value in np.unique(outline_sizes[full_fit_indices]):
                size_mask = outline_full_fit & (outline_sizes == size_value)
                size_indices = np.nonzero(size_mask)[0]
                _draw_rect_outline_batch_rgb(
                    canvas,
                    outline_x0[size_indices],
                    outline_y0[size_indices],
                    int(size_value),
                    outline_colors[size_indices],
                    width_px,
                    thickness=2,
                )

        clipped_indices = np.nonzero(~outline_full_fit)[0]
        for idx in clipped_indices:
            _draw_rect_outline_rgb(
                canvas,
                int(x0_clamped[outline_indices[idx]]),
                int(y0_clamped[outline_indices[idx]]),
                int(x1_clamped[outline_indices[idx]]),
                int(y1_clamped[outline_indices[idx]]),
                outline_colors[idx],
                thickness=2,
            )

    # Overlay a rotated square (diamond) for bright stars to emphasize stars above 2nd magnitude.
    bright_indices = np.nonzero(valid_base & (vmag <= 2.0))[0]
    if bright_indices.size > 0:
        bright_scale = np.power(10.0, 0.12 * np.clip(2.0 - vmag[bright_indices], 0.0, 4.0))
        for local_i, idx in enumerate(bright_indices):
            cx = float(ix[idx])
            cy = float(iy[idx])
            diamond_scale = (
                _OUTLINE_DIAMOND_SCALE if outline_bright_bodies else _DIAMOND_OVERLAY_SCALE
            )
            half_diag = max(
                0.5,
                0.5 * float(size_px[idx]) * float(bright_scale[local_i]) * diamond_scale,
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

    for target_star_index, target_alpha in twinkle_targets:
        if float(target_alpha) <= 0.0:
            continue
        selected_indices = np.flatnonzero(source_index == int(target_star_index))
        if selected_indices.size > 0:
            selected = int(selected_indices[0])
            radius = max(
                0.5,
                0.5 * float(size_px[selected]),
            )
            cx = float(x[selected])
            cy = float(y[selected])
            xmin = max(0, int(math.floor(cx - radius)))
            xmax = min(width_px, int(math.ceil(cx + radius + 1.0)))
            ymin = max(0, int(math.floor(cy - radius)))
            ymax = min(height_px, int(math.ceil(cy + radius + 1.0)))
            if xmin < xmax and ymin < ymax:
                xs = np.arange(xmin, xmax, dtype=np.float32)[None, :]
                ys = np.arange(ymin, ymax, dtype=np.float32)[:, None]
                mask = (xs - cx) ** 2 + (ys - cy) ** 2 <= radius**2
                canvas[ymin:ymax, xmin:xmax, :][mask] *= max(
                    0.0, 1.0 - float(np.clip(target_alpha, 0.0, 1.0))
                )

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


def _draw_stars_fast_impl(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    star_base_radius: float,
    *,
    visibility_boost: float = 1.0,
    outline_bright_bodies: bool = False,
    outline_render_scale: float = 1.0,
    light_background_outline: bool = False,
    draw_vmag_limit: float | None = None,
    draw_vmag_min_exclusive: float | None = None,
    viewport_size: tuple[int, int] | None = None,
    content_fov_deg: float | None = None,
    twinkle_targets: tuple[tuple[int, float], ...] = (),
) -> None:
    _draw_stars_render(
        painter,
        geometry,
        celestial_data,
        viewer_data,
        star_base_radius,
        visibility_boost=visibility_boost,
        outline_bright_bodies=outline_bright_bodies,
        outline_render_scale=outline_render_scale,
        light_background_outline=light_background_outline,
        draw_vmag_limit=draw_vmag_limit,
        draw_vmag_min_exclusive=draw_vmag_min_exclusive,
        fast_mode=True,
        viewport_size=viewport_size,
        content_fov_deg=content_fov_deg,
        twinkle_targets=twinkle_targets,
    )


def _draw_stars_normal_impl(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    star_base_radius: float,
    *,
    visibility_boost: float = 1.0,
    outline_bright_bodies: bool = False,
    outline_render_scale: float = 1.0,
    light_background_outline: bool = False,
    draw_vmag_limit: float | None = None,
    draw_vmag_min_exclusive: float | None = None,
    viewport_size: tuple[int, int] | None = None,
    content_fov_deg: float | None = None,
    twinkle_targets: tuple[tuple[int, float], ...] = (),
) -> None:
    _draw_stars_render(
        painter,
        geometry,
        celestial_data,
        viewer_data,
        star_base_radius,
        visibility_boost=visibility_boost,
        outline_bright_bodies=outline_bright_bodies,
        outline_render_scale=outline_render_scale,
        light_background_outline=light_background_outline,
        draw_vmag_limit=draw_vmag_limit,
        draw_vmag_min_exclusive=draw_vmag_min_exclusive,
        fast_mode=False,
        viewport_size=viewport_size,
        content_fov_deg=content_fov_deg,
        twinkle_targets=twinkle_targets,
    )


def draw_stars_fast(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    star_base_radius: float,
    *,
    visibility_boost: float = 1.0,
    outline_bright_bodies: bool = False,
    outline_render_scale: float = 1.0,
    light_background_outline: bool = False,
    draw_vmag_limit: float | None = None,
    draw_vmag_min_exclusive: float | None = None,
    viewport_size: tuple[int, int] | None = None,
    content_fov_deg: float | None = None,
    twinkle_targets: tuple[tuple[int, float], ...] = (),
) -> None:
    """Draw stars using the fast-mode star simplifications."""
    _draw_stars_fast_impl(
        painter,
        geometry,
        celestial_data,
        viewer_data,
        star_base_radius,
        visibility_boost=visibility_boost,
        outline_bright_bodies=outline_bright_bodies,
        outline_render_scale=outline_render_scale,
        light_background_outline=light_background_outline,
        draw_vmag_limit=draw_vmag_limit,
        draw_vmag_min_exclusive=draw_vmag_min_exclusive,
        viewport_size=viewport_size,
        content_fov_deg=content_fov_deg,
        twinkle_targets=twinkle_targets,
    )


def draw_stars_normal(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    star_base_radius: float,
    *,
    visibility_boost: float = 1.0,
    outline_bright_bodies: bool = False,
    outline_render_scale: float = 1.0,
    light_background_outline: bool = False,
    draw_vmag_limit: float | None = None,
    draw_vmag_min_exclusive: float | None = None,
    viewport_size: tuple[int, int] | None = None,
    content_fov_deg: float | None = None,
    twinkle_targets: tuple[tuple[int, float], ...] = (),
) -> None:
    """Draw stars using the full normal-mode star renderer."""
    _draw_stars_normal_impl(
        painter,
        geometry,
        celestial_data,
        viewer_data,
        star_base_radius,
        visibility_boost=visibility_boost,
        outline_bright_bodies=outline_bright_bodies,
        outline_render_scale=outline_render_scale,
        light_background_outline=light_background_outline,
        draw_vmag_limit=draw_vmag_limit,
        draw_vmag_min_exclusive=draw_vmag_min_exclusive,
        viewport_size=viewport_size,
        content_fov_deg=content_fov_deg,
        twinkle_targets=twinkle_targets,
    )


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
    light_background_outline: bool = False,
    draw_vmag_limit: float | None = None,
    draw_vmag_min_exclusive: float | None = None,
    fast_mode: bool = False,
    viewport_size: tuple[int, int] | None = None,
    content_fov_deg: float | None = None,
    twinkle_targets: tuple[tuple[int, float], ...] = (),
) -> None:
    """Compatibility wrapper kept for existing callers and tests."""
    if fast_mode:
        _draw_stars_fast_impl(
            painter,
            geometry,
            celestial_data,
            viewer_data,
            star_base_radius,
            visibility_boost=visibility_boost,
            outline_bright_bodies=outline_bright_bodies,
            outline_render_scale=outline_render_scale,
            light_background_outline=light_background_outline,
            draw_vmag_limit=draw_vmag_limit,
            draw_vmag_min_exclusive=draw_vmag_min_exclusive,
            viewport_size=viewport_size,
            content_fov_deg=content_fov_deg,
            twinkle_targets=twinkle_targets,
        )
        return
    _draw_stars_normal_impl(
        painter,
        geometry,
        celestial_data,
        viewer_data,
        star_base_radius,
        visibility_boost=visibility_boost,
        outline_bright_bodies=outline_bright_bodies,
        outline_render_scale=outline_render_scale,
        light_background_outline=light_background_outline,
        draw_vmag_limit=draw_vmag_limit,
        draw_vmag_min_exclusive=draw_vmag_min_exclusive,
        viewport_size=viewport_size,
        content_fov_deg=content_fov_deg,
        twinkle_targets=twinkle_targets,
    )
