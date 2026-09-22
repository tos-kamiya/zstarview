"""Compare the current NumPy star rectangle rasterizer with a Numba kernel.

This benchmark intentionally measures only the regular-star RGB canvas stage.
It excludes coordinate projection, QImage construction, and QPainter drawing.
The Numba implementation is a prototype and is not part of the production
render path.
"""

from __future__ import annotations

import argparse
import time

import numpy as np
from numba import njit
from PySide6.QtCore import Qt
from PySide6.QtGui import QImage, QPainter

from zstarview.render.stars import _canvas_to_argb32_premultiplied_image


@njit(cache=True)
def _rasterize_rectangles_numba(
    x0: np.ndarray,
    y0: np.ndarray,
    size_px: np.ndarray,
    colors: np.ndarray,
    width: int,
    height: int,
) -> np.ndarray:
    canvas = np.zeros((height, width, 3), dtype=np.float32)
    for index in range(size_px.shape[0]):
        if size_px[index] == 2 and (
            x0[index] < 0
            or y0[index] < 0
            or x0[index] + size_px[index] > width
            or y0[index] + size_px[index] > height
        ):
            continue
        left = max(0, x0[index])
        top = max(0, y0[index])
        right = min(width, x0[index] + size_px[index])
        bottom = min(height, y0[index] + size_px[index])
        for y in range(top, bottom):
            for x in range(left, right):
                canvas[y, x, 0] += colors[index, 0]
                canvas[y, x, 1] += colors[index, 1]
                canvas[y, x, 2] += colors[index, 2]
    return canvas


def _rasterize_rectangles_numpy(
    x0: np.ndarray,
    y0: np.ndarray,
    size_px: np.ndarray,
    colors: np.ndarray,
    width: int,
    height: int,
) -> np.ndarray:
    """Mirror the regular-star portion of render/stars.py."""
    canvas = np.zeros((height, width, 3), dtype=np.float32)
    x1 = x0 + size_px
    y1 = y0 + size_px
    x0_clamped = np.clip(x0, 0, width)
    y0_clamped = np.clip(y0, 0, height)
    x1_clamped = np.clip(x1, 0, width)
    y1_clamped = np.clip(y1, 0, height)
    valid_base = (
        (x1_clamped > x0_clamped)
        & (y1_clamped > y0_clamped)
        & (size_px > 0)
    )

    single_indices = np.nonzero(valid_base & (size_px == 1))[0]
    if single_indices.size > 0:
        single_layer = np.zeros_like(canvas)
        flat_single = single_layer.reshape(-1, 3)
        flat_idx = y0_clamped[single_indices] * width + x0_clamped[single_indices]
        np.add.at(flat_single, flat_idx, colors[single_indices])
        canvas += single_layer

    size2_indices = np.nonzero(
        valid_base
        & (size_px == 2)
        & (x0 >= 0)
        & (y0 >= 0)
        & (x1 <= width)
        & (y1 <= height)
    )[0]
    if size2_indices.size > 0:
        size2_layer = np.zeros_like(canvas)
        flat_size2 = size2_layer.reshape(-1, 3)
        base_idx = y0[size2_indices] * width + x0[size2_indices]
        colors_size2 = colors[size2_indices]
        np.add.at(flat_size2, base_idx, colors_size2)
        np.add.at(flat_size2, base_idx + 1, colors_size2)
        np.add.at(flat_size2, base_idx + width, colors_size2)
        np.add.at(flat_size2, base_idx + width + 1, colors_size2)
        canvas += size2_layer

    rectangle_indices = np.nonzero(valid_base & (size_px >= 3) & (size_px <= 6))[0]
    for index in rectangle_indices:
        canvas[
            y0_clamped[index] : y1_clamped[index],
            x0_clamped[index] : x1_clamped[index],
            :,
        ] += colors[index]
    return canvas


def _make_input(star_count: int, width: int, height: int, seed: int) -> tuple[np.ndarray, ...]:
    rng = np.random.default_rng(seed)
    x = rng.integers(-20, width + 20, size=star_count, dtype=np.int32)
    y = rng.integers(-20, height + 20, size=star_count, dtype=np.int32)
    size_px = rng.integers(1, 7, size=star_count, dtype=np.int32)
    colors = rng.random((star_count, 3), dtype=np.float32)
    return x - (size_px // 2), y - (size_px // 2), size_px, colors


def _median_seconds(function, args: tuple[object, ...], repeats: int) -> float:
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        function(*args)
        samples.append(time.perf_counter() - start)
    return float(np.median(np.asarray(samples, dtype=np.float64)))


def _canvas_to_qimage_old(canvas: np.ndarray, width: int, height: int) -> QImage:
    canvas_uint8 = np.ascontiguousarray(np.clip(canvas, 0.0, 255.0).astype(np.uint8))
    rgba = np.zeros((height, width, 4), dtype=np.uint8)
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
    return QImage(rgba.data, width, height, width * 4, QImage.Format_RGBA8888).copy()


def _canvas_to_qimage_new(canvas: np.ndarray, width: int, height: int) -> QImage:
    del width, height
    return _canvas_to_argb32_premultiplied_image(canvas)


def _compose_qimage(image: QImage, width: int, height: int) -> None:
    target = QImage(width, height, QImage.Format_ARGB32_Premultiplied)
    target.fill(Qt.GlobalColor.transparent)
    painter = QPainter(target)
    painter.drawImage(0, 0, image)
    painter.end()


def _median_full_seconds(
    rasterizer,
    image_converter,
    args: tuple[object, ...],
    repeats: int,
) -> float:
    samples = []
    width = int(args[4])
    height = int(args[5])
    for _ in range(repeats):
        start = time.perf_counter()
        canvas = rasterizer(*args)
        image = image_converter(canvas, width, height)
        _compose_qimage(image, width, height)
        samples.append(time.perf_counter() - start)
    return float(np.median(np.asarray(samples, dtype=np.float64)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stars", type=int, default=20000)
    parser.add_argument("--width", type=int, default=1920)
    parser.add_argument("--height", type=int, default=1080)
    parser.add_argument("--repeats", type=int, default=7)
    args = parser.parse_args()

    x0, y0, size_px, colors = _make_input(args.stars, args.width, args.height, seed=20260922)
    kernel_args = (x0, y0, size_px, colors, args.width, args.height)

    # Exclude one-time LLVM compilation from the measured Numba timings.
    numba_result = _rasterize_rectangles_numba(*kernel_args)
    numpy_result = _rasterize_rectangles_numpy(*kernel_args)
    np.testing.assert_allclose(numba_result, numpy_result, rtol=0.0, atol=1.0e-6)

    numpy_seconds = _median_seconds(_rasterize_rectangles_numpy, kernel_args, args.repeats)
    numba_seconds = _median_seconds(_rasterize_rectangles_numba, kernel_args, args.repeats)
    numpy_full_old_seconds = _median_full_seconds(
        _rasterize_rectangles_numpy, _canvas_to_qimage_old, kernel_args, args.repeats
    )
    numba_full_old_seconds = _median_full_seconds(
        _rasterize_rectangles_numba, _canvas_to_qimage_old, kernel_args, args.repeats
    )
    numpy_full_new_seconds = _median_full_seconds(
        _rasterize_rectangles_numpy,
        _canvas_to_qimage_new,
        kernel_args,
        args.repeats,
    )
    numba_full_new_seconds = _median_full_seconds(
        _rasterize_rectangles_numba,
        _canvas_to_qimage_new,
        kernel_args,
        args.repeats,
    )
    print(f"stars={args.stars} viewport={args.width}x{args.height} repeats={args.repeats}")
    print(f"numpy_median_s={numpy_seconds:.6f}")
    print(f"numba_median_s={numba_seconds:.6f}")
    print(f"speedup={numpy_seconds / numba_seconds:.2f}x")
    print(f"numpy_full_old_median_s={numpy_full_old_seconds:.6f}")
    print(f"numba_full_old_median_s={numba_full_old_seconds:.6f}")
    print(f"numpy_full_new_median_s={numpy_full_new_seconds:.6f}")
    print(f"numba_full_new_median_s={numba_full_new_seconds:.6f}")
    print(f"old_full_speedup={numpy_full_old_seconds / numba_full_old_seconds:.2f}x")
    print(f"new_full_speedup={numpy_full_new_seconds / numba_full_new_seconds:.2f}x")
    print(f"numpy_format_speedup={numpy_full_old_seconds / numpy_full_new_seconds:.2f}x")
    print(f"numba_format_speedup={numba_full_old_seconds / numba_full_new_seconds:.2f}x")
    print("scope=regular rectangle rasterization; full includes RGBA, QImage, QPainter")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
