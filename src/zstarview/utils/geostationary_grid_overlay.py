from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from .geostationary_image_mapping import build_projection, estimate_pixel_to_lonlat


@dataclass(frozen=True, slots=True)
class GeostationaryGridData:
    lon_deg: np.ndarray
    lat_deg: np.ndarray


def load_lonlat_grid_from_npz(path: Path) -> GeostationaryGridData:
    with np.load(path, allow_pickle=False) as data:
        if "lon_deg" not in data or "lat_deg" not in data:
            raise ValueError(f"{path} does not contain lon_deg/lat_deg arrays")
        lon_deg = np.asarray(data["lon_deg"], dtype=np.float64)
        lat_deg = np.asarray(data["lat_deg"], dtype=np.float64)
    if lon_deg.shape != lat_deg.shape:
        raise ValueError("lon_deg and lat_deg arrays must have the same shape")
    return GeostationaryGridData(lon_deg=lon_deg, lat_deg=lat_deg)


def load_lonlat_grid_for_image(
    image_path: Path,
    *,
    grid_npz: Path | None = None,
    control_map: Path | None = None,
    longitude_of_projection_origin: float = 0.0,
    perspective_point_height: float = 35785831.0,
    sweep_angle_axis: str = "x",
) -> GeostationaryGridData:
    if grid_npz is not None:
        return load_lonlat_grid_from_npz(grid_npz)
    if control_map is None:
        raise ValueError("A grid .npz file or a control map is required")
    with Image.open(image_path) as image:
        _, _, lon_deg, lat_deg = estimate_pixel_to_lonlat(
            image,
            control_map,
            projection=build_projection(
                longitude_of_projection_origin=longitude_of_projection_origin,
                perspective_point_height=perspective_point_height,
                sweep_angle_axis=sweep_angle_axis,
            ),
        )
    return GeostationaryGridData(lon_deg=np.asarray(lon_deg, dtype=np.float64), lat_deg=np.asarray(lat_deg, dtype=np.float64))


def _wrap_lon_delta(lon_deg: np.ndarray, target_deg: float) -> np.ndarray:
    return ((lon_deg - target_deg + 180.0) % 360.0) - 180.0


def _target_levels(step_deg: int, major_step_deg: int, *, min_deg: int, max_deg: int) -> list[tuple[float, bool]]:
    levels: list[tuple[float, bool]] = []
    for level in range(min_deg, max_deg + 1, int(step_deg)):
        levels.append((float(level), (level % int(major_step_deg)) == 0))
    return levels


def _row_crossings(values: np.ndarray, target: float, *, wrap: bool) -> list[float]:
    if wrap:
        diff = _wrap_lon_delta(values, target)
    else:
        diff = values - target
    finite = np.isfinite(diff)
    if int(finite.sum()) < 2:
        return []
    d0 = diff[:-1]
    d1 = diff[1:]
    valid = np.isfinite(d0) & np.isfinite(d1)
    crossing = valid & ((d0 == 0.0) | (d1 == 0.0) | (d0 * d1 < 0.0))
    idxs = np.flatnonzero(crossing)
    if idxs.size == 0:
        return []
    points: list[float] = []
    for idx in idxs:
        left = float(diff[idx])
        right = float(diff[idx + 1])
        denom = right - left
        if denom == 0.0:
            continue
        t = -left / denom
        if 0.0 <= t <= 1.0:
            points.append(float(idx + t))
    return points


def _col_crossings(values: np.ndarray, target: float) -> list[float]:
    diff = values - target
    finite = np.isfinite(diff)
    if int(finite.sum()) < 2:
        return []
    d0 = diff[:-1]
    d1 = diff[1:]
    valid = np.isfinite(d0) & np.isfinite(d1)
    crossing = valid & ((d0 == 0.0) | (d1 == 0.0) | (d0 * d1 < 0.0))
    idxs = np.flatnonzero(crossing)
    if idxs.size == 0:
        return []
    points: list[float] = []
    for idx in idxs:
        left = float(diff[idx])
        right = float(diff[idx + 1])
        denom = right - left
        if denom == 0.0:
            continue
        t = -left / denom
        if 0.0 <= t <= 1.0:
            points.append(float(idx + t))
    return points


def _extract_lon_line_points(lon_deg: np.ndarray, target: float) -> list[tuple[float, float]]:
    points: list[tuple[float, float]] = []
    active: list[tuple[float, float]] = []
    previous_x: float | None = None
    for y_idx in range(lon_deg.shape[0]):
        candidates = _row_crossings(lon_deg[y_idx], target, wrap=True)
        if not candidates:
            if len(active) >= 2:
                points.extend(active)
                points.append((float("nan"), float("nan")))
            active = []
            previous_x = None
            continue
        x_val = min(candidates, key=lambda value: abs(value - previous_x)) if previous_x is not None else min(
            candidates,
            key=lambda value: abs(value - (lon_deg.shape[1] - 1) * 0.5),
        )
        active.append((x_val, float(y_idx)))
        previous_x = x_val
    if len(active) >= 2:
        points.extend(active)
    return points


def _extract_lat_line_points(lat_deg: np.ndarray, target: float) -> list[tuple[float, float]]:
    points: list[tuple[float, float]] = []
    active: list[tuple[float, float]] = []
    previous_y: float | None = None
    for x_idx in range(lat_deg.shape[1]):
        candidates = _col_crossings(lat_deg[:, x_idx], target)
        if not candidates:
            if len(active) >= 2:
                points.extend(active)
                points.append((float("nan"), float("nan")))
            active = []
            previous_y = None
            continue
        y_val = min(candidates, key=lambda value: abs(value - previous_y)) if previous_y is not None else min(
            candidates,
            key=lambda value: abs(value - (lat_deg.shape[0] - 1) * 0.5),
        )
        active.append((float(x_idx), y_val))
        previous_y = y_val
    if len(active) >= 2:
        points.extend(active)
    return points


def _draw_polyline(draw: ImageDraw.ImageDraw, points: list[tuple[float, float]], *, color: tuple[int, int, int], width: int) -> None:
    run: list[tuple[float, float]] = []
    for x, y in points:
        if not (np.isfinite(x) and np.isfinite(y)):
            if len(run) >= 2:
                draw.line(run, fill=color, width=width)
            run = []
            continue
        run.append((float(x), float(y)))
    if len(run) >= 2:
        draw.line(run, fill=color, width=width)


def draw_geostationary_latlon_grid(
    image: Image.Image,
    lon_deg: np.ndarray,
    lat_deg: np.ndarray,
    *,
    step_deg: int = 10,
    major_step_deg: int = 30,
    minor_color: tuple[int, int, int] = (0, 0, 0),
    major_color: tuple[int, int, int] = (255, 0, 0),
    minor_width: int = 1,
    major_width: int = 2,
) -> Image.Image:
    if lon_deg.shape != lat_deg.shape:
        raise ValueError("lon_deg and lat_deg arrays must have the same shape")
    overlay = image.convert("RGBA")
    draw = ImageDraw.Draw(overlay)

    lon_levels = _target_levels(step_deg, major_step_deg, min_deg=-180, max_deg=180)
    lat_levels = _target_levels(step_deg, major_step_deg, min_deg=-90, max_deg=90)

    for level, is_major in lon_levels:
        points = _extract_lon_line_points(lon_deg, level)
        _draw_polyline(
            draw,
            points,
            color=major_color if is_major else minor_color,
            width=major_width if is_major else minor_width,
        )

    for level, is_major in lat_levels:
        points = _extract_lat_line_points(lat_deg, level)
        _draw_polyline(
            draw,
            points,
            color=major_color if is_major else minor_color,
            width=major_width if is_major else minor_width,
        )

    return overlay.convert(image.mode)
