from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image
from pyproj import CRS, Transformer

DEFAULT_GEOSTATIONARY_PROJECTION_ATTRS: dict[str, float | str] = {
    "grid_mapping_name": "geostationary",
    "perspective_point_height": 35785831.0,
    "latitude_of_projection_origin": 0.0,
    "longitude_of_projection_origin": 0.0,
    "sweep_angle_axis": "x",
    "semi_major_axis": 6378137.0,
    "semi_minor_axis": 6356752.31414,
}


@dataclass(frozen=True, slots=True)
class ColorControlPoint:
    lat_deg: float
    lon_deg: float
    rgb: tuple[int, int, int]


@dataclass(frozen=True, slots=True)
class PixelControlPoint(ColorControlPoint):
    x_px: float
    y_px: float
    match_count: int = 1


@dataclass(frozen=True, slots=True)
class GeostationaryPixelMap:
    projection: CRS
    pixel_to_proj_x_coeffs: np.ndarray
    pixel_to_proj_y_coeffs: np.ndarray

    def pixel_to_proj(self, x_px: np.ndarray | float, y_px: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        x_arr = np.asarray(x_px, dtype=np.float64)
        y_arr = np.asarray(y_px, dtype=np.float64)
        a = self.pixel_to_proj_x_coeffs
        b = self.pixel_to_proj_y_coeffs
        x_proj = a[2] + a[0] * x_arr + a[1] * y_arr
        y_proj = b[2] + b[0] * x_arr + b[1] * y_arr
        return x_proj, y_proj

    def pixel_to_lonlat(self, x_px: np.ndarray | float, y_px: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        x_proj, y_proj = self.pixel_to_proj(x_px, y_px)
        to_lonlat = Transformer.from_crs(self.projection, "EPSG:4326", always_xy=True)
        lon_deg, lat_deg = to_lonlat.transform(x_proj, y_proj)
        return np.asarray(lon_deg, dtype=np.float64), np.asarray(lat_deg, dtype=np.float64)

    def pixel_grid_to_lonlat(self, width: int, height: int) -> tuple[np.ndarray, np.ndarray]:
        xs = np.arange(width, dtype=np.float64)
        ys = np.arange(height, dtype=np.float64)
        grid_x, grid_y = np.meshgrid(xs, ys)
        lon_deg, lat_deg = self.pixel_to_lonlat(grid_x, grid_y)
        return lon_deg, lat_deg


@dataclass(frozen=True, slots=True)
class GeostationaryCameraMap:
    projection: CRS
    pixel_to_proj_x_coeffs: np.ndarray
    pixel_to_proj_y_coeffs: np.ndarray
    longitude_of_projection_origin: float
    perspective_point_height: float
    sweep_angle_axis: str
    rms_pixel_error: float
    max_pixel_error: float

    def pixel_to_proj(self, x_px: np.ndarray | float, y_px: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        x_arr = np.asarray(x_px, dtype=np.float64)
        y_arr = np.asarray(y_px, dtype=np.float64)
        a = self.pixel_to_proj_x_coeffs
        b = self.pixel_to_proj_y_coeffs
        x_proj = a[2] + a[0] * x_arr + a[1] * y_arr
        y_proj = b[2] + b[0] * x_arr + b[1] * y_arr
        return x_proj, y_proj

    def pixel_to_lonlat(self, x_px: np.ndarray | float, y_px: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        x_proj, y_proj = self.pixel_to_proj(x_px, y_px)
        to_lonlat = Transformer.from_crs(self.projection, "EPSG:4326", always_xy=True)
        lon_deg, lat_deg = to_lonlat.transform(x_proj, y_proj)
        return np.asarray(lon_deg, dtype=np.float64), np.asarray(lat_deg, dtype=np.float64)

    def pixel_grid_to_lonlat(self, width: int, height: int) -> tuple[np.ndarray, np.ndarray]:
        xs = np.arange(width, dtype=np.float64)
        ys = np.arange(height, dtype=np.float64)
        grid_x, grid_y = np.meshgrid(xs, ys)
        lon_deg, lat_deg = self.pixel_to_lonlat(grid_x, grid_y)
        return lon_deg, lat_deg


@dataclass(frozen=True, slots=True)
class PixelLonLatMap:
    lon_coeffs: np.ndarray
    lat_coeffs: np.ndarray
    x_center: float
    y_center: float
    x_scale: float
    y_scale: float

    def _features(self, x_px: np.ndarray | float, y_px: np.ndarray | float) -> np.ndarray:
        x_arr = (np.asarray(x_px, dtype=np.float64) - self.x_center) / self.x_scale
        y_arr = (np.asarray(y_px, dtype=np.float64) - self.y_center) / self.y_scale
        return np.stack(
            [
                np.ones_like(x_arr),
                x_arr,
                y_arr,
                x_arr * x_arr,
                x_arr * y_arr,
                y_arr * y_arr,
            ],
            axis=0,
        )

    def pixel_to_lonlat(self, x_px: np.ndarray | float, y_px: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        feats = self._features(x_px, y_px)
        lon_deg = np.tensordot(self.lon_coeffs, feats, axes=(0, 0))
        lat_deg = np.tensordot(self.lat_coeffs, feats, axes=(0, 0))
        return np.asarray(lon_deg, dtype=np.float64), np.asarray(lat_deg, dtype=np.float64)

    def pixel_grid_to_lonlat(self, width: int, height: int) -> tuple[np.ndarray, np.ndarray]:
        xs = np.arange(width, dtype=np.float64)
        ys = np.arange(height, dtype=np.float64)
        grid_x, grid_y = np.meshgrid(xs, ys)
        return self.pixel_to_lonlat(grid_x, grid_y)


@dataclass(frozen=True, slots=True)
class PixelThinPlateSplineMap:
    lon_weights: np.ndarray
    lat_weights: np.ndarray
    x_samples: np.ndarray
    y_samples: np.ndarray
    x_center: float
    y_center: float
    x_scale: float
    y_scale: float

    def _normalized_coords(self, x_px: np.ndarray | float, y_px: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        x_arr = (np.asarray(x_px, dtype=np.float64) - self.x_center) / self.x_scale
        y_arr = (np.asarray(y_px, dtype=np.float64) - self.y_center) / self.y_scale
        return x_arr, y_arr

    def _kernel(self, x_px: np.ndarray | float, y_px: np.ndarray | float) -> np.ndarray:
        x_arr, y_arr = self._normalized_coords(x_px, y_px)
        dx = np.asarray(x_arr, dtype=np.float64)[..., None] - self.x_samples[None, ...]
        dy = np.asarray(y_arr, dtype=np.float64)[..., None] - self.y_samples[None, ...]
        r2 = dx * dx + dy * dy
        with np.errstate(divide="ignore", invalid="ignore"):
            out = np.where(r2 > 0.0, r2 * np.log(r2), 0.0)
        return out

    def _features(self, x_px: np.ndarray | float, y_px: np.ndarray | float) -> np.ndarray:
        x_arr, y_arr = self._normalized_coords(x_px, y_px)
        return np.stack([np.ones_like(x_arr), x_arr, y_arr], axis=0)

    def pixel_to_lonlat(self, x_px: np.ndarray | float, y_px: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        kernel = self._kernel(x_px, y_px)
        feats = self._features(x_px, y_px)
        lon_deg = np.tensordot(self.lon_weights["affine"], feats, axes=(0, 0)) + np.tensordot(
            self.lon_weights["radial"], kernel, axes=(0, -1)
        )
        lat_deg = np.tensordot(self.lat_weights["affine"], feats, axes=(0, 0)) + np.tensordot(
            self.lat_weights["radial"], kernel, axes=(0, -1)
        )
        return np.asarray(lon_deg, dtype=np.float64), np.asarray(lat_deg, dtype=np.float64)

    def pixel_grid_to_lonlat(self, width: int, height: int) -> tuple[np.ndarray, np.ndarray]:
        xs = np.arange(width, dtype=np.float64)
        ys = np.arange(height, dtype=np.float64)
        grid_x, grid_y = np.meshgrid(xs, ys)
        return self.pixel_to_lonlat(grid_x, grid_y)


@dataclass(frozen=True, slots=True)
class EquidistantConicFit:
    projection: CRS
    pixel_map: GeostationaryPixelMap
    longitude_of_projection_origin: float
    latitude_of_projection_origin: float
    standard_parallel_1: float
    standard_parallel_2: float
    rms_pixel_error: float
    max_pixel_error: float


_COLOR_RE = re.compile(r"#([0-9a-fA-F]{6})")


def build_projection(
    *,
    longitude_of_projection_origin: float = 0.0,
    perspective_point_height: float = 35785831.0,
    sweep_angle_axis: str = "x",
    latitude_of_projection_origin: float = 0.0,
) -> CRS:
    return CRS.from_cf(
        {
            "grid_mapping_name": "geostationary",
            "longitude_of_projection_origin": float(longitude_of_projection_origin),
            "latitude_of_projection_origin": float(latitude_of_projection_origin),
            "perspective_point_height": float(perspective_point_height),
            "sweep_angle_axis": sweep_angle_axis,
            "semi_major_axis": 6378137.0,
            "semi_minor_axis": 6356752.31414,
        }
    )


def build_equidistant_conic_projection(
    *,
    longitude_of_projection_origin: float,
    latitude_of_projection_origin: float,
    standard_parallel_1: float,
    standard_parallel_2: float,
) -> CRS:
    return CRS.from_proj4(
        f"+proj=eqdc +lat_0={float(latitude_of_projection_origin)} "
        f"+lon_0={float(longitude_of_projection_origin)} "
        f"+lat_1={float(standard_parallel_1)} +lat_2={float(standard_parallel_2)} "
        "+datum=WGS84 +units=m +no_defs"
    )


def parse_latlonmap(path: Path) -> list[ColorControlPoint]:
    control_points: list[ColorControlPoint] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        match = _COLOR_RE.search(text)
        if match is None:
            raise ValueError(f"Missing RGB color in line: {line}")
        left = text[: match.start()].rstrip(", \t")
        parts = [part.strip() for part in left.split(",")]
        if len(parts) < 2:
            raise ValueError(f"Missing latitude/longitude in line: {line}")
        try:
            lat_deg = float(parts[0])
            lon_deg = float(parts[1])
        except ValueError as exc:
            raise ValueError(f"Invalid latitude/longitude in line: {line}") from exc
        rgb = tuple(int(match.group(1)[i : i + 2], 16) for i in (0, 2, 4))
        control_points.append(ColorControlPoint(lat_deg=lat_deg, lon_deg=lon_deg, rgb=rgb))
    if not control_points:
        raise ValueError(f"No control points found in {path}")
    return control_points


def _find_exact_color_points(image_rgb: np.ndarray, rgb: tuple[int, int, int]) -> tuple[float, float, int]:
    mask = np.all(image_rgb == np.asarray(rgb, dtype=np.uint8), axis=-1)
    count = int(mask.sum())
    if count <= 0:
        raise ValueError(f"Could not find exact RGB match for {rgb}")
    ys, xs = np.nonzero(mask)
    return float(xs.mean()), float(ys.mean()), count


def locate_control_points(image: Image.Image, control_points: Sequence[ColorControlPoint]) -> list[PixelControlPoint]:
    image_rgb = np.asarray(image.convert("RGB"), dtype=np.uint8)
    located: list[PixelControlPoint] = []
    for point in control_points:
        x_px, y_px, count = _find_exact_color_points(image_rgb, point.rgb)
        located.append(
            PixelControlPoint(
                lat_deg=point.lat_deg,
                lon_deg=point.lon_deg,
                rgb=point.rgb,
                x_px=x_px,
                y_px=y_px,
                match_count=count,
            )
        )
    return located


def fit_pixel_map(
    control_points: Sequence[PixelControlPoint],
    *,
    projection: CRS | None = None,
) -> GeostationaryPixelMap:
    if len(control_points) < 3:
        raise ValueError("At least 3 control points are required")
    crs = projection or build_projection()
    to_proj = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    proj_x, proj_y = to_proj.transform(
        np.asarray([point.lon_deg for point in control_points], dtype=np.float64),
        np.asarray([point.lat_deg for point in control_points], dtype=np.float64),
    )
    proj_x = np.asarray(proj_x, dtype=np.float64)
    proj_y = np.asarray(proj_y, dtype=np.float64)
    pixels_x = np.asarray([point.x_px for point in control_points], dtype=np.float64)
    pixels_y = np.asarray([point.y_px for point in control_points], dtype=np.float64)

    design = np.column_stack([pixels_x, pixels_y, np.ones_like(pixels_x)])
    coeffs_x, *_ = np.linalg.lstsq(design, proj_x, rcond=None)
    coeffs_y, *_ = np.linalg.lstsq(design, proj_y, rcond=None)
    return GeostationaryPixelMap(
        projection=crs,
        pixel_to_proj_x_coeffs=np.asarray(coeffs_x, dtype=np.float64),
        pixel_to_proj_y_coeffs=np.asarray(coeffs_y, dtype=np.float64),
    )


def fit_pixel_map_from_image(
    image: Image.Image,
    control_map: Path,
    *,
    projection: CRS | None = None,
) -> tuple[GeostationaryPixelMap, list[PixelControlPoint]]:
    control_points = parse_latlonmap(control_map)
    located = locate_control_points(image, control_points)
    return fit_pixel_map(located, projection=projection), located


def _fit_projection_affine(
    control_points: Sequence[PixelControlPoint],
    *,
    projection: CRS,
) -> tuple[GeostationaryPixelMap, float, float]:
    to_proj = Transformer.from_crs("EPSG:4326", projection, always_xy=True)
    proj_x, proj_y = to_proj.transform(
        np.asarray([point.lon_deg for point in control_points], dtype=np.float64),
        np.asarray([point.lat_deg for point in control_points], dtype=np.float64),
    )
    proj_x = np.asarray(proj_x, dtype=np.float64)
    proj_y = np.asarray(proj_y, dtype=np.float64)
    pixels_x = np.asarray([point.x_px for point in control_points], dtype=np.float64)
    pixels_y = np.asarray([point.y_px for point in control_points], dtype=np.float64)

    design = np.column_stack([pixels_x, pixels_y, np.ones_like(pixels_x)])
    coeffs_x, *_ = np.linalg.lstsq(design, proj_x, rcond=None)
    coeffs_y, *_ = np.linalg.lstsq(design, proj_y, rcond=None)

    pred_x = design @ coeffs_x
    pred_y = design @ coeffs_y
    residual_x = pred_x - proj_x
    residual_y = pred_y - proj_y
    rms = float(np.sqrt(np.mean(residual_x * residual_x + residual_y * residual_y)))
    max_error = float(np.sqrt(np.max(residual_x * residual_x + residual_y * residual_y)))
    return (
        GeostationaryPixelMap(
            projection=projection,
            pixel_to_proj_x_coeffs=np.asarray(coeffs_x, dtype=np.float64),
            pixel_to_proj_y_coeffs=np.asarray(coeffs_y, dtype=np.float64),
        ),
        rms,
        max_error,
    )


def _fit_projection_affine_from_candidates(
    control_points: Sequence[PixelControlPoint],
    *,
    longitude_of_projection_origin: float,
    latitude_of_projection_origin: float,
    standard_parallel_1: float,
    standard_parallel_2: float,
) -> tuple[GeostationaryPixelMap, float, float]:
    projection = build_equidistant_conic_projection(
        longitude_of_projection_origin=longitude_of_projection_origin,
        latitude_of_projection_origin=latitude_of_projection_origin,
        standard_parallel_1=standard_parallel_1,
        standard_parallel_2=standard_parallel_2,
    )
    return _fit_projection_affine(control_points, projection=projection)


def _candidate_values(center: float, radius: float, step: float) -> np.ndarray:
    start = center - radius
    stop = center + radius
    values = np.arange(start, stop + step * 0.5, step, dtype=np.float64)
    values = np.unique(np.concatenate([values, np.asarray([center], dtype=np.float64)]))
    values.sort()
    return values


def _candidate_height_values(center: float, radius: float, step: float) -> np.ndarray:
    return _candidate_values(center, radius, step)


def _candidate_center_radius_step(values: np.ndarray) -> tuple[float, float, float]:
    unique = np.unique(np.sort(np.asarray(values, dtype=np.float64)))
    center = float(np.median(unique))
    if unique.size <= 1:
        return center, 0.0, 1.0
    radius = float(np.max(np.abs(unique - center)))
    diffs = np.diff(unique)
    step = float(np.min(diffs)) if diffs.size > 0 else 1.0
    return center, radius, step


def _default_equidistant_conic_centers(control_points: Sequence[PixelControlPoint]) -> tuple[float, float, float, float]:
    latitudes = np.asarray([point.lat_deg for point in control_points], dtype=np.float64)
    longitudes = _unwrap_longitudes(np.asarray([point.lon_deg for point in control_points], dtype=np.float64))
    lon0_center = float(np.median(longitudes))
    lat0_center = float(np.median(latitudes))
    lat1_center = float(np.percentile(latitudes, 30.0))
    lat2_center = float(np.percentile(latitudes, 70.0))
    if lat1_center > lat2_center:
        lat1_center, lat2_center = lat2_center, lat1_center
    return lon0_center, lat0_center, lat1_center, lat2_center


def _refine_equidistant_conic_projection_fit(
    control_points: Sequence[PixelControlPoint],
    *,
    longitude_of_projection_origin_center: float,
    latitude_of_projection_origin_center: float,
    standard_parallel_1_center: float,
    standard_parallel_2_center: float,
    longitude_of_projection_origin_radius: float,
    longitude_of_projection_origin_step: float,
    latitude_of_projection_origin_radius: float,
    latitude_of_projection_origin_step: float,
    standard_parallel_1_radius: float,
    standard_parallel_1_step: float,
    standard_parallel_2_radius: float,
    standard_parallel_2_step: float,
) -> EquidistantConicFit:
    best_fit: EquidistantConicFit | None = None
    best_rms = float("inf")
    best_max = float("inf")

    for lon0 in _candidate_values(
        longitude_of_projection_origin_center,
        longitude_of_projection_origin_radius,
        longitude_of_projection_origin_step,
    ):
        for lat0 in _candidate_values(
            latitude_of_projection_origin_center,
            latitude_of_projection_origin_radius,
            latitude_of_projection_origin_step,
        ):
            for lat1 in _candidate_values(
                standard_parallel_1_center,
                standard_parallel_1_radius,
                standard_parallel_1_step,
            ):
                for lat2 in _candidate_values(
                    standard_parallel_2_center,
                    standard_parallel_2_radius,
                    standard_parallel_2_step,
                ):
                    if abs(float(lat2) - float(lat1)) < 1e-9:
                        continue
                    if float(lat1) > float(lat2):
                        continue
                    pixel_map, rms, max_error = _fit_projection_affine_from_candidates(
                        control_points,
                        longitude_of_projection_origin=float(lon0),
                        latitude_of_projection_origin=float(lat0),
                        standard_parallel_1=float(lat1),
                        standard_parallel_2=float(lat2),
                    )
                    if rms < best_rms:
                        best_fit = EquidistantConicFit(
                            projection=pixel_map.projection,
                            pixel_map=pixel_map,
                            longitude_of_projection_origin=float(lon0),
                            latitude_of_projection_origin=float(lat0),
                            standard_parallel_1=float(lat1),
                            standard_parallel_2=float(lat2),
                            rms_pixel_error=float(rms),
                            max_pixel_error=float(max_error),
                        )
                        best_rms = float(rms)
                        best_max = float(max_error)

    if best_fit is None:
        raise ValueError("No equidistant conic projection candidate found")
    _ = best_max
    return best_fit


def _refine_geostationary_projection_fit(
    control_points: Sequence[PixelControlPoint],
    *,
    lon0_center: float,
    height_center: float,
    sweep_angle_axis: str,
    lon0_radius: float,
    lon0_step: float,
    height_radius: float,
    height_step: float,
) -> tuple[GeostationaryPixelMap, float, float]:
    best_map: GeostationaryPixelMap | None = None
    best_rms = float("inf")
    best_max = float("inf")
    for lon0 in _candidate_values(lon0_center, lon0_radius, lon0_step):
        for height in _candidate_height_values(height_center, height_radius, height_step):
            projection = build_projection(
                longitude_of_projection_origin=float(lon0),
                perspective_point_height=float(height),
                sweep_angle_axis=sweep_angle_axis,
            )
            pixel_map, rms, max_error = _fit_projection_affine(control_points, projection=projection)
            if rms < best_rms:
                best_map = pixel_map
                best_rms = rms
                best_max = max_error
    if best_map is None:
        raise ValueError("No geostationary projection candidate found")
    return best_map, best_rms, best_max


def _search_equidistant_conic_projection_fit(
    control_points: Sequence[PixelControlPoint],
) -> EquidistantConicFit:
    lon0_center, lat0_center, lat1_center, lat2_center = _default_equidistant_conic_centers(control_points)
    stage1 = _refine_equidistant_conic_projection_fit(
        control_points,
        longitude_of_projection_origin_center=lon0_center,
        latitude_of_projection_origin_center=lat0_center,
        standard_parallel_1_center=lat1_center,
        standard_parallel_2_center=lat2_center,
        longitude_of_projection_origin_radius=25.0,
        longitude_of_projection_origin_step=5.0,
        latitude_of_projection_origin_radius=15.0,
        latitude_of_projection_origin_step=5.0,
        standard_parallel_1_radius=15.0,
        standard_parallel_1_step=5.0,
        standard_parallel_2_radius=15.0,
        standard_parallel_2_step=5.0,
    )
    stage2 = _refine_equidistant_conic_projection_fit(
        control_points,
        longitude_of_projection_origin_center=stage1.longitude_of_projection_origin,
        latitude_of_projection_origin_center=stage1.latitude_of_projection_origin,
        standard_parallel_1_center=stage1.standard_parallel_1,
        standard_parallel_2_center=stage1.standard_parallel_2,
        longitude_of_projection_origin_radius=6.0,
        longitude_of_projection_origin_step=2.0,
        latitude_of_projection_origin_radius=6.0,
        latitude_of_projection_origin_step=2.0,
        standard_parallel_1_radius=6.0,
        standard_parallel_1_step=2.0,
        standard_parallel_2_radius=6.0,
        standard_parallel_2_step=2.0,
    )
    stage3 = _refine_equidistant_conic_projection_fit(
        control_points,
        longitude_of_projection_origin_center=stage2.longitude_of_projection_origin,
        latitude_of_projection_origin_center=stage2.latitude_of_projection_origin,
        standard_parallel_1_center=stage2.standard_parallel_1,
        standard_parallel_2_center=stage2.standard_parallel_2,
        longitude_of_projection_origin_radius=1.0,
        longitude_of_projection_origin_step=0.5,
        latitude_of_projection_origin_radius=1.0,
        latitude_of_projection_origin_step=0.5,
        standard_parallel_1_radius=1.0,
        standard_parallel_1_step=0.5,
        standard_parallel_2_radius=1.0,
        standard_parallel_2_step=0.5,
    )
    return stage3


def fit_geostationary_camera_map(
    control_points: Sequence[PixelControlPoint],
    *,
    lon0_hint: float = 0.0,
    height_hint: float = 35785831.0,
    sweep_hint: str = "x",
) -> GeostationaryCameraMap:
    if len(control_points) < 4:
        raise ValueError("At least 4 control points are required for geostationary fitting")
    search_specs = [
        (float(lon0_hint), float(height_hint), str(sweep_hint)),
        (float(lon0_hint), float(height_hint), "y" if str(sweep_hint) == "x" else "x"),
    ]
    best_map: GeostationaryPixelMap | None = None
    best_rms = float("inf")
    best_max = float("inf")
    best_lon0 = 0.0
    best_height = 35785831.0
    best_sweep = "x"

    for lon0_center, height_center, sweep_axis in search_specs:
        stage1_map, stage1_rms, stage1_max = _refine_geostationary_projection_fit(
            control_points,
            lon0_center=lon0_center,
            height_center=height_center,
            sweep_angle_axis=sweep_axis,
            lon0_radius=30.0,
            lon0_step=2.5,
            height_radius=750000.0,
            height_step=50000.0,
        )
        lon0_guess = float(stage1_map.projection.to_cf()["longitude_of_projection_origin"])
        height_guess = float(stage1_map.projection.to_cf()["perspective_point_height"])
        stage2_map, stage2_rms, stage2_max = _refine_geostationary_projection_fit(
            control_points,
            lon0_center=lon0_guess,
            height_center=height_guess,
            sweep_angle_axis=sweep_axis,
            lon0_radius=5.0,
            lon0_step=0.5,
            height_radius=150000.0,
            height_step=10000.0,
        )
        lon0_guess = float(stage2_map.projection.to_cf()["longitude_of_projection_origin"])
        height_guess = float(stage2_map.projection.to_cf()["perspective_point_height"])
        stage3_map, stage3_rms, stage3_max = _refine_geostationary_projection_fit(
            control_points,
            lon0_center=lon0_guess,
            height_center=height_guess,
            sweep_angle_axis=sweep_axis,
            lon0_radius=1.0,
            lon0_step=0.1,
            height_radius=30000.0,
            height_step=2000.0,
        )
        candidate_map = stage3_map
        candidate_rms = stage3_rms
        candidate_max = stage3_max
        candidate_lon0 = float(candidate_map.projection.to_cf()["longitude_of_projection_origin"])
        candidate_height = float(candidate_map.projection.to_cf()["perspective_point_height"])
        if candidate_rms < best_rms:
            best_map = candidate_map
            best_rms = candidate_rms
            best_max = candidate_max
            best_lon0 = candidate_lon0
            best_height = candidate_height
            best_sweep = sweep_axis

        # Keep stage1/stage2 variables alive for readability; they may be useful in debugging.
        _ = (stage1_rms, stage1_max, stage2_rms, stage2_max)

    if best_map is None:
        raise ValueError("Failed to fit geostationary projection")
    return GeostationaryCameraMap(
        projection=best_map.projection,
        pixel_to_proj_x_coeffs=best_map.pixel_to_proj_x_coeffs,
        pixel_to_proj_y_coeffs=best_map.pixel_to_proj_y_coeffs,
        longitude_of_projection_origin=best_lon0,
        perspective_point_height=best_height,
        sweep_angle_axis=best_sweep,
        rms_pixel_error=best_rms,
        max_pixel_error=best_max,
    )


def fit_geostationary_camera_map_from_image(
    image: Image.Image,
    control_map: Path,
    *,
    lon0_hint: float = 0.0,
    height_hint: float = 35785831.0,
    sweep_hint: str = "x",
) -> tuple[GeostationaryCameraMap, list[PixelControlPoint]]:
    control_points = parse_latlonmap(control_map)
    located = locate_control_points(image, control_points)
    return (
        fit_geostationary_camera_map(
            located,
            lon0_hint=lon0_hint,
            height_hint=height_hint,
            sweep_hint=sweep_hint,
        ),
        located,
    )


def fit_equidistant_conic_projection(
    control_points: Sequence[PixelControlPoint],
    *,
    longitude_of_projection_origin_candidates: Sequence[float] | None = None,
    latitude_of_projection_origin_candidates: Sequence[float] | None = None,
    standard_parallel_1_candidates: Sequence[float] | None = None,
    standard_parallel_2_candidates: Sequence[float] | None = None,
) -> EquidistantConicFit:
    if (
        longitude_of_projection_origin_candidates is None
        and latitude_of_projection_origin_candidates is None
        and standard_parallel_1_candidates is None
        and standard_parallel_2_candidates is None
    ):
        return _search_equidistant_conic_projection_fit(control_points)

    lon0_values = np.asarray(
        list(longitude_of_projection_origin_candidates or []),
        dtype=np.float64,
    )
    lat0_values = np.asarray(
        list(latitude_of_projection_origin_candidates or []),
        dtype=np.float64,
    )
    lat1_values = np.asarray(
        list(standard_parallel_1_candidates or []),
        dtype=np.float64,
    )
    lat2_values = np.asarray(
        list(standard_parallel_2_candidates or []),
        dtype=np.float64,
    )
    if lon0_values.size == 0 or lat0_values.size == 0 or lat1_values.size == 0 or lat2_values.size == 0:
        raise ValueError("All Equidistant Conic candidate lists must be non-empty when provided")
    lon0_center, lon0_radius, lon0_step = _candidate_center_radius_step(lon0_values)
    lat0_center, lat0_radius, lat0_step = _candidate_center_radius_step(lat0_values)
    lat1_center, lat1_radius, lat1_step = _candidate_center_radius_step(lat1_values)
    lat2_center, lat2_radius, lat2_step = _candidate_center_radius_step(lat2_values)
    return _refine_equidistant_conic_projection_fit(
        control_points,
        longitude_of_projection_origin_center=lon0_center,
        latitude_of_projection_origin_center=lat0_center,
        standard_parallel_1_center=lat1_center,
        standard_parallel_2_center=lat2_center,
        longitude_of_projection_origin_radius=lon0_radius,
        longitude_of_projection_origin_step=lon0_step,
        latitude_of_projection_origin_radius=lat0_radius,
        latitude_of_projection_origin_step=lat0_step,
        standard_parallel_1_radius=lat1_radius,
        standard_parallel_1_step=lat1_step,
        standard_parallel_2_radius=lat2_radius,
        standard_parallel_2_step=lat2_step,
    )


def fit_equidistant_conic_projection_from_image(
    image: Image.Image,
    control_map: Path,
    *,
    longitude_of_projection_origin_candidates: Sequence[float] | None = None,
    latitude_of_projection_origin_candidates: Sequence[float] | None = None,
    standard_parallel_1_candidates: Sequence[float] | None = None,
    standard_parallel_2_candidates: Sequence[float] | None = None,
) -> tuple[EquidistantConicFit, list[PixelControlPoint]]:
    control_points = parse_latlonmap(control_map)
    located = locate_control_points(image, control_points)
    return (
        fit_equidistant_conic_projection(
            located,
            longitude_of_projection_origin_candidates=longitude_of_projection_origin_candidates,
            latitude_of_projection_origin_candidates=latitude_of_projection_origin_candidates,
            standard_parallel_1_candidates=standard_parallel_1_candidates,
            standard_parallel_2_candidates=standard_parallel_2_candidates,
        ),
        located,
    )


def _unwrap_longitudes(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return values
    reference = float(np.median(values))
    unwrapped = values.copy()
    for index, value in enumerate(unwrapped):
        while value - reference > 180.0:
            value -= 360.0
        while value - reference < -180.0:
            value += 360.0
        unwrapped[index] = value
    return unwrapped


def fit_pixel_lonlat_polynomial(
    control_points: Sequence[PixelControlPoint],
) -> PixelLonLatMap:
    if len(control_points) < 6:
        raise ValueError("At least 6 control points are required for a quadratic fit")
    x = np.asarray([point.x_px for point in control_points], dtype=np.float64)
    y = np.asarray([point.y_px for point in control_points], dtype=np.float64)
    lon = _unwrap_longitudes(np.asarray([point.lon_deg for point in control_points], dtype=np.float64))
    lat = np.asarray([point.lat_deg for point in control_points], dtype=np.float64)
    x_center = float(np.mean(x))
    y_center = float(np.mean(y))
    x_scale = float(np.std(x) or 1.0)
    y_scale = float(np.std(y) or 1.0)
    xn = (x - x_center) / x_scale
    yn = (y - y_center) / y_scale
    design = np.column_stack([np.ones_like(xn), xn, yn, xn * xn, xn * yn, yn * yn])
    lon_coeffs, *_ = np.linalg.lstsq(design, lon, rcond=None)
    lat_coeffs, *_ = np.linalg.lstsq(design, lat, rcond=None)
    return PixelLonLatMap(
        lon_coeffs=np.asarray(lon_coeffs, dtype=np.float64),
        lat_coeffs=np.asarray(lat_coeffs, dtype=np.float64),
        x_center=x_center,
        y_center=y_center,
        x_scale=x_scale,
        y_scale=y_scale,
    )


def fit_pixel_lonlat_polynomial_from_image(
    image: Image.Image,
    control_map: Path,
) -> tuple[PixelLonLatMap, list[PixelControlPoint]]:
    control_points = parse_latlonmap(control_map)
    located = locate_control_points(image, control_points)
    return fit_pixel_lonlat_polynomial(located), located


def _tps_kernel_matrix(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    r2 = dx * dx + dy * dy
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.where(r2 > 0.0, r2 * np.log(r2), 0.0)
    return k


def _fit_tps_weights(x: np.ndarray, y: np.ndarray, values: np.ndarray, *, regularization: float) -> dict[str, np.ndarray]:
    n = x.size
    if n < 3:
        raise ValueError("At least 3 control points are required for TPS fitting")
    k = _tps_kernel_matrix(x, y)
    if regularization > 0.0:
        k = k + np.eye(n, dtype=np.float64) * float(regularization)
    p = np.column_stack([np.ones(n, dtype=np.float64), x, y])
    system = np.block(
        [
            [k, p],
            [p.T, np.zeros((3, 3), dtype=np.float64)],
        ]
    )
    rhs = np.concatenate([values, np.zeros(3, dtype=np.float64)])
    weights = np.linalg.solve(system, rhs)
    return {
        "radial": np.asarray(weights[:n], dtype=np.float64),
        "affine": np.asarray(weights[n:], dtype=np.float64),
    }


def fit_pixel_lonlat_tps(
    control_points: Sequence[PixelControlPoint],
    *,
    regularization: float = 1e-6,
) -> PixelThinPlateSplineMap:
    if len(control_points) < 3:
        raise ValueError("At least 3 control points are required for TPS fitting")
    x = np.asarray([point.x_px for point in control_points], dtype=np.float64)
    y = np.asarray([point.y_px for point in control_points], dtype=np.float64)
    lon = _unwrap_longitudes(np.asarray([point.lon_deg for point in control_points], dtype=np.float64))
    lat = np.asarray([point.lat_deg for point in control_points], dtype=np.float64)
    x_center = float(np.mean(x))
    y_center = float(np.mean(y))
    x_scale = float(np.std(x) or 1.0)
    y_scale = float(np.std(y) or 1.0)
    xn = (x - x_center) / x_scale
    yn = (y - y_center) / y_scale
    lon_weights = _fit_tps_weights(xn, yn, lon, regularization=regularization)
    lat_weights = _fit_tps_weights(xn, yn, lat, regularization=regularization)
    return PixelThinPlateSplineMap(
        lon_weights=lon_weights,
        lat_weights=lat_weights,
        x_samples=xn,
        y_samples=yn,
        x_center=x_center,
        y_center=y_center,
        x_scale=x_scale,
        y_scale=y_scale,
    )


def fit_pixel_lonlat_tps_from_image(
    image: Image.Image,
    control_map: Path,
    *,
    regularization: float = 1e-6,
) -> tuple[PixelThinPlateSplineMap, list[PixelControlPoint]]:
    control_points = parse_latlonmap(control_map)
    located = locate_control_points(image, control_points)
    return fit_pixel_lonlat_tps(located, regularization=regularization), located


def estimate_pixel_to_lonlat(
    image: Image.Image,
    control_map: Path,
    *,
    projection: CRS | None = None,
) -> tuple[GeostationaryPixelMap, list[PixelControlPoint], np.ndarray, np.ndarray]:
    pixel_map, located = fit_pixel_map_from_image(image, control_map, projection=projection)
    lon_deg, lat_deg = pixel_map.pixel_grid_to_lonlat(image.width, image.height)
    return pixel_map, located, lon_deg, lat_deg
