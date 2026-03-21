from __future__ import annotations

import argparse
import datetime as dt
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from PIL import Image
from pyproj import CRS

from zstarview.clouddisc.config import CloudDiscConfig
from zstarview.clouddisc.geo_area import GeoArea
from zstarview.clouddisc.projectors.az import az_project_lonlat_grid
from zstarview.clouddisc.providers._goes_abi import DATA_VAR as GOES_DATA_VAR
from zstarview.clouddisc.providers._goes_abi import GRID_VAR as GOES_GRID_VAR
from zstarview.clouddisc.providers._goes_abi import load_cmi_with_area
from zstarview.clouddisc.providers._hima_isatss import DATA_VAR as HIMA_DATA_VAR
from zstarview.clouddisc.providers._hima_isatss import GRID_VAR as HIMA_GRID_VAR
from zstarview.clouddisc.providers._hima_isatss import stitch_tiles_from_paths
from zstarview.clouddisc.providers.hima import HimaProvider
from zstarview.clouddisc.render.grayscale import convert_bt_to_la_image
from zstarview.clouddisc.sampling.bt_sampler import build_bt_sampler
from zstarview.clouddisc.sampling.estimate_bt_warm_cold import estimate_bt_cold_hybrid, estimate_bt_warm_from_equator_band

try:
    from pyresample.geometry import AreaDefinition as PyresampleAreaDefinition
except Exception:  # pragma: no cover - optional dependency for manual comparison only
    PyresampleAreaDefinition = None


@dataclass
class RenderStats:
    label: str
    bt_warm: float
    bt_cold: float
    coverage_ratio: float
    inside_pixels: int
    inside_valid_pixels: int
    eq_samples: int
    bt_min: float
    bt_max: float
    bt_mean: float


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare current cloud-disc output with a pyresample-backed area definition.",
    )
    parser.add_argument("--source", choices=("goes", "hima"), required=True, help="Input source type.")
    parser.add_argument("--input", type=Path, required=True, help="GOES NetCDF file or Himawari tile directory.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Directory to write PNGs and summary.json.")
    parser.add_argument("--lat", type=float, required=True, help="Observer latitude in degrees.")
    parser.add_argument("--lon", type=float, required=True, help="Observer longitude in degrees.")
    parser.add_argument("--alt", type=float, required=True, help="View altitude in degrees.")
    parser.add_argument("--az", type=float, required=True, help="View azimuth in degrees.")
    parser.add_argument("--radius-px", type=int, default=512, help="Output cloud-disc radius in pixels.")
    parser.add_argument("--edge-fov-deg", type=float, default=90.0, help="Render edge field of view.")
    parser.add_argument("--mask-fov-deg", type=float, default=90.0, help="Render mask field of view.")
    parser.add_argument("--cloud-shell-km", type=float, default=6376.0, help="Cloud shell radius in km.")
    parser.add_argument(
        "--used-time",
        type=str,
        default=None,
        help="UTC timestamp for Himawari local tile runs. If omitted, infer it from the tile directory path.",
    )
    return parser.parse_args()


def _build_geo_area(
    *,
    area_id: str,
    description: str,
    crs: CRS,
    width: int,
    height: int,
    area_extent: tuple[float, float, float, float],
) -> GeoArea:
    return GeoArea(
        area_id=area_id,
        description=description,
        proj_id="geos",
        projection=crs,
        width=width,
        height=height,
        area_extent=area_extent,
    )


def _build_pyresample_area(
    *,
    area_id: str,
    description: str,
    crs: CRS,
    width: int,
    height: int,
    area_extent: tuple[float, float, float, float],
) -> Any:
    if PyresampleAreaDefinition is None:
        raise RuntimeError("pyresample is not installed in this environment")
    return PyresampleAreaDefinition(
        area_id=area_id,
        description=description,
        proj_id="geos",
        projection=crs,
        width=width,
        height=height,
        area_extent=area_extent,
    )


def _goes_area_params(ds: xr.Dataset) -> tuple[CRS, int, int, tuple[float, float, float, float]]:
    proj_var = ds[GOES_GRID_VAR]
    crs = CRS.from_cf(dict(proj_var.attrs))
    scale = float(proj_var.attrs["perspective_point_height"])
    x = np.asarray(ds.x.values, dtype=np.float64)
    y = np.asarray(ds.y.values, dtype=np.float64)
    x_step = float(x[1] - x[0])
    y_step = float(y[1] - y[0])
    area_extent = (
        float(x[0] * scale),
        float((y[0] + y.size * y_step) * scale),
        float((x[0] + x.size * x_step) * scale),
        float(y[0] * scale),
    )
    return crs, int(x.size), int(y.size), area_extent


def _hima_area_params(ds: xr.Dataset) -> tuple[CRS, int, int, tuple[float, float, float, float]]:
    proj_var = ds[HIMA_GRID_VAR]
    crs = CRS.from_cf(dict(proj_var.attrs))
    x = np.asarray(ds.x.values, dtype=np.float64)
    y = np.asarray(ds.y.values, dtype=np.float64)
    geos_scale = float(proj_var.attrs["perspective_point_height"]) * 1e-6
    x_step = float(x[1] - x[0])
    y_step = float(y[1] - y[0])
    area_extent = (
        float(x[0] * geos_scale),
        float((y[0] + y.size * y_step) * geos_scale),
        float((x[0] + x.size * x_step) * geos_scale),
        float(y[0] * geos_scale),
    )
    return crs, int(x.size), int(y.size), area_extent


def _load_goes_pair(path: Path) -> tuple[xr.DataArray, xr.DataArray]:
    current = load_cmi_with_area(path)
    with xr.open_dataset(path) as ds:
        if GOES_DATA_VAR not in ds.variables:
            raise ValueError(f"{path.name} does not contain {GOES_DATA_VAR}")
        crs, width, height, area_extent = _goes_area_params(ds)
        legacy_area = _build_pyresample_area(
            area_id="goes_abi_cmi_c13",
            description="GOES ABI CMIPF C13",
            crs=crs,
            width=width,
            height=height,
            area_extent=area_extent,
        )
        legacy = ds[GOES_DATA_VAR].astype(np.float32).compute()
    legacy.attrs = dict(legacy.attrs)
    legacy.attrs["area"] = legacy_area
    return current, legacy


def _infer_hima_used_time_from_path(path: Path) -> dt.datetime | None:
    parts = path.expanduser().resolve().parts
    for i in range(len(parts) - 4):
        if parts[i] != "AHI-L2-FLDK-ISatSS":
            continue
        year_s, month_s, day_s, hm_s = parts[i + 1 : i + 5]
        if len(hm_s) != 4 or not hm_s.isdigit():
            continue
        try:
            year = int(year_s)
            month = int(month_s)
            day = int(day_s)
            hour = int(hm_s[:2])
            minute = int(hm_s[2:])
            return dt.datetime(year, month, day, hour, minute, tzinfo=dt.timezone.utc)
        except ValueError:
            continue
    return None


def _parse_used_time(value: str | None, *, source: str, input_path: Path) -> dt.datetime:
    if value is None:
        if source == "hima":
            inferred = _infer_hima_used_time_from_path(input_path)
            if inferred is not None:
                return inferred
            raise ValueError(
                "Could not infer Himawari used time from path. Pass --used-time explicitly or use a directory path containing AHI-L2-FLDK-ISatSS/YYYY/MM/DD/HHMM."
            )
        return dt.datetime.now(dt.timezone.utc)
    parsed = dt.datetime.fromisoformat(value)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed.astimezone(dt.timezone.utc)


def _load_hima_pair(
    tile_dir: Path,
    *,
    used_time: dt.datetime,
    observer_lat: float,
    observer_lon: float,
    cloud_shell_km: float,
) -> tuple[xr.DataArray, xr.DataArray]:
    provider = HimaProvider(CloudDiscConfig(cache_dir=tile_dir / ".compare-cache"))
    current, _used_time, src_paths = provider.fetch_bt_c13_from_local_dir(
        tile_dir,
        used_time=used_time,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
        cloud_shell_km=cloud_shell_km,
    )

    stitched = stitch_tiles_from_paths(src_paths, source_label=str(tile_dir))
    try:
        crs, width, height, area_extent = _hima_area_params(stitched)
        legacy_area = _build_pyresample_area(
            area_id="himawari_isatss_m1c13",
            description="Himawari ISatSS stitched M1C13",
            crs=crs,
            width=width,
            height=height,
            area_extent=area_extent,
        )
        legacy = stitched[HIMA_DATA_VAR].astype(np.float32).copy(deep=False)
        legacy.attrs = dict(legacy.attrs)
        legacy.attrs["area"] = legacy_area
        legacy.attrs["source_key_count"] = len(src_paths)
        legacy.attrs["observer_lat"] = float(observer_lat)
        legacy.attrs["observer_lon"] = float(observer_lon)
        return current, legacy
    finally:
        stitched.close()


def _render_metrics(
    da: xr.DataArray,
    *,
    lat: float,
    lon: float,
    alt: float,
    az: float,
    radius_px: int,
    edge_fov_deg: float,
    mask_fov_deg: float,
    cloud_shell_km: float,
    label: str,
) -> tuple[Image.Image, np.ndarray, np.ndarray, RenderStats]:
    sampler = build_bt_sampler(da)
    lon_grid, lat_grid, mask_inside = az_project_lonlat_grid(
        lat0_deg=lat,
        lon0_deg=lon,
        alt0_deg=alt,
        az0_deg=az,
        radius_px=radius_px + 1,
        cloud_shell_km=cloud_shell_km,
        alt_min_deg=0.0,
        edge_fov_deg=edge_fov_deg,
        mask_fov_deg=mask_fov_deg,
    )
    bt = sampler(lon_grid, lat_grid)
    finite_bt = np.isfinite(bt)
    inside_count = int(np.count_nonzero(mask_inside))
    valid_inside = int(np.count_nonzero(mask_inside & finite_bt))
    coverage_ratio = float(valid_inside) / float(inside_count) if inside_count else 1.0
    bt_warm, eq_samples = estimate_bt_warm_from_equator_band(
        da,
        lon_center_deg=lon,
        delta_lon=60.0,
        equator_lat=0.0,
        warm_p=97.0,
        half=5,
    )
    bt_cold = estimate_bt_cold_hybrid(bt, mask_inside, eq_samples, bt_warm, cold_local_p=5.0, cold_eq_p=3.0)
    img = convert_bt_to_la_image(bt, mask_inside, bt_warm, bt_cold)
    inside_vals = bt[mask_inside & finite_bt].astype(np.float64)
    if inside_vals.size:
        bt_min = float(np.min(inside_vals))
        bt_max = float(np.max(inside_vals))
        bt_mean = float(np.mean(inside_vals))
    else:
        bt_min = float("nan")
        bt_max = float("nan")
        bt_mean = float("nan")
    stats = RenderStats(
        label=label,
        bt_warm=float(bt_warm),
        bt_cold=float(bt_cold),
        coverage_ratio=coverage_ratio,
        inside_pixels=inside_count,
        inside_valid_pixels=valid_inside,
        eq_samples=int(eq_samples.size),
        bt_min=bt_min,
        bt_max=bt_max,
        bt_mean=bt_mean,
    )
    return img, bt, mask_inside, stats


def _save_diff_image(a: np.ndarray, b: np.ndarray, dst: Path) -> dict[str, float]:
    diff = np.abs(a.astype(np.int16) - b.astype(np.int16)).astype(np.uint8)
    if diff.ndim == 2:
        Image.fromarray(diff, mode="L").save(dst)
    else:
        Image.fromarray(diff, mode="LA").save(dst)
    return {
        "max_abs_diff": float(np.max(diff)),
        "mean_abs_diff": float(np.mean(diff)),
        "nonzero_pixels": int(np.count_nonzero(np.any(diff != 0, axis=-1) if diff.ndim == 3 else diff != 0)),
    }


def main() -> int:
    args = _parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if PyresampleAreaDefinition is None:
        print("pyresample is not installed. This script can render the current path only after pyresample is installed.")
        return 2

    if args.source == "goes":
        current_da, legacy_da = _load_goes_pair(args.input)
    else:
        current_da, legacy_da = _load_hima_pair(
            args.input,
            used_time=_parse_used_time(args.used_time, source=args.source, input_path=args.input),
            observer_lat=float(args.lat),
            observer_lon=float(args.lon),
            cloud_shell_km=float(args.cloud_shell_km),
        )

    current_img, current_bt, current_mask, current_stats = _render_metrics(
        current_da,
        lat=float(args.lat),
        lon=float(args.lon),
        alt=float(args.alt),
        az=float(args.az),
        radius_px=int(args.radius_px),
        edge_fov_deg=float(args.edge_fov_deg),
        mask_fov_deg=float(args.mask_fov_deg),
        cloud_shell_km=float(args.cloud_shell_km),
        label="current",
    )
    legacy_img, legacy_bt, legacy_mask, legacy_stats = _render_metrics(
        legacy_da,
        lat=float(args.lat),
        lon=float(args.lon),
        alt=float(args.alt),
        az=float(args.az),
        radius_px=int(args.radius_px),
        edge_fov_deg=float(args.edge_fov_deg),
        mask_fov_deg=float(args.mask_fov_deg),
        cloud_shell_km=float(args.cloud_shell_km),
        label="pyresample",
    )

    current_img.save(out_dir / "current.png")
    legacy_img.save(out_dir / "pyresample.png")

    current_arr = np.asarray(current_img, dtype=np.uint8)
    legacy_arr = np.asarray(legacy_img, dtype=np.uint8)
    image_diff = _save_diff_image(current_arr, legacy_arr, out_dir / "image_diff.png")

    bt_diff = np.abs(current_bt - legacy_bt)
    bt_diff = bt_diff[np.isfinite(bt_diff)]
    mask_mismatch = int(np.count_nonzero(current_mask ^ legacy_mask))

    summary = {
        "source": args.source,
        "input": str(args.input),
        "current": asdict(current_stats),
        "pyresample": asdict(legacy_stats),
        "delta": {
            "bt_warm": current_stats.bt_warm - legacy_stats.bt_warm,
            "bt_cold": current_stats.bt_cold - legacy_stats.bt_cold,
            "coverage_ratio": current_stats.coverage_ratio - legacy_stats.coverage_ratio,
            "eq_samples": current_stats.eq_samples - legacy_stats.eq_samples,
            "mask_mismatch_pixels": mask_mismatch,
            "bt_max_abs_diff": float(np.max(bt_diff)) if bt_diff.size else 0.0,
            "bt_mean_abs_diff": float(np.mean(bt_diff)) if bt_diff.size else 0.0,
            "image_max_abs_diff": image_diff["max_abs_diff"],
            "image_mean_abs_diff": image_diff["mean_abs_diff"],
            "image_nonzero_pixels": image_diff["nonzero_pixels"],
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
