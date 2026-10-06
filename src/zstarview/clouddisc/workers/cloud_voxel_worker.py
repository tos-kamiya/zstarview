"""One-shot subprocess renderer for native-pixel cloud voxels."""

from __future__ import annotations

import argparse
import json
import os
import pickle
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

import numpy as np

from ...geosatellite.types import GeoSatelliteVoxelSource
from ...render.cloud_voxels import (
    CLOUD_AMOUNT_SUBTRACTION,
    shade_geo_satellite_voxels,
    shade_native_voxels,
)
from ...render.ground_mask import inverse_project_disc
from ...types import ScreenGeometry

PROTOCOL_VERSION = 4


def _worker_main(input_path: Path, output_path: Path, result_path: Path) -> int:
    try:
        with input_path.open("rb") as fp:
            request = pickle.load(fp)
        if request.get("protocol") != PROTOCOL_VERSION:
            raise ValueError("unsupported cloud voxel protocol")
        width = int(request["width"])
        height = int(request["height"])
        radius = int(request["radius_px"])
        center = (radius, radius)
        geometry = ScreenGeometry(center, radius)
        altitudes, azimuths, inside = inverse_project_disc(
            width,
            height,
            geometry,
            tuple(request["view_center"]),
            # The image radius covers content_fov, including the overscan
            # beyond the screen geometry's edge radius.
            edge_fov_deg=float(request["content_fov_deg"]),
            content_fov_deg=float(request["content_fov_deg"]),
        )
        base = np.zeros((int(np.count_nonzero(inside)), 3), dtype=np.float32)
        source = request["source"]
        shade = (
            shade_geo_satellite_voxels
            if isinstance(source, GeoSatelliteVoxelSource)
            else shade_native_voxels
        )
        shaded, info, transmission = shade(
            source,
            float(request["lat"]),
            float(request["lon"]),
            altitudes,
            azimuths,
            float(request["sun_alt_deg"]),
            float(request["sun_az_deg"]),
            base,
            sunlight_mix=float(request["sunlight_mix"]),
            night_color_rgb=tuple(request["night_color_rgb"]),
            opacity=1.0,
            height_layer_transform=not isinstance(source, GeoSatelliteVoxelSource),
            return_transmission=True,
            cloud_amount_subtract=CLOUD_AMOUNT_SUBTRACTION,
        )
        rgba = np.zeros((height, width, 4), dtype=np.uint8)
        alpha = np.clip(1.0 - transmission, 0.0, 1.0)
        straight_rgb = np.zeros_like(shaded)
        nonzero = alpha > 1.0e-6
        straight_rgb[nonzero] = shaded[nonzero] / alpha[nonzero, None]
        rgba[inside, :3] = np.rint(
            np.clip(straight_rgb, 0.0, 1.0) * 255.0
        ).astype(np.uint8)
        rgba[inside, 3] = np.rint(alpha * 255.0).astype(np.uint8)
        with output_path.open("wb") as fp:
            np.save(fp, rgba, allow_pickle=False)
        result_path.write_text(
            json.dumps(
                {
                    "protocol": PROTOCOL_VERSION,
                    "request_id": int(request["request_id"]),
                    "source_key": str(request.get("source_key", "")),
                    "status": "ok",
                    "shape": list(rgba.shape),
                    "dtype": str(rgba.dtype),
                    "voxel_info": info,
                },
                ensure_ascii=True,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        return 0
    except Exception as exc:
        result_path.write_text(
            json.dumps(
                {
                    "protocol": PROTOCOL_VERSION,
                    "status": "error",
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "traceback": traceback.format_exc(),
                },
                ensure_ascii=True,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        return 1


def render_cloud_voxels_in_subprocess(
    source,
    *,
    lat: float,
    lon: float,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    content_fov_deg: float,
    radius_px: int,
    sun_alt_deg: float,
    sun_az_deg: float,
    sunlight_mix: float,
    night_color_rgb: tuple[float, float, float],
    request_id: int,
    timeout_s: float = 120.0,
) -> np.ndarray:
    """Render the content disc; radius_px is its image-space radius."""
    radius = max(1, int(radius_px))
    input_data = {
        "protocol": PROTOCOL_VERSION,
        "request_id": int(request_id),
        "source_key": str(getattr(source, "source_key", "")),
        "source": source,
        "lat": float(lat),
        "lon": float(lon),
        "view_center": (float(view_center[0]), float(view_center[1])),
        "edge_fov_deg": float(edge_fov_deg),
        "content_fov_deg": max(float(edge_fov_deg), float(content_fov_deg)),
        "radius_px": radius,
        "width": radius * 2 + 1,
        "height": radius * 2 + 1,
        "sun_alt_deg": float(sun_alt_deg),
        "sun_az_deg": float(sun_az_deg),
        "sunlight_mix": float(sunlight_mix),
        "night_color_rgb": tuple(float(value) for value in night_color_rgb),
    }
    with tempfile.TemporaryDirectory(prefix="zstarview-cloud-voxel-") as temp_dir:
        work = Path(temp_dir)
        input_path = work / "request.pkl"
        output_path = work / "cloud.npy"
        result_path = work / "result.json"
        with input_path.open("wb") as fp:
            pickle.dump(input_data, fp, protocol=pickle.HIGHEST_PROTOCOL)
        env = os.environ.copy()
        root = str(Path(__file__).resolve().parents[3])
        old_pythonpath = env.get("PYTHONPATH")
        env["PYTHONPATH"] = root if not old_pythonpath else root + os.pathsep + old_pythonpath
        command = [
            sys.executable,
            "-m",
            "zstarview.clouddisc.workers.cloud_voxel_worker",
            "--worker",
            str(input_path),
            str(output_path),
            str(result_path),
        ]
        try:
            completed = subprocess.run(
                command,
                check=False,
                timeout=max(0.1, float(timeout_s)),
                capture_output=True,
                text=True,
                env=env,
            )
        except subprocess.TimeoutExpired as exc:
            raise TimeoutError("cloud voxel worker timed out") from exc
        if not result_path.is_file():
            raise RuntimeError(
                "cloud voxel worker exited without a result manifest: "
                + completed.stderr[-1000:]
            )
        manifest = json.loads(result_path.read_text(encoding="utf-8"))
        if manifest.get("protocol") != PROTOCOL_VERSION:
            raise RuntimeError("cloud voxel worker protocol mismatch")
        if manifest.get("status") != "ok":
            raise RuntimeError(
                "cloud voxel worker failed: "
                + str(manifest.get("error", completed.stderr[-1000:]))
            )
        if int(manifest.get("request_id", -1)) != int(request_id):
            raise RuntimeError("cloud voxel worker returned a stale request")
        if str(manifest.get("source_key", "")) != input_data["source_key"]:
            raise RuntimeError("cloud voxel worker returned a different source")
        with output_path.open("rb") as fp:
            rgba = np.load(fp, allow_pickle=False)
        expected_shape = (radius * 2 + 1, radius * 2 + 1, 4)
        if rgba.shape != expected_shape or rgba.dtype != np.uint8:
            raise RuntimeError("cloud voxel worker returned an invalid image")
        return rgba


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("result", type=Path)
    args = parser.parse_args()
    return _worker_main(args.input, args.output, args.result)


if __name__ == "__main__":
    raise SystemExit(main())
