#!/usr/bin/env python3
"""Benchmark the low-level star renderer and its renderer-owned image cache.

The script measures `zstarview.render.stars.draw_stars_normal/fast` directly.
Catalog loading, star position calculation, and input construction are done
before the timed region.  It does not exercise the GUI-level star-surface or
frame caches; see `docs/design.md` for how those wrap this renderer.

Scenarios (each is timed separately):

- ``hit``:   identical inputs with a warm `StarRenderCache`.
- ``miss``:  view azimuth advances by 0.1 deg per call, so every call
             regenerates the image (warm, after one untimed call).
- ``pos``:   explicit `screen_positions`; the image cache is never reused.
- ``split``: scenic split sequence (faint pass on the low-res surface,
             then bright pass at viewport size) sharing one cache, with
             identical inputs.  Shows whether the two passes evict each other.

Examples:
  uv run -p .venv/bin/python dev-samples/benchmark_star_render.py
  uv run -p .venv/bin/python dev-samples/benchmark_star_render.py \
      --sizes 960x720,2560x1440 --vmag 7,9 --repeat 3 --json out.json
  uv run -p .venv/bin/python dev-samples/benchmark_star_render.py \
      --dump-images base.npz      # pixel snapshots for regression comparison
  uv run -p .venv/bin/python dev-samples/benchmark_star_render.py \
      --compare-images base.npz   # compare current output with a snapshot
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import time
import tracemalloc
from dataclasses import dataclass

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import astropy.time  # noqa: E402
import numpy as np  # noqa: E402
from PySide6.QtGui import QGuiApplication, QImage, QPainter  # noqa: E402

from zstarview.astro import (  # noqa: E402
    calculate_visible_stars,
    prepare_star_catalog_arrays,
)
from zstarview.catalog import load_star_catalog  # noqa: E402
from zstarview.paths import STARS_CSV_FILE  # noqa: E402
from zstarview.render import stars as render_stars  # noqa: E402
from zstarview.render.geometry import (  # noqa: E402
    _altaz_to_normalized_xy_vectorized,
    _normalized_to_screen_xy_vectorized,
    get_screen_geometry,
)
from zstarview.render.qt_image import qimage_to_np_rgba  # noqa: E402
from zstarview.types import CelestialData, ScreenGeometry, ViewerData  # noqa: E402

FIXED_TIME_UTC = "2026-01-15T12:00:00"
FIXED_LOCATION = (35.47, 133.05)
VIEW_CENTER = (45.0, 180.0)
EDGE_FOV_DEG = 95.0
CONTENT_FOV_DEG = 115.0
STAR_BASE_RADIUS = 3.0
SCENIC_SURFACE_SCALE = 0.75


@dataclass(frozen=True)
class Case:
    width: int
    height: int
    vmag: float
    fast: bool
    light: bool
    outline: bool

    @property
    def name(self) -> str:
        return (
            f"{self.width}x{self.height}-v{self.vmag:g}"
            f"-{'fast' if self.fast else 'normal'}"
            f"-{'light' if self.light else 'dark'}"
            f"{'-outline' if self.outline else ''}"
        )


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sizes", default="960x720,2560x1440")
    parser.add_argument("--vmag", default="7,9")
    parser.add_argument("--modes", default="normal,fast", help="normal,fast")
    parser.add_argument("--backgrounds", default="dark,light", help="dark,light")
    parser.add_argument("--outline", action="store_true", help="also run outline bright-body mode")
    parser.add_argument("--scenarios", default="hit,miss,pos,split")
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument("--warmup", type=int, default=3)
    parser.add_argument("--repeat", type=int, default=2, help="rounds; case order is reversed on odd rounds")
    parser.add_argument("--json", help="write raw samples and summary to this JSON file")
    parser.add_argument("--dump-images", help="write pixel snapshots (npz) and exit")
    parser.add_argument("--compare-images", help="compare pixel snapshots (npz) and exit")
    return parser.parse_args(argv)


_STAR_CACHE: dict[float, tuple[dict, int]] = {}


def _stars_for(vmag: float) -> tuple[dict, int]:
    if vmag not in _STAR_CACHE:
        df = load_star_catalog(STARS_CSV_FILE, vmag_threshold=vmag)
        arrays = prepare_star_catalog_arrays(df)
        stars, _ = calculate_visible_stars(
            arrays,
            FIXED_LOCATION[0],
            FIXED_LOCATION[1],
            0.0,
            astropy.time.Time(FIXED_TIME_UTC, scale="utc"),
        )
        _STAR_CACHE[vmag] = (stars, len(df))
    return _STAR_CACHE[vmag]


def _celestial(stars: dict) -> CelestialData:
    return CelestialData(
        time=astropy.time.Time(FIXED_TIME_UTC, scale="utc"),
        planets=[],
        stars=stars,  # type: ignore[arg-type]
        deep_sky_objects={},  # type: ignore[arg-type]
        celestial_equator_points=[],
        ecliptic_points=[],
        horizon_points=[],
    )


def _viewer(az: float = VIEW_CENTER[1]) -> ViewerData:
    return ViewerData(
        location=FIXED_LOCATION,
        timezone_name="UTC",
        city_name="bench",
        view_center=(VIEW_CENTER[0], az),
        edge_fov_deg=EDGE_FOV_DEG,
        content_fov_deg=CONTENT_FOV_DEG,
    )


def _geometry(width: int, height: int) -> ScreenGeometry:
    return get_screen_geometry(
        width,
        height,
        VIEW_CENTER[0],
        edge_fov_deg=EDGE_FOV_DEG,
        content_fov_deg=CONTENT_FOV_DEG,
    )


def _draw(
    case: Case,
    image: QImage,
    geometry: ScreenGeometry,
    data: CelestialData,
    viewer: ViewerData,
    cache: render_stars.StarRenderCache | None,
    *,
    size: tuple[int, int],
    vmag_limit: float,
    vmag_min_exclusive: float | None = None,
    screen_positions: np.ndarray | None = None,
) -> None:
    painter = QPainter(image)
    try:
        draw = render_stars.draw_stars_fast if case.fast else render_stars.draw_stars_normal
        draw(
            painter,
            geometry,
            data,
            viewer,
            STAR_BASE_RADIUS,
            outline_bright_bodies=case.outline,
            light_background_outline=case.light,
            draw_vmag_limit=vmag_limit,
            draw_vmag_min_exclusive=vmag_min_exclusive,
            viewport_size=size,
            screen_positions=screen_positions,
            render_cache=cache,
        )
    finally:
        painter.end()


def _new_image(width: int, height: int) -> QImage:
    image = QImage(width, height, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    return image


def _screen_positions(stars: dict, viewer: ViewerData, geometry: ScreenGeometry, vmag: float) -> np.ndarray:
    mask = stars["vmag"] <= vmag
    nx, ny = _altaz_to_normalized_xy_vectorized(
        stars["alt"][mask], stars["az"][mask], viewer.view_center, edge_fov_deg=float(viewer.edge_fov_deg)
    )
    x, y = _normalized_to_screen_xy_vectorized(nx, ny, geometry)
    return np.column_stack([x, y])


def _run_scenario(case: Case, scenario: str, iterations: int, warmup: int) -> dict[str, object]:
    stars, catalog_rows = _stars_for(case.vmag)
    data = _celestial(stars)
    geometry = _geometry(case.width, case.height)
    size = (case.width, case.height)
    image = _new_image(case.width, case.height)
    cache = render_stars.StarRenderCache()
    viewer = _viewer()

    if scenario == "hit":
        def step(i: int) -> None:
            _draw(case, image, geometry, data, viewer, cache, size=size, vmag_limit=case.vmag)
    elif scenario == "miss":
        viewers = [_viewer(VIEW_CENTER[1] + 0.1 * (i + 1)) for i in range(iterations + warmup + 2)]

        def step(i: int) -> None:
            _draw(case, image, geometry, data, viewers[i], cache, size=size, vmag_limit=case.vmag)
    elif scenario == "pos":
        positions = _screen_positions(stars, viewer, geometry, case.vmag)

        def step(i: int) -> None:
            _draw(
                case, image, geometry, data, viewer, cache,
                size=size, vmag_limit=case.vmag, screen_positions=positions,
            )
    elif scenario == "split":
        low_w = max(1, int(case.width * SCENIC_SURFACE_SCALE))
        low_h = max(1, int(case.height * SCENIC_SURFACE_SCALE))
        low_image = _new_image(low_w, low_h)
        low_geometry = ScreenGeometry(
            center=(
                int(round(geometry.center[0] * SCENIC_SURFACE_SCALE)),
                int(round(geometry.center[1] * SCENIC_SURFACE_SCALE)),
            ),
            radius=max(1, int(round(geometry.radius * SCENIC_SURFACE_SCALE))),
        )

        def step(i: int) -> None:
            _draw(
                case, low_image, low_geometry, data, viewer, cache,
                size=(low_w, low_h), vmag_limit=case.vmag, vmag_min_exclusive=4.0,
            )
            _draw(case, image, geometry, data, viewer, cache, size=size, vmag_limit=4.0)
    else:
        raise ValueError(scenario)

    # Prime once (first generation is excluded from warm timings).
    first_ns = time.perf_counter_ns()
    step(0)
    first_ns = time.perf_counter_ns() - first_ns
    for i in range(warmup):
        step(i + 1)
    samples: list[int] = []
    cpu_samples: list[int] = []
    base = 1 + warmup
    for i in range(iterations):
        t0 = time.perf_counter_ns()
        c0 = time.process_time_ns()
        step(base + i)
        cpu_samples.append(time.process_time_ns() - c0)
        samples.append(time.perf_counter_ns() - t0)
    # tracemalloc slows Python-level loops heavily, so the allocation peak is
    # measured in one separate, untimed call (NumPy reports its buffers).
    tracemalloc.start()
    tracemalloc.reset_peak()
    step(base + iterations)
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    rendered = int(np.count_nonzero(stars["vmag"] <= case.vmag))
    return {
        "case": case.name,
        "scenario": scenario,
        "catalog_rows": catalog_rows,
        "input_stars": int(stars["alt"].size),
        "draw_target_stars": rendered,
        "first_ms": first_ns / 1e6,
        "samples_ns": samples,
        "cpu_ns": cpu_samples,
        "tracemalloc_peak_mb": peak / 1e6,
    }


def _summary(rec: dict[str, object]) -> str:
    samples = sorted(rec["samples_ns"])  # type: ignore[arg-type]
    med = statistics.median(samples) / 1e6
    p95 = samples[min(len(samples) - 1, int(round(0.95 * (len(samples) - 1))))] / 1e6
    cpu = statistics.median(rec["cpu_ns"]) / 1e6  # type: ignore[arg-type]
    return (
        f"{rec['case']:<36} {rec['scenario']:<5} n={len(samples):<3} "
        f"med={med:8.3f}ms p95={p95:8.3f}ms cpu={cpu:8.3f}ms first={rec['first_ms']:8.2f}ms "
        f"peak={rec['tracemalloc_peak_mb']:7.1f}MB stars={rec['draw_target_stars']}/{rec['input_stars']}"
    )


def _cases(args: argparse.Namespace) -> list[Case]:
    sizes = [tuple(int(v) for v in s.split("x")) for s in args.sizes.split(",") if s]
    vmags = [float(v) for v in args.vmag.split(",") if v]
    modes = [m == "fast" for m in args.modes.split(",") if m]
    backgrounds = [b == "light" for b in args.backgrounds.split(",") if b]
    outlines = [False, True] if args.outline else [False]
    return [
        Case(w, h, v, f, lt, o)
        for (w, h), v, f, lt, o in itertools.product(sizes, vmags, modes, backgrounds, outlines)
    ]


def _snapshot_cases() -> list[Case]:
    base = [Case(960, 720, 9.0, fast, light, outline)
            for fast in (False, True) for light in (False, True) for outline in (False, True)]
    return base + [Case(640, 960, 7.0, False, False, False)]


def _snapshot(case: Case) -> dict[str, np.ndarray]:
    """Render a set of variants, including cache-hit and twinkle paths."""
    out: dict[str, np.ndarray] = {}
    stars, _ = _stars_for(case.vmag)
    data = _celestial(stars)
    geometry = _geometry(case.width, case.height)
    viewer = _viewer()
    size = (case.width, case.height)
    cache = render_stars.StarRenderCache()
    for label, kwargs in (
        ("plain", {}),
        ("plain-hit", {}),
        ("faint", {"vmag_min_exclusive": 4.0}),
        ("bright", {"vmag_limit": 4.0}),
    ):
        image = _new_image(*size)
        merged: dict = {"vmag_limit": case.vmag}
        merged.update(kwargs)
        _draw(case, image, geometry, data, viewer, cache, size=size, **merged)
        out[f"{case.name}/{label}"] = qimage_to_np_rgba(image)
    image = _new_image(*size)
    _draw(
        case, image, geometry, data, viewer, cache, size=size, vmag_limit=case.vmag,
        screen_positions=_screen_positions(stars, viewer, geometry, case.vmag) + 0.3,
    )
    out[f"{case.name}/pos"] = qimage_to_np_rgba(image)
    # Twinkle on a few bright rows.
    rows = np.flatnonzero(stars["vmag"] <= 3.0)[:5]
    targets = tuple((int(stars["star_index"][r]), 0.6) for r in rows)
    image = _new_image(*size)
    painter = QPainter(image)
    (render_stars.draw_stars_fast if case.fast else render_stars.draw_stars_normal)(
        painter, geometry, data, viewer, STAR_BASE_RADIUS,
        outline_bright_bodies=case.outline, light_background_outline=case.light,
        draw_vmag_limit=case.vmag, viewport_size=size, twinkle_targets=targets,
        render_cache=cache,
    )
    painter.end()
    out[f"{case.name}/twinkle"] = qimage_to_np_rgba(image)
    return out


def _git_rev() -> str:
    try:
        rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "src"], capture_output=True, text=True, check=True).stdout.strip()
        return rev + ("+dirty" if dirty else "")
    except Exception:
        return "unknown"


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    _app = QGuiApplication.instance() or QGuiApplication(sys.argv[:1])
    if args.dump_images or args.compare_images:
        snaps: dict[str, np.ndarray] = {}
        for case in _snapshot_cases():
            snaps.update(_snapshot(case))
        if args.dump_images:
            np.savez_compressed(args.dump_images, **snaps)
            print(f"wrote {len(snaps)} snapshots to {args.dump_images}")
            return 0
        ref = np.load(args.compare_images)
        failures = 0
        for key, arr in snaps.items():
            if key not in ref:
                print(f"MISSING {key}")
                failures += 1
                continue
            diff = np.abs(arr.astype(int) - ref[key].astype(int))
            if diff.max() != 0:
                failures += 1
                print(f"DIFF {key}: max={diff.max()} pixels={int(np.count_nonzero(diff.max(axis=2)))}")
        print(f"compared {len(snaps)} snapshots: {failures} differ")
        return 1 if failures else 0

    env = {
        "commit": _git_rev(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "fixed": {
            "time": FIXED_TIME_UTC, "location": FIXED_LOCATION, "view_center": VIEW_CENTER,
            "edge_fov": EDGE_FOV_DEG, "content_fov": CONTENT_FOV_DEG, "star_base_radius": STAR_BASE_RADIUS,
            "dpr": 1.0,
        },
        "argv": sys.argv[1:] if argv is None else argv,
    }
    print(json.dumps(env))
    cases = _cases(args)
    scenarios = [s for s in args.scenarios.split(",") if s]
    for vmag in sorted({c.vmag for c in cases}):
        _stars_for(vmag)
    records: list[dict[str, object]] = []
    for rnd in range(args.repeat):
        order = list(itertools.product(cases, scenarios))
        if rnd % 2 == 1:
            order.reverse()
        for case, scenario in order:
            if scenario == "split" and case.light:
                continue
            rec = _run_scenario(case, scenario, args.iterations, args.warmup)
            rec["round"] = rnd
            records.append(rec)
            print(f"[r{rnd}] " + _summary(rec), flush=True)
    rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    print(f"max RSS: {rss_mb:.1f} MB")
    if args.json:
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump({"env": env, "max_rss_mb": rss_mb, "records": records}, handle)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
