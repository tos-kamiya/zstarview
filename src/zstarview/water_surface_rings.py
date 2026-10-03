"""Connect existing water samples along their scan-distance circles."""

from __future__ import annotations

import math
import threading
from collections import defaultdict
from collections.abc import Sequence

from .water_overlay import (
    DEFAULT_WATER_AZIMUTH_STEP_DEG,
    WaterOverlayPoint,
    WaterOverlayPolyline,
    _cooperative_yield,
)


def build_water_surface_ring_polylines(
    points: Sequence[WaterOverlayPoint],
    *,
    azimuth_step_deg: float = DEFAULT_WATER_AZIMUTH_STEP_DEG,
    abort_event: threading.Event | None = None,
) -> tuple[WaterOverlayPolyline, ...]:
    """Join adjacent wet samples without crossing missing azimuth samples.

    Use the full dot input, before display thinning. Preserve sample heights,
    distances and occlusion information by reusing the original point objects.
    Samples without a scan distance cannot identify a circle and are omitted.
    """
    if not math.isfinite(azimuth_step_deg) or azimuth_step_deg <= 0:
        raise ValueError("azimuth_step_deg must be finite and positive")
    groups: dict[tuple[float, str, str], list[WaterOverlayPoint]] = defaultdict(list)
    for index, point in enumerate(points):
        _cooperative_yield(abort_event, interval=128, iteration_index=index)
        radius_m = point.scan_distance_m
        if radius_m is None or not math.isfinite(radius_m) or radius_m <= 0:
            continue
        if not math.isfinite(point.az_deg) or not math.isfinite(point.alt_deg):
            continue
        groups[(round(radius_m, 6), point.water_category, point.water_id)].append(point)

    lines: list[WaterOverlayPolyline] = []
    max_gap_deg = azimuth_step_deg * 1.5
    for index, ((radius, category, water_id), samples) in enumerate(
        sorted(groups.items())
    ):
        _cooperative_yield(abort_event, interval=16, iteration_index=index)
        ordered = sorted(samples, key=lambda point: point.az_deg % 360)
        # Overlapping inputs may include the same azimuth sample twice.
        ordered = list(
            {round(point.az_deg % 360, 8): point for point in ordered}.values()
        )
        run: list[WaterOverlayPoint] = []

        def append_run() -> None:
            if len(run) >= 2:
                lines.append(
                    WaterOverlayPolyline(
                        water_id=f"water-ring-{radius:g}-{water_id}",
                        water_category=category,
                        points=tuple(run),
                    )
                )
            run.clear()

        for point in ordered:
            if run:
                previous = run[-1]
                gap = (point.az_deg - previous.az_deg) % 360
                missing_index = (
                    previous.scan_azimuth_index is not None
                    and point.scan_azimuth_index is not None
                    and point.scan_azimuth_index != previous.scan_azimuth_index + 1
                )
                if gap > max_gap_deg or missing_index:
                    append_run()
            run.append(point)
        append_run()
        if len(ordered) >= 2:
            seam_gap = (ordered[0].az_deg - ordered[-1].az_deg) % 360
            if 0 < seam_gap <= max_gap_deg:
                lines.append(
                    WaterOverlayPolyline(
                        water_id=f"water-ring-{radius:g}-{water_id}-seam",
                        water_category=category,
                        points=(ordered[-1], ordered[0]),
                    )
                )
    return tuple(lines)
