from __future__ import annotations

import threading
from concurrent.futures import Future, ThreadPoolExecutor

from PySide6.QtCore import QObject, Signal

from ..road_night_lights import (
    ROAD_NIGHT_LIGHT_MAX_DISTANCE_KM,
    clip_road_night_lights_to_annulus,
    load_or_fetch_road_night_lights,
    load_road_night_lights_cache,
    project_road_night_lights,
    road_night_lights_scope_key,
    simplify_road_night_light_way_for_observer,
)
from ..types import ViewerData


class RoadNightLightsController(QObject):
    """Fetch and project road geometry without blocking the GUI thread."""

    road_started = Signal(object)
    road_ready = Signal(object)
    road_failed = Signal(object)

    def __init__(self, *, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._executor = ThreadPoolExecutor(
            max_workers=1, thread_name_prefix="road-lights"
        )
        self._abort_event = threading.Event()
        self._future: Future[None] | None = None
        self._key: tuple[float, float] | None = None
        self._stopping = False

    def has_in_flight_update(self) -> bool:
        return self._future is not None and not self._future.done()

    def shutdown(self) -> None:
        self._stopping = True
        self._abort_event.set()
        self._executor.shutdown(wait=False, cancel_futures=True)

    def update(self, *, viewer_data: ViewerData, reason: str = "manual") -> bool:
        if self._stopping or self.has_in_flight_update():
            return False
        key = (
            round(float(viewer_data.lat_deg), 4),
            round(float(viewer_data.lon_deg), 4),
        )
        if self._key == key:
            return False
        self._key = key
        self._abort_event.clear()
        self.road_started.emit(
            {"banner": "Road night lights: loading...", "reason": reason}
        )
        self._future = self._executor.submit(self._run, viewer_data)
        self._future.add_done_callback(self._finished)
        return True

    def _run(self, viewer_data: ViewerData) -> None:
        scope_key = road_night_lights_scope_key(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            radius_km=ROAD_NIGHT_LIGHT_MAX_DISTANCE_KM,
        )
        cache_hit = load_road_night_lights_cache(scope_key) is not None
        snapshot = load_or_fetch_road_night_lights(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            radius_km=ROAD_NIGHT_LIGHT_MAX_DISTANCE_KM,
            abort_event=self._abort_event,
        )
        simplified = tuple(
            simplify_road_night_light_way_for_observer(
                way,
                observer_lat_deg=float(viewer_data.lat_deg),
                observer_lon_deg=float(viewer_data.lon_deg),
            )
            for way in snapshot.ways
        )
        clipped = clip_road_night_lights_to_annulus(
            simplified,
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
        )
        polylines = project_road_night_lights(
            clipped,
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            observer_height_m=float(viewer_data.observer_height_m),
        )
        self.road_ready.emit(
            {
                "polylines": list(polylines),
                "source": "Road: cache" if cache_hit else "Road: API",
            }
        )

    def _finished(self, future: Future[None]) -> None:
        self._future = None
        if self._stopping or self._abort_event.is_set():
            return
        try:
            future.result()
        except (
            OSError,
            RuntimeError,
            TypeError,
            ValueError,
        ) as exc:  # pragma: no cover - Qt callback boundary
            self.road_failed.emit({"banner": f"Road night lights: unavailable ({exc})"})
