# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
import threading
from pathlib import Path
from typing import Optional

from PySide6.QtCore import QObject, Signal

from ..data.import_overture_buildings import (
    DEFAULT_FETCH_RADIUS_KM,
    derive_dataset_name,
    import_overture_buildings,
)
from ..urban_debug_layer import resolve_urban_debug_layer_for_viewer
from ..types import ViewerData

logger = logging.getLogger(__name__)


class UrbanOutlineController(QObject):
    urban_started = Signal(object)
    urban_ready = Signal(object)
    urban_failed = Signal(object)

    def __init__(
        self,
        *,
        derived_root_dir: Path,
        min_building_height_m: float = 0.0,
        radius_km: float = DEFAULT_FETCH_RADIUS_KM,
        feature_type: str = "building",
        overturemaps_bin: str = "overturemaps",
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._derived_root_dir = Path(derived_root_dir)
        self._min_building_height_m = float(min_building_height_m)
        self._radius_km = float(radius_km)
        self._feature_type = str(feature_type)
        self._overturemaps_bin = str(overturemaps_bin)
        self._running = False
        self._stopping = False
        self._completed_key: Optional[str] = None
        self._lock = threading.Lock()

    def shutdown(self) -> None:
        with self._lock:
            self._stopping = True

    def update(
        self,
        *,
        viewer_data: ViewerData,
        reason: str = "manual",
    ) -> bool:
        dataset_name = derive_dataset_name(
            float(viewer_data.lat_deg),
            float(viewer_data.lon_deg),
            self._radius_km,
            self._feature_type,
            self._min_building_height_m,
        )
        with self._lock:
            if self._stopping or self._running:
                return False
            if self._completed_key == dataset_name:
                return False
            self._running = True

        derived_dir = self._derived_root_dir / dataset_name / "bldg"
        if derived_dir.exists():
            self.urban_started.emit({"banner": "Urban outline: loading cache..."})
        else:
            self.urban_started.emit({"banner": "Urban outline: downloading..."})

        worker = threading.Thread(
            target=self._run_update,
            kwargs={
                "viewer_data": viewer_data,
                "dataset_name": dataset_name,
                "reason": reason,
            },
            daemon=True,
        )
        worker.start()
        return True

    def _run_update(
        self,
        *,
        viewer_data: ViewerData,
        dataset_name: str,
        reason: str,
    ) -> None:
        try:
            derived_dir = self._derived_root_dir / dataset_name / "bldg"
            source = "cache"
            if not derived_dir.exists():
                if reason == "initial":
                    logger.info("Downloading initial urban-outline data from Overture...")
                else:
                    logger.info("Downloading urban-outline data from Overture...")
                import_overture_buildings(
                    lat_deg=float(viewer_data.lat_deg),
                    lon_deg=float(viewer_data.lon_deg),
                    radius_km=self._radius_km,
                    derived_root_dir=self._derived_root_dir,
                    min_building_height_m=self._min_building_height_m,
                    feature_type=self._feature_type,
                    fmt="geojsonseq",
                    overturemaps_bin=self._overturemaps_bin,
                    dataset_name=dataset_name,
                    keep_download=None,
                    no_stac=False,
                )
                source = "overture"

            outlines = resolve_urban_debug_layer_for_viewer(
                viewer_data,
                derived_root_dir=self._derived_root_dir,
            )
            with self._lock:
                if not self._stopping:
                    self._completed_key = dataset_name
            if not self._stopping:
                self.urban_ready.emit(
                    {
                        "outlines": outlines,
                        "source": source,
                    }
                )
        except Exception as exc:
            logger.warning("Urban outline update failed: %s", exc, exc_info=True)
            if not self._stopping:
                self.urban_failed.emit({"banner": f"Urban outline: {exc}"})
        finally:
            with self._lock:
                self._running = False
