# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from PySide6.QtCore import QObject, Signal

from ..data.skyscraper_tiles import (
    SKYSCRAPER_OUTER_RADIUS_KM,
    SkyscraperSeedTile,
    select_skyscraper_seed_tiles_for_viewer,
    skyscraper_tile_derived_dir,
)
from ..data.import_overture_buildings import (
    DEFAULT_FETCH_RADIUS_KM,
    OVERTURE_CACHE_TTL_DAYS,
    derive_dataset_name,
    is_derived_dataset_stale,
    import_overture_buildings_for_bbox,
    import_overture_buildings,
    resolve_overture_release_for_cache_root,
)
from ..urban_outline_layer import resolve_urban_outline_layer_for_viewer
from ..paths import OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR, SKYSCRAPER_TILES_FILE
from ..paths import CACHE_PATH
from ..types import UrbanOutlinePolyline, ViewerData

logger = logging.getLogger(__name__)

SKYSCRAPER_MIN_HEIGHT_M = 150.0


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
        skyscraper_outer_radius_km: float = SKYSCRAPER_OUTER_RADIUS_KM,
        feature_type: str = "both",
        overturemaps_bin: str = "overturemaps",
        skyscraper_only: bool = False,
        skyscraper_seed_file: Path | None = None,
        skyscraper_derived_root_dir: Path | None = None,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._derived_root_dir = Path(derived_root_dir)
        self._min_building_height_m = float(min_building_height_m)
        self._radius_km = float(radius_km)
        self._skyscraper_outer_radius_km = float(skyscraper_outer_radius_km)
        self._feature_type = str(feature_type)
        self._overturemaps_bin = str(overturemaps_bin)
        self._skyscraper_only = bool(skyscraper_only)
        self._skyscraper_seed_file = Path(skyscraper_seed_file or SKYSCRAPER_TILES_FILE)
        self._skyscraper_derived_root_dir = Path(skyscraper_derived_root_dir or OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR)
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
        required_dirs = () if self._skyscraper_only else self._required_derived_dirs(viewer_data)
        skyscraper_tiles = self._selected_skyscraper_tiles(viewer_data)
        skyscraper_dirs = tuple(
            skyscraper_tile_derived_dir(
                tile,
                derived_root_dir=self._skyscraper_derived_root_dir,
            )
            for tile in skyscraper_tiles
        )
        with self._lock:
            if self._stopping or self._running:
                return False
            if self._completed_key == dataset_name:
                return False
            self._running = True

        now = datetime.now(timezone.utc)
        if all(
            path.exists() and not is_derived_dataset_stale(path, ttl_days=OVERTURE_CACHE_TTL_DAYS, now_utc=now)
            for _, path in required_dirs
        ) and all(
            path.exists() and not is_derived_dataset_stale(path, ttl_days=OVERTURE_CACHE_TTL_DAYS, now_utc=now)
            for path in skyscraper_dirs
        ):
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
            now = datetime.now(timezone.utc)
            current_overture_release = resolve_overture_release_for_cache_root(
                cache_root_dir=Path(CACHE_PATH),
                now_utc=now,
            )
            source = "cache"
            required_dirs = () if self._skyscraper_only else self._required_derived_dirs(viewer_data)
            for overture_feature_type, derived_dir in required_dirs:
                derived_exists = derived_dir.exists()
                derived_is_stale = derived_exists and is_derived_dataset_stale(
                    derived_dir,
                    ttl_days=OVERTURE_CACHE_TTL_DAYS,
                    now_utc=now,
                    expected_overture_release=current_overture_release,
                )
                if derived_exists and not derived_is_stale:
                    continue
                if derived_is_stale:
                    logger.info(
                        "Refreshing expired urban-outline cache from Overture: %s (ttl=%s days)",
                        derived_dir,
                        int(OVERTURE_CACHE_TTL_DAYS),
                    )
                elif reason == "initial":
                    logger.info("Downloading initial urban-outline data from Overture...")
                else:
                    logger.info("Downloading urban-outline data from Overture...")
                try:
                    import_overture_buildings(
                        lat_deg=float(viewer_data.lat_deg),
                        lon_deg=float(viewer_data.lon_deg),
                        radius_km=self._radius_km,
                        derived_root_dir=self._derived_root_dir,
                        min_building_height_m=self._min_building_height_m,
                        feature_type=overture_feature_type,
                        fmt="geojsonseq",
                        overturemaps_bin=self._overturemaps_bin,
                        dataset_name=derived_dir.parent.name,
                        keep_download=None,
                        no_stac=False,
                        overture_release=current_overture_release,
                        skip_release_lookup=True,
                        now_utc=now,
                    )
                    source = "overture"
                except Exception:
                    if derived_dir.exists():
                        logger.warning(
                            "Urban outline refresh failed; using stale cache: %s",
                            derived_dir,
                            exc_info=True,
                        )
                        source = "cache-stale"
                        continue
                    raise

            outlines = None
            if required_dirs:
                outlines = resolve_urban_outline_layer_for_viewer(
                    viewer_data,
                    derived_root_dir=self._derived_root_dir,
                    derived_dirs=tuple(path for _, path in required_dirs),
                )
            merged_outlines = outlines
            try:
                skyscraper_tiles = self._selected_skyscraper_tiles(viewer_data)
                skyscraper_dirs: list[Path] = []
                for tile in skyscraper_tiles:
                    derived_dir = skyscraper_tile_derived_dir(
                        tile,
                        derived_root_dir=self._skyscraper_derived_root_dir,
                    )
                    skyscraper_dirs.append(derived_dir)
                    derived_exists = derived_dir.exists()
                    derived_is_stale = derived_exists and is_derived_dataset_stale(
                        derived_dir,
                        ttl_days=OVERTURE_CACHE_TTL_DAYS,
                        now_utc=now,
                        expected_overture_release=current_overture_release,
                    )
                    if derived_exists and not derived_is_stale:
                        continue
                    if derived_is_stale:
                        logger.info(
                            "Refreshing expired skyscraper urban-outline cache z=%s x=%s y=%s (ttl=%s days)",
                            tile.zoom,
                            tile.x,
                            tile.y,
                            int(OVERTURE_CACHE_TTL_DAYS),
                        )
                    else:
                        logger.info(
                            "Downloading skyscraper urban-outline tile z=%s x=%s y=%s...",
                            tile.zoom,
                            tile.x,
                            tile.y,
                        )
                    try:
                        import_overture_buildings_for_bbox(
                            bbox=(
                                tile.envelope.min_lon_deg,
                                tile.envelope.min_lat_deg,
                                tile.envelope.max_lon_deg,
                                tile.envelope.max_lat_deg,
                            ),
                            derived_root_dir=self._skyscraper_derived_root_dir,
                            min_building_height_m=SKYSCRAPER_MIN_HEIGHT_M,
                            feature_type="building",
                            fmt="geojsonseq",
                            overturemaps_bin=self._overturemaps_bin,
                            dataset_name=tile.cache_key,
                            keep_download=None,
                            no_stac=False,
                            overture_release=current_overture_release,
                            skip_release_lookup=True,
                            now_utc=now,
                        )
                        source = "overture"
                    except Exception:
                        if derived_dir.exists():
                            logger.warning(
                                "Skyscraper urban-outline refresh failed; using stale cache: %s",
                                derived_dir,
                                exc_info=True,
                            )
                            source = "cache-stale"
                            continue
                        raise
                skyscraper_outlines = resolve_urban_outline_layer_for_viewer(
                    viewer_data,
                    derived_root_dir=self._skyscraper_derived_root_dir,
                    derived_dirs=tuple(skyscraper_dirs),
                    radius_km=self._skyscraper_outer_radius_km,
                    min_distance_km=self._radius_km,
                    min_height_m=max(SKYSCRAPER_MIN_HEIGHT_M, self._min_building_height_m),
                )
                merged_outlines = self._merge_outline_layers(outlines, skyscraper_outlines)
            except Exception as exc:
                logger.warning("Skyscraper urban-outline update failed: %s", exc, exc_info=True)
            with self._lock:
                if not self._stopping:
                    self._completed_key = dataset_name
            if not self._stopping:
                self.urban_ready.emit(
                    {
                        "outlines": merged_outlines,
                        "source": source,
                    }
                )
        except Exception as exc:
            missing_overturemaps = isinstance(exc, FileNotFoundError) and "overturemaps" in str(exc)
            if missing_overturemaps:
                logger.info("Urban outline unavailable: overturemaps CLI not found")
            else:
                logger.warning("Urban outline update failed: %s", exc, exc_info=True)
            if not self._stopping:
                if missing_overturemaps:
                    self.urban_failed.emit({"banner": "Urban outline: Could not find overturemaps"})
                else:
                    self.urban_failed.emit({"banner": f"Urban outline: {exc}"})
        finally:
            with self._lock:
                self._running = False

    def _required_feature_types(self) -> tuple[str, ...]:
        if self._feature_type == "both":
            return ("building", "building_part")
        return (self._feature_type,)

    def _required_derived_dirs(self, viewer_data: ViewerData) -> tuple[tuple[str, Path], ...]:
        return tuple(
            (
                overture_feature_type,
                self._derived_root_dir
                / derive_dataset_name(
                    float(viewer_data.lat_deg),
                    float(viewer_data.lon_deg),
                    self._radius_km,
                    overture_feature_type,
                    self._min_building_height_m,
                )
                / "bldg",
            )
            for overture_feature_type in self._required_feature_types()
        )

    def _selected_skyscraper_tiles(self, viewer_data: ViewerData) -> tuple[SkyscraperSeedTile, ...]:
        if self._skyscraper_outer_radius_km <= 0.0:
            return ()
        try:
            return select_skyscraper_seed_tiles_for_viewer(
                observer_lat_deg=float(viewer_data.lat_deg),
                observer_lon_deg=float(viewer_data.lon_deg),
                inner_radius_km=self._radius_km,
                outer_radius_km=self._skyscraper_outer_radius_km,
                seed_file=self._skyscraper_seed_file,
            )
        except Exception as exc:
            logger.warning("Could not load skyscraper seed tiles: %s", exc, exc_info=True)
            return ()

    @staticmethod
    def _merge_outline_layers(
        base_outlines,
        extra_outlines,
    ):
        merged = []
        if base_outlines:
            merged.extend(base_outlines)
        if extra_outlines:
            merged.extend(
                UrbanOutlinePolyline(
                    points=list(outline.points),
                    height_m=float(outline.height_m),
                    source="skyscraper",
                )
                for outline in extra_outlines
            )
        return merged or None
