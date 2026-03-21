from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from pyproj import CRS, Transformer


@dataclass(slots=True)
class GeoArea:
    """Minimal geostationary area definition used by cloud-disc sampling."""

    area_id: str
    description: str
    proj_id: str
    projection: CRS
    width: int
    height: int
    area_extent: tuple[float, float, float, float]
    _lonlat_to_proj: Transformer = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self._lonlat_to_proj = Transformer.from_crs("EPSG:4326", self.projection, always_xy=True)

    @property
    def crs(self) -> CRS:
        return self.projection

    @property
    def proj_dict(self) -> dict[str, Any]:
        return dict(self.projection.to_dict())

    @property
    def shape(self) -> tuple[int, int]:
        return (self.height, self.width)

    def get_xy_from_lonlat(self, lon: Any, lat: Any) -> tuple[Any, Any]:
        projected_x, projected_y = self._lonlat_to_proj.transform(lon, lat)
        min_x, min_y, max_x, max_y = self.area_extent
        pixel_size_x = (max_x - min_x) / float(self.width)
        pixel_size_y = (max_y - min_y) / float(self.height)
        x = (np.asarray(projected_x, dtype=np.float64) - min_x) / pixel_size_x
        y = (max_y - np.asarray(projected_y, dtype=np.float64)) / pixel_size_y
        if x.ndim == 0 and y.ndim == 0:
            return float(x), float(y)
        return x, y
