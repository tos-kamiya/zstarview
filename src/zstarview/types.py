from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union

import astropy
import numpy as np


@dataclass
class LunarEclipseInfo:
    """Contains data about a lunar eclipse for rendering."""

    is_eclipse: bool = False
    eclipse_type: Optional[str] = None
    # Alt/Az of the center of the Earth's shadow
    shadow_center_alt: float = 0.0
    shadow_center_az: float = 0.0
    # Angular radii in degrees
    umbra_radius_deg: float = 0.0
    penumbra_radius_deg: float = 0.0
    moon_radius_deg: float = 0.0


@dataclass
class SolarEclipseInfo:
    """Contains data about a solar eclipse for rendering / dimming."""

    is_eclipse: bool = False
    eclipse_type: Optional[str] = None  # "partial" / "annular" / "total" / None
    sep_deg: float = 0.0  # Sun–Moon center separation [deg]
    obscuration: float = 0.0  # Fraction of the Sun's disk obscured by the Moon [0.0, 1.0]


@dataclass
class PlanetBody:
    name: str
    alt: float
    az: float
    symbol: str
    is_visible: bool
    phase_angle: Optional[float] = None  # moon only
    lunar_eclipse_info: Optional[LunarEclipseInfo] = None  # moon only
    solar_eclipse_info: Optional[SolarEclipseInfo] = None  # sun only


@dataclass
class ViewerData:
    """Contains information about the observer."""

    location: Tuple[float, float]  # (lat, lon)
    timezone_name: str
    city_name: str
    view_center: Tuple[float, float] = (90.0, 180.0)


@dataclass
class CelestialData:
    """Container for all calculated celestial data for a specific time."""

    time: astropy.time.Time
    planets: List[PlanetBody]
    stars: Dict[str, np.ndarray]
    celestial_equator_points: List[Tuple[float, float]]
    ecliptic_points: List[Tuple[float, float]]
    horizon_points: List[Tuple[float, float]]


@dataclass
class ScreenGeometry:
    """Screen geometry for drawing."""

    center: Tuple[int, int]
    radius: int


CelestialObject = Union[PlanetBody, Dict[str, Any]]


@dataclass(frozen=True, slots=True)
class HatchConfig:
    """雲ハッチの固定パラメータセット."""
    tile_px: int           # 縞のピッチ（px）
    line_px: int           # 縞の線幅（px）
    angle_deg: float       # 縞の傾き（deg）
    strength: int          # α強度（0..255）
    phase_px: float        # 位相シフト（px）
    edge_px: float         # 縞エッジのソフト化幅（px）

    def as_key(self) -> tuple:
        """キャッシュ用の厳密キー（浮動小数は丸める）"""
        return (
            int(self.tile_px),
            int(self.line_px),
            round(float(self.angle_deg), 6),
            int(max(0, min(255, self.strength))),
            round(float(self.phase_px) % max(1.0, float(self.tile_px)), 6),
            max(0.0, float(self.edge_px)),
        )
