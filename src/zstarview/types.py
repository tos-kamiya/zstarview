from dataclasses import dataclass
from typing import List, Tuple, Optional


@dataclass
class StarData:
    name: str
    alt: float
    az: float
    vmag: float
    bv: float


@dataclass
class PlanetBody:
    name: str
    alt: float
    az: float
    symbol: str
    phase_angle: Optional[float] = None  # moon only


@dataclass
class SkyData:
    """Container for all calculated sky data for a specific time and location."""

    location: Tuple[float, float]  # (lat, lon)
    time: object
    planets: List[PlanetBody]
    stars: List[StarData]
    celestial_equator_points: List[Tuple[float, float]]
    ecliptic_points: List[Tuple[float, float]]
    horizon_points: List[Tuple[float, float]]
    timezone_name: str
    city_name: str
    view_center: Tuple[float, float] = (90.0, 180.0)
