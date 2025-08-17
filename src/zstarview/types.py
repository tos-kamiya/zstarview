from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import astropy
import numpy as np


@dataclass
class EclipseInfo:
    """Contains data about a lunar eclipse for rendering."""

    is_eclipse: bool = False
    # Alt/Az of the center of the Earth's shadow
    shadow_center_alt: float = 0.0
    shadow_center_az: float = 0.0
    # Angular radii in degrees
    umbra_radius_deg: float = 0.0
    penumbra_radius_deg: float = 0.0
    moon_radius_deg: float = 0.0


@dataclass
class PlanetBody:
    name: str
    alt: float
    az: float
    symbol: str
    is_visible: bool
    phase_angle: Optional[float] = None  # moon only
    eclipse_info: Optional[EclipseInfo] = None  # moon only


@dataclass
class ViewerData:
    """Contains information about the observer."""

    location: Tuple[float, float]  # (lat, lon)
    timezone_name: str
    city_name: str
    view_center: Tuple[float, float] = (90.0, 180.0)


@dataclass
class SkyData:
    """Container for all calculated sky data for a specific time."""

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