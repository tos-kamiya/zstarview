from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union, TypeAlias, TypedDict

import astropy
import astropy.time
import numpy as np

from .paths import FIELD_OF_VIEW_DEG

# --- Unit-oriented type aliases (contract clarity) ---
Deg: TypeAlias = float
AltDeg: TypeAlias = Deg
AzDeg: TypeAlias = Deg
LatDeg: TypeAlias = Deg
LonDeg: TypeAlias = Deg
LocationLatLon: TypeAlias = Tuple[LatDeg, LonDeg]
ViewCenterAltAz: TypeAlias = Tuple[AltDeg, AzDeg]


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
    alt: AltDeg  # altitude in degrees
    az: AzDeg  # azimuth in degrees (0=N, 90=E)
    symbol: str
    is_visible: bool
    vmag: Optional[float] = None  # planets (except sun/moon) when available
    phase_angle: Optional[float] = None  # moon only
    lunar_eclipse_info: Optional[LunarEclipseInfo] = None  # moon only
    solar_eclipse_info: Optional[SolarEclipseInfo] = None  # sun only

    @property
    def alt_deg(self) -> AltDeg:
        """Altitude in degrees (alias for `alt`)."""
        return self.alt

    @alt_deg.setter
    def alt_deg(self, value: AltDeg) -> None:
        self.alt = value

    @property
    def az_deg(self) -> AzDeg:
        """Azimuth in degrees (alias for `az`, 0=N, 90=E)."""
        return self.az

    @az_deg.setter
    def az_deg(self, value: AzDeg) -> None:
        self.az = value


@dataclass
class ViewerData:
    """Contains information about the observer."""

    location: LocationLatLon  # (lat_deg, lon_deg)
    timezone_name: str
    city_name: str
    view_center: ViewCenterAltAz = (90.0, 180.0)  # (alt_deg, az_deg)
    content_fov_deg: float = float(FIELD_OF_VIEW_DEG)
    observer_height_m: float = 1.7
    location_height_label: str | None = None
    location_height_m: float | None = None
    show_observer_height: bool = False

    @property
    def lat_deg(self) -> LatDeg:
        return self.location[0]

    @property
    def lon_deg(self) -> LonDeg:
        return self.location[1]

    @property
    def view_alt_deg(self) -> AltDeg:
        return self.view_center[0]

    @property
    def view_az_deg(self) -> AzDeg:
        return self.view_center[1]


@dataclass
class CelestialData:
    """Container for all calculated celestial data for a specific time."""

    time: astropy.time.Time
    planets: List[PlanetBody]
    stars: "StarsTable"
    deep_sky_objects: "DeepSkyTable"
    celestial_equator_points: List[Tuple[float, float]]  # (alt_deg, az_deg)
    ecliptic_points: List[Tuple[float, float]]  # (alt_deg, az_deg)
    horizon_points: List[Tuple[float, float]]  # (alt_deg, az_deg)


@dataclass
class ScreenGeometry:
    """Screen geometry for drawing."""

    center: Tuple[int, int]
    radius: int


@dataclass(frozen=True)
class UrbanOutlinePolyline:
    """Renderable urban outline polyline with source building height."""

    points: List[Tuple[float, float]]
    height_m: float
    source: str = "base"


CelestialObject = Union[PlanetBody, Dict[str, Any]]


class StarsTable(TypedDict):
    """Vectorized star table contract used across astro -> render.

    All arrays must be aligned by index and share the same length.
    - name: star display names
    - alt: altitude [deg]
    - az: azimuth [deg] (0=N, 90=E)
    - vmag: visual magnitude
    - bv: B-V color index
    """

    name: np.ndarray
    source_id: np.ndarray
    alt: np.ndarray
    az: np.ndarray
    vmag: np.ndarray
    bv: np.ndarray
    size_factor: np.ndarray
    color_factor_base: np.ndarray


class DeepSkyTable(TypedDict):
    """Vectorized deep-sky table contract used across astro -> render."""

    id: np.ndarray
    name: np.ndarray
    type: np.ndarray
    alt: np.ndarray
    az: np.ndarray
    vmag: np.ndarray
    major_arcmin: np.ndarray
    minor_arcmin: np.ndarray
    pa_deg: np.ndarray
