from dataclasses import dataclass
from typing import Any, TypeAlias, TypedDict, Union

import astropy
import astropy.time
import numpy as np

# --- Unit-oriented type aliases (contract clarity) ---
Deg: TypeAlias = float
AltDeg: TypeAlias = Deg
AzDeg: TypeAlias = Deg
LatDeg: TypeAlias = Deg
LonDeg: TypeAlias = Deg
LocationLatLon: TypeAlias = tuple[LatDeg, LonDeg]
ViewCenterAltAz: TypeAlias = tuple[AltDeg, AzDeg]


def _normalize_edge_and_content_fov(edge_fov_deg: float, content_fov_deg: float) -> tuple[float, float]:
    """Return a FOV pair that satisfies content >= edge."""
    edge = float(edge_fov_deg)
    content = float(content_fov_deg)
    content = max(content, edge)
    return edge, content


@dataclass
class LunarEclipseInfo:
    """Contains data about a lunar eclipse for rendering."""

    is_eclipse: bool = False
    eclipse_type: str | None = None
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
    eclipse_type: str | None = None  # "partial" / "annular" / "total" / None
    sep_deg: float = 0.0  # Sun-Moon center separation [deg]
    obscuration: float = 0.0  # Fraction of the Sun's disk obscured by the Moon [0.0, 1.0]


@dataclass
class PlanetBody:
    name: str
    alt: AltDeg  # altitude in degrees
    az: AzDeg  # azimuth in degrees (0=N, 90=E)
    symbol: str
    is_visible: bool
    vmag: float | None = None  # planets (except sun/moon) when available
    phase_angle: float | None = None  # moon only
    lunar_eclipse_info: LunarEclipseInfo | None = None  # moon only
    solar_eclipse_info: SolarEclipseInfo | None = None  # sun only

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


@dataclass(frozen=True)
class ViewerData:
    """Contains information about the observer."""

    location: LocationLatLon  # (lat_deg, lon_deg)
    timezone_name: str
    city_name: str
    view_center: ViewCenterAltAz = (90.0, 180.0)  # (alt_deg, az_deg)
    edge_fov_deg: float = 90.0
    content_fov_deg: float = 115.0
    observer_height_m: float = 1.7
    height_add_m: float = 1.7
    ground_elevation_m: float = 0.0
    location_height_label: str | None = None
    location_height_m: float = 0.0

    def __post_init__(self) -> None:
        edge_fov_deg, content_fov_deg = _normalize_edge_and_content_fov(
            self.edge_fov_deg,
            self.content_fov_deg,
        )
        object.__setattr__(self, "edge_fov_deg", edge_fov_deg)
        object.__setattr__(self, "content_fov_deg", content_fov_deg)

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


@dataclass(frozen=True)
class ViewProjection:
    """Projection parameters for mapping alt/az into the render disc."""

    view_center: ViewCenterAltAz = (90.0, 180.0)  # (alt_deg, az_deg)
    edge_fov_deg: float = 90.0
    content_fov_deg: float = 115.0

    def __post_init__(self) -> None:
        edge_fov_deg, content_fov_deg = _normalize_edge_and_content_fov(
            self.edge_fov_deg,
            self.content_fov_deg,
        )
        object.__setattr__(self, "edge_fov_deg", edge_fov_deg)
        object.__setattr__(self, "content_fov_deg", content_fov_deg)


@dataclass
class CelestialData:
    """Container for all calculated celestial data for a specific time."""

    time: astropy.time.Time
    planets: list[PlanetBody]
    stars: "StarsTable"
    deep_sky_objects: "DeepSkyTable"
    celestial_equator_points: list[tuple[float, float]]  # (alt_deg, az_deg)
    ecliptic_points: list[tuple[float, float]]  # (alt_deg, az_deg)
    horizon_points: list[tuple[float, float]]  # (alt_deg, az_deg)
    star_catalog_meta: "StarCatalogMeta | None" = None
    star_time: astropy.time.Time | None = None


@dataclass
class ScreenGeometry:
    """Screen geometry for drawing."""

    center: tuple[int, int]
    radius: int


@dataclass(frozen=True)
class UrbanOutlinePolyline:
    """Renderable urban outline polyline with source building height."""

    points: list[tuple[float, float]]
    height_m: float
    distance_km: float = float("inf")
    source: str = "base"


CelestialObject = Union[PlanetBody, dict[str, Any]]


@dataclass(frozen=True)
class StarCatalogMeta:
    """Sparse metadata lookup tables keyed by the full star catalog index."""

    name_indices: np.ndarray
    names: np.ndarray
    source_id_indices: np.ndarray
    source_ids: np.ndarray


class StarsTable(TypedDict):
    """Vectorized star table contract used across astro -> render.

    All arrays must be aligned by index and share the same length.
    - star_index: row index in the full prepared star catalog
    - alt: altitude [deg]
    - az: azimuth [deg] (0=N, 90=E)
    - vmag: visual magnitude
    - bv: B-V color index
    """

    star_index: np.ndarray
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
