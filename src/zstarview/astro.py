import io
import math
import sys
import threading
from contextlib import contextmanager
from pathlib import Path
from typing import Any, List, Optional, Tuple, TypedDict, cast
from urllib.request import Request, urlopen

import astropy
import astropy.time
import astropy.units as u
import numpy as np
import polars as pl
import skyfield.api
import skyfield.iokit
from astropy.coordinates import AltAz, EarthLocation, GeocentricTrueEcliptic, SkyCoord
from skyfield.api import Loader, Topos
from skyfield.magnitudelib import planetary_magnitude

from .asterisms import ASTERISM_REQUIRED_SOURCE_IDS
from .paths import (
    CACHE_PATH,
    EPHEMERIS_FILENAME,
    EPHEMERIS_URL,
    FIELD_OF_VIEW_DEG,
    PLANET_IDS,
    PLANET_SYMBOLS,
    STAR_FIELD_OF_VIEW_DEG,
)
from .user_agent import build_user_agent
from .types import (
    DeepSkyTable,
    LunarEclipseInfo,
    PlanetBody,
    SolarEclipseInfo,
    StarCatalogMeta,
    StarsTable,
)

# Skyfield ephemeris cache loader (separate from UI)
_cache_path = Path(CACHE_PATH)
_cache_path.mkdir(parents=True, exist_ok=True)


def _skyfield_urlopen(url: object, *args: object, **kwargs: object):
    if isinstance(url, str):
        url = Request(url, headers={"User-Agent": build_user_agent("skyfield-loader")})
    return urlopen(url, *args, **kwargs)


# Skyfield downloads ephemeris files through this module-level urlopen.
skyfield.iokit.urlopen = _skyfield_urlopen
_starfield_load = Loader(str(_cache_path))
_starfield_load.urls[EPHEMERIS_FILENAME] = EPHEMERIS_URL.rpartition("/")[0] + "/"
_ephemeris_lock = threading.Lock()
_cached_ephemeris: Any | None = None

_ICRS_UNIT_BASIS = SkyCoord(
    ra=np.array([0.0, 90.0, 0.0]) * u.deg,
    dec=np.array([0.0, 0.0, 90.0]) * u.deg,
    frame="icrs",
)


@contextmanager
def _temporary_standard_stream_fallback():
    """Provide dummy stdio streams for pythonw/gui-script downloads on Windows."""

    original_stdout = sys.stdout
    original_stderr = sys.stderr
    fallback_stdout = io.StringIO() if original_stdout is None else None
    fallback_stderr = io.StringIO() if original_stderr is None else None
    try:
        if fallback_stdout is not None:
            sys.stdout = fallback_stdout
        if fallback_stderr is not None:
            sys.stderr = fallback_stderr
        yield
    finally:
        sys.stdout = original_stdout
        sys.stderr = original_stderr


def load_ephemeris() -> Any:
    """Load the configured ephemeris, even under Windows gui-script launches."""
    global _cached_ephemeris
    if _cached_ephemeris is not None:
        return _cached_ephemeris
    with _ephemeris_lock:
        if _cached_ephemeris is None:
            with _temporary_standard_stream_fallback():
                _cached_ephemeris = _starfield_load(EPHEMERIS_FILENAME)
    return _cached_ephemeris


class StarCatalogArrays(TypedDict):
    """Pre-normalized star catalog arrays for fast repeated sky updates."""

    catalog_index: np.ndarray
    unit_vectors: np.ndarray
    vmag: np.ndarray
    bv: np.ndarray
    size_scale: np.ndarray
    color_base: np.ndarray
class DeepSkyCatalogArrays(TypedDict):
    """Pre-normalized deep-sky catalog arrays for fast repeated sky updates."""

    id: np.ndarray
    name: np.ndarray
    type: np.ndarray
    ra_h: np.ndarray
    dec: np.ndarray
    unit_vectors: np.ndarray
    vmag: np.ndarray
    major_arcmin: np.ndarray
    minor_arcmin: np.ndarray
    pa_deg: np.ndarray


def prepare_star_catalog_arrays(
    star_df: pl.DataFrame,
    *,
    max_vmag: float | None = None,
    vmag_brightness_scale: float = -0.39,
) -> StarCatalogArrays:
    """Normalize a Polars star catalog to NumPy arrays once at startup."""
    catalog_index = np.arange(star_df.height, dtype=np.int32)
    ra_rad = np.radians(star_df["RAh"].cast(pl.Float64, strict=False).to_numpy() * 15.0)
    dec_rad = np.radians(star_df["Dec"].cast(pl.Float64, strict=False).to_numpy())
    sin_dec = np.sin(dec_rad)
    cos_dec = np.cos(dec_rad)
    unit_vectors = np.column_stack(
        (
            cos_dec * np.cos(ra_rad),
            cos_dec * np.sin(ra_rad),
            sin_dec,
        )
    )
    vmag = star_df["Vmag"].cast(pl.Float64, strict=False).to_numpy()
    bv = star_df["B-V"].cast(pl.Float64, strict=False).fill_null(np.nan).to_numpy()
    v_ref = 1.0
    scale = float(vmag_brightness_scale)
    L_raw = 10.0 ** (scale * (vmag - v_ref))
    size_scale = np.power(np.clip(L_raw, 0.0, 1.0), 0.3)
    color_base = np.clip(np.power(L_raw, 0.6), 0.0, 1.0)

    if max_vmag is not None:
        mask = vmag <= float(max_vmag)
        catalog_index = catalog_index[mask]
        unit_vectors = unit_vectors[mask]
        vmag = vmag[mask]
        bv = bv[mask]
        size_scale = size_scale[mask]
        color_base = color_base[mask]

    return {
        "catalog_index": catalog_index,
        "unit_vectors": unit_vectors,
        "vmag": vmag,
        "bv": bv,
        "size_scale": size_scale,
        "color_base": color_base,
    }


def prepare_star_catalog_meta(star_df: pl.DataFrame) -> StarCatalogMeta:
    """Build sparse name/source-id lookup tables keyed by full catalog index."""
    names_raw = star_df["Name"].fill_null("").to_numpy() if "Name" in star_df.columns else np.array([], dtype=object)
    if names_raw.size > 0:
        names_str = np.array([str(value).strip() for value in names_raw], dtype=object)
        name_mask = np.array([bool(value) for value in names_str], dtype=bool)
        name_indices = np.nonzero(name_mask)[0].astype(np.int32, copy=False)
        names = names_str[name_mask]
    else:
        name_indices = np.array([], dtype=np.int32)
        names = np.array([], dtype=object)

    source_raw = (
        star_df["SourceId"].fill_null("").to_numpy()
        if "SourceId" in star_df.columns
        else np.array([], dtype=object)
    )
    if source_raw.size > 0:
        source_str = np.array([str(value).strip() for value in source_raw], dtype=object)
        source_mask = np.array([value in ASTERISM_REQUIRED_SOURCE_IDS for value in source_str], dtype=bool)
        source_id_indices = np.nonzero(source_mask)[0].astype(np.int32, copy=False)
        source_ids = source_str[source_mask]
    else:
        source_id_indices = np.array([], dtype=np.int32)
        source_ids = np.array([], dtype=object)

    return StarCatalogMeta(
        name_indices=name_indices,
        names=names,
        source_id_indices=source_id_indices,
        source_ids=source_ids,
    )


def _lookup_sparse_star_meta(
    catalog_indices: np.ndarray,
    meta_indices: np.ndarray,
    values: np.ndarray,
) -> np.ndarray:
    if catalog_indices.size == 0 or meta_indices.size == 0:
        return np.full(catalog_indices.shape, "", dtype=object)
    pos = np.searchsorted(meta_indices, catalog_indices)
    valid = pos < meta_indices.size
    out = np.full(catalog_indices.shape, "", dtype=object)
    if np.any(valid):
        matched = valid.copy()
        matched[valid] = meta_indices[pos[valid]] == catalog_indices[valid]
        out[matched] = values[pos[matched]]
    return out


def lookup_star_name(meta: StarCatalogMeta | None, catalog_index: int) -> str:
    if meta is None or meta.name_indices.size == 0:
        return ""
    pos = int(np.searchsorted(meta.name_indices, int(catalog_index)))
    if pos >= meta.name_indices.size or int(meta.name_indices[pos]) != int(catalog_index):
        return ""
    return str(meta.names[pos]).strip()


def lookup_star_source_id(meta: StarCatalogMeta | None, catalog_index: int) -> str:
    if meta is None or meta.source_id_indices.size == 0:
        return ""
    pos = int(np.searchsorted(meta.source_id_indices, int(catalog_index)))
    if pos >= meta.source_id_indices.size or int(meta.source_id_indices[pos]) != int(catalog_index):
        return ""
    return str(meta.source_ids[pos]).strip()


def resolve_star_names(stars: StarsTable, meta: StarCatalogMeta | None) -> np.ndarray:
    catalog_indices = np.asarray(stars["star_index"], dtype=np.int32)
    return _lookup_sparse_star_meta(catalog_indices, meta.name_indices, meta.names) if meta is not None else np.full(catalog_indices.shape, "", dtype=object)


def resolve_star_source_ids(stars: StarsTable, meta: StarCatalogMeta | None) -> np.ndarray:
    catalog_indices = np.asarray(stars["star_index"], dtype=np.int32)
    return (
        _lookup_sparse_star_meta(catalog_indices, meta.source_id_indices, meta.source_ids)
        if meta is not None
        else np.full(catalog_indices.shape, "", dtype=object)
    )


def prepare_deep_sky_catalog_arrays(dso_df: pl.DataFrame) -> DeepSkyCatalogArrays:
    """Normalize a DSO catalog to NumPy arrays once at startup."""
    ids = dso_df["Id"].to_numpy()
    names = dso_df["Name"].to_numpy()
    types = dso_df["Type"].to_numpy()
    ra_h = dso_df["RAh"].cast(pl.Float64, strict=False).to_numpy()
    dec = dso_df["Dec"].cast(pl.Float64, strict=False).to_numpy()
    ra_rad = np.radians(ra_h * 15.0)
    dec_rad = np.radians(dec)
    sin_dec = np.sin(dec_rad)
    cos_dec = np.cos(dec_rad)
    unit_vectors = np.column_stack(
        (
            cos_dec * np.cos(ra_rad),
            cos_dec * np.sin(ra_rad),
            sin_dec,
        )
    )
    vmag = dso_df["Vmag"].cast(pl.Float64, strict=False).fill_null(np.nan).to_numpy()
    major_arcmin = dso_df["MajorArcmin"].cast(pl.Float64, strict=False).fill_null(np.nan).to_numpy()
    minor_arcmin = dso_df["MinorArcmin"].cast(pl.Float64, strict=False).fill_null(np.nan).to_numpy()
    pa_deg = dso_df["PAdeg"].cast(pl.Float64, strict=False).fill_null(np.nan).to_numpy()
    return {
        "id": ids,
        "name": names,
        "type": types,
        "ra_h": ra_h,
        "dec": dec,
        "unit_vectors": unit_vectors,
        "vmag": vmag,
        "major_arcmin": major_arcmin,
        "minor_arcmin": minor_arcmin,
        "pa_deg": pa_deg,
    }


def altaz_to_normalized_xy(
    alt: float,
    az: float,
    view_center: Tuple[float, float],
    *,
    edge_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> Tuple[float, float]:
    """Convert alt/az to normalized screen coordinates relative to a view center.

    Returns (nx, ny) where 1.0 equals `edge_fov_deg` degrees of angular distance.
    """
    center_alt, center_az = view_center
    alt1 = math.radians(center_alt)
    az1 = math.radians(center_az)
    alt2 = math.radians(alt)
    az2 = math.radians(az)

    cos_theta = math.sin(alt1) * math.sin(alt2) + math.cos(alt1) * math.cos(alt2) * math.cos(az2 - az1)
    theta = math.acos(max(-1.0, min(1.0, cos_theta)))

    edge_fov_rad = math.radians(max(1.0e-6, float(edge_fov_deg)))
    r = theta / edge_fov_rad

    dx = math.cos(alt2) * math.sin(az2 - az1)
    dy = math.cos(alt1) * math.sin(alt2) - math.sin(alt1) * math.cos(alt2) * math.cos(az2 - az1)
    length = math.hypot(dx, dy)
    if length != 0:
        dx /= length
        dy /= length
    nx = r * dx
    ny = -r * dy
    return (nx, ny)


def is_in_fov(alt: float, az: float, view_center: Tuple[float, float], *, fov_deg: float = FIELD_OF_VIEW_DEG) -> bool:
    """Check if a target at (alt, az) is within the field of view relative to view_center."""
    center_alt, center_az = view_center
    alt1, az1 = math.radians(center_alt), math.radians(center_az)
    alt2, az2 = math.radians(alt), math.radians(az)
    cos_theta = math.sin(alt1) * math.sin(alt2) + math.cos(alt1) * math.cos(alt2) * math.cos(az2 - az1)
    theta = math.acos(min(1.0, max(-1.0, cos_theta)))
    return math.degrees(theta) <= fov_deg


def is_in_fov_vectorized(
    alt: np.ndarray,
    az: np.ndarray,
    view_center: Tuple[float, float],
    *,
    fov_deg: float = FIELD_OF_VIEW_DEG,
) -> np.ndarray:
    """Vectorized check if targets at (alt, az) are within the field of view."""
    center_alt, center_az = view_center
    alt1, az1 = np.radians(center_alt), np.radians(center_az)
    alt2, az2 = np.radians(alt), np.radians(az)

    cos_theta = np.sin(alt1) * np.sin(alt2) + np.cos(alt1) * np.cos(alt2) * np.cos(az2 - az1)

    # Clip to avoid domain errors with arccos
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)

    return np.degrees(theta) <= fov_deg


def build_icrs_to_altaz_matrix(time_obj: astropy.time.Time, location: EarthLocation) -> np.ndarray:
    """Return a 3x3 rotation matrix from ICRS to the AltAz frame."""
    altaz_frame = AltAz(obstime=time_obj, location=location)
    transformed = _ICRS_UNIT_BASIS.transform_to(altaz_frame)
    return transformed.cartesian.xyz.to_value(u.one)


def apply_icrs_to_altaz_matrix(unit_vectors: np.ndarray, matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Rotate ICRS unit vectors with the precomputed matrix and return alt/az."""
    altaz_cart = (matrix @ unit_vectors.T).T  # Shape: (N, 3)
    north = altaz_cart[:, 0]
    east = altaz_cart[:, 1]
    up = altaz_cart[:, 2]

    alt_rad = np.arcsin(np.clip(up, -1.0, 1.0))
    az_rad = np.arctan2(east, north)
    az_deg = np.degrees(np.mod(az_rad, math.tau))
    return np.degrees(alt_rad), az_deg


def calculate_visible_stars(
    star_source: pl.DataFrame | StarCatalogArrays,
    lat: float,
    lon: float,
    observer_height_m: float,
    time_obj: astropy.time.Time,
    view_center: Tuple[float, float],
    content_fov_deg: float = STAR_FIELD_OF_VIEW_DEG,
    max_vmag: float | None = None,
    subset_indices: np.ndarray | None = None,
    star_data_policy: str = "scenic_view_scoped",
) -> Tuple[StarsTable, EarthLocation]:
    """Compute visible stars and return them with the observer location."""
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=observer_height_m * u.m)

    source_is_df = isinstance(star_source, pl.DataFrame)
    if source_is_df:
        cat = prepare_star_catalog_arrays(cast(pl.DataFrame, star_source), max_vmag=max_vmag)
    else:
        cat = cast(StarCatalogArrays, star_source)

    # Get pre-normalized arrays
    catalog_index = cat["catalog_index"]
    unit_vectors = cat["unit_vectors"]
    vmag = cat["vmag"]
    bv = cat["bv"]
    size_scale = cat["size_scale"]
    color_base = cat["color_base"]

    if subset_indices is not None:
        idx = np.asarray(subset_indices, dtype=np.int32)
        if idx.size == 0:
            empty: StarsTable = {
                "star_index": catalog_index[:0],
                "alt": np.array([], dtype=float),
                "az": np.array([], dtype=float),
                "vmag": np.array([], dtype=float),
                "bv": np.array([], dtype=float),
                "size_factor": np.array([], dtype=float),
                "color_factor_base": np.array([], dtype=float),
            }
            return (empty, location)
        catalog_index = catalog_index[idx]
        unit_vectors = unit_vectors[idx]
        vmag = vmag[idx]
        bv = bv[idx]
        size_scale = size_scale[idx]
        color_base = color_base[idx]

    if max_vmag is not None and not source_is_df:
        vlim = float(max_vmag)
        mag_mask = vmag <= vlim
        if not np.any(mag_mask):
            empty: StarsTable = {
                "star_index": catalog_index[:0],
                "alt": np.array([], dtype=float),
                "az": np.array([], dtype=float),
                "vmag": np.array([], dtype=float),
                "bv": np.array([], dtype=float),
                "size_factor": np.array([], dtype=float),
                "color_factor_base": np.array([], dtype=float),
            }
            return (empty, location)
        catalog_index = catalog_index[mag_mask]
        unit_vectors = unit_vectors[mag_mask]
        vmag = vmag[mag_mask]
        bv = bv[mag_mask]
        size_scale = size_scale[mag_mask]
        color_base = color_base[mag_mask]

    matrix = build_icrs_to_altaz_matrix(time_obj, location)
    alt, az = apply_icrs_to_altaz_matrix(unit_vectors, matrix)
    star_data_policy = str(star_data_policy).strip().lower()
    if star_data_policy == "positional_static":
        in_view_mask = np.ones(alt.shape, dtype=bool)
    else:
        in_view_mask = is_in_fov_vectorized(
            alt,
            az,
            view_center,
            fov_deg=content_fov_deg,
        )

    # Filter the results using the boolean mask
    visible_stars: StarsTable = {
        "star_index": catalog_index[in_view_mask],
        "alt": alt[in_view_mask],
        "az": az[in_view_mask],
        "vmag": vmag[in_view_mask],
        "bv": bv[in_view_mask],
        "size_factor": size_scale[in_view_mask],
        "color_factor_base": color_base[in_view_mask],
    }

    return (visible_stars, location)


def calculate_visible_deep_sky_objects(
    dso_catalog: DeepSkyCatalogArrays,
    lat: float,
    lon: float,
    observer_height_m: float,
    time_obj: astropy.time.Time,
    view_center: Tuple[float, float],
    content_fov_deg: float = STAR_FIELD_OF_VIEW_DEG,
) -> DeepSkyTable:
    """Compute visible deep-sky objects and return vectorized rows."""
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=observer_height_m * u.m)
    matrix = build_icrs_to_altaz_matrix(time_obj, location)
    alt, az = apply_icrs_to_altaz_matrix(dso_catalog["unit_vectors"], matrix)
    in_view_mask = is_in_fov_vectorized(alt, az, view_center, fov_deg=content_fov_deg)
    return {
        "id": dso_catalog["id"][in_view_mask],
        "name": dso_catalog["name"][in_view_mask],
        "type": dso_catalog["type"][in_view_mask],
        "alt": alt[in_view_mask],
        "az": az[in_view_mask],
        "vmag": dso_catalog["vmag"][in_view_mask],
        "major_arcmin": dso_catalog["major_arcmin"][in_view_mask],
        "minor_arcmin": dso_catalog["minor_arcmin"][in_view_mask],
        "pa_deg": dso_catalog["pa_deg"][in_view_mask],
    }


def radec_to_altaz(
    ra_hours: float,
    dec_deg: float,
    lat: float,
    lon: float,
    observer_height_m: float,
    time_obj: astropy.time.Time,
) -> Tuple[float, float]:
    """Convert ICRS RA/Dec to topocentric Alt/Az for a given observer/time."""
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=observer_height_m * u.m)
    altaz_frame = AltAz(obstime=time_obj, location=location)
    coords = SkyCoord(ra=(float(ra_hours) * 15.0) * u.deg, dec=float(dec_deg) * u.deg, frame="icrs")
    altaz_coords = coords.transform_to(altaz_frame)
    return float(altaz_coords.alt.deg), float(altaz_coords.az.deg)


def calculate_moon_phase_angle(observer: Topos, t: skyfield.timelib.Time, planets: Any) -> float:
    """Calculate the phase angle of the Moon (Sun-Moon separation at observer)."""
    moon = planets["moon"]
    sun = planets["sun"]
    e_to_moon = observer.at(t).observe(moon).apparent()
    e_to_sun = observer.at(t).observe(sun).apparent()
    return e_to_moon.separation_from(e_to_sun).degrees


def calculate_lunar_eclipse_data(
    t: astropy.time.Time,
    observer,
    planets: Any,
) -> LunarEclipseInfo:
    """Calculate lunar eclipse parameters for a given time and observer.

    Returns an `EclipseInfo` object containing:
    - whether an eclipse occurs,
    - the eclipse type (penumbral, partial, total),
    - shadow center alt/az,
    - umbra and penumbra angular radii,
    - apparent moon radius.
    """
    earth = planets["earth"]
    sun = planets["sun"]
    moon = planets["moon"]

    # Sun and Moon positions from Earth's center in GCRS
    sun_pos = earth.at(t).observe(sun)
    moon_pos = earth.at(t).observe(moon)

    # Exclude if Sun-Moon separation is not close to 180 deg
    separation_deg = sun_pos.separation_from(moon_pos).degrees
    if abs(separation_deg - 180.0) > 3.0:
        return LunarEclipseInfo(is_eclipse=False)

    # Astronomical constants
    R_earth_km = 6371.0
    R_sun_km = 696340.0
    R_moon_km = 1737.4

    # Distances from Earth's center
    D_sun_km = sun_pos.distance().km
    D_moon_km = moon_pos.distance().km

    # Angular radii
    earth_ang_rad = math.asin(R_earth_km / D_moon_km)
    sun_ang_rad = math.asin(R_sun_km / D_sun_km)
    moon_ang_rad = math.asin(R_moon_km / D_moon_km)

    # Umbra and penumbra radii
    umbra_radius_rad = max(0.0, earth_ang_rad - sun_ang_rad)
    penumbra_radius_rad = earth_ang_rad + sun_ang_rad

    # Angular distance between shadow center and Moon
    d_rad = math.radians(180.0 - separation_deg)

    # Classify eclipse type
    if d_rad > (penumbra_radius_rad + moon_ang_rad):
        eclipse_type = "none"
    elif d_rad > (umbra_radius_rad + moon_ang_rad):
        eclipse_type = "penumbral"
    elif d_rad > abs(umbra_radius_rad - moon_ang_rad):
        eclipse_type = "partial"
    else:
        eclipse_type = "total"

    # Shadow center as seen from the observer (opposite Sun direction)
    sun_apparent = observer.at(t).observe(sun).apparent()
    s_alt, s_az, _ = sun_apparent.altaz()
    shadow_center_alt = -s_alt.degrees
    shadow_center_az = (s_az.degrees + 180.0) % 360.0

    return LunarEclipseInfo(
        is_eclipse=(eclipse_type != "none"),
        eclipse_type=eclipse_type,
        shadow_center_alt=shadow_center_alt,
        shadow_center_az=shadow_center_az,
        umbra_radius_deg=math.degrees(umbra_radius_rad),
        penumbra_radius_deg=math.degrees(penumbra_radius_rad),
        moon_radius_deg=math.degrees(moon_ang_rad),
    )


SUN_RADIUS_KM = 695700.0
MOON_RADIUS_KM = 1737.4


def _angular_radius_rad(radius_km: float, distance_km: float) -> float:
    """Apparent angular radius [rad] = asin(R / D) (more precise than small-angle approximation)."""
    x = radius_km / max(1e-9, distance_km)
    return math.asin(max(-1.0, min(1.0, x)))


def _circle_overlap_area_fraction(R: float, r: float, d: float) -> float:
    """
    Calculates the overlapping area of two circles (radii R, r; center distance d)
    divided by the area of the solar disk (pi R^2).
    R, r, and d are in the same units (here, angular radius in radians).
    """
    # Separated
    if d >= R + r:
        return 0.0
    # Contained (the smaller circle is completely inside the larger one)
    if d <= abs(R - r):
        return (min(R, r) ** 2) / (R**2)

    R2, r2 = R * R, r * r
    # Angles (segment angles)
    alpha = math.acos((d * d + R2 - r2) / (2.0 * d * R))
    beta = math.acos((d * d + r2 - R2) / (2.0 * d * r))
    # Lens area (Brahmagupta's formula-like square root term)
    lens = R2 * alpha + r2 * beta - 0.5 * math.sqrt(max(0.0, (-d + R + r) * (d + R - r) * (d - R + r) * (d + R + r)))
    return lens / (math.pi * R2)


def calculate_solar_eclipse_data(
    t: astropy.time.Time,
    observer,
    planets: Any,
) -> SolarEclipseInfo:
    sun = planets["sun"]
    moon = planets["moon"]

    # Observer's topocentric apparent position
    obs = observer.at(t)
    sun_app = obs.observe(sun).apparent()
    moon_app = obs.observe(moon).apparent()

    # Angular distance between centers
    sep = sun_app.separation_from(moon_app)
    sep_deg = sep.degrees
    sep_rad = sep.radians

    # Apparent angular radius (calculated from distance each time)
    D_sun_km = sun_app.distance().km
    D_moon_km = moon_app.distance().km
    sun_rad = _angular_radius_rad(SUN_RADIUS_KM, D_sun_km)  # [rad]
    moon_rad = _angular_radius_rad(MOON_RADIUS_KM, D_moon_km)  # [rad]

    # Obscuration (area fraction of the Sun obscured by the Moon)
    obscuration = _circle_overlap_area_fraction(sun_rad, moon_rad, sep_rad)

    # Determine eclipse type
    if sep_rad >= sun_rad + moon_rad or obscuration <= 0.0:
        eclipse_type = "none"
    elif sep_rad <= abs(sun_rad - moon_rad):
        # The smaller circle is completely contained: total if Moon > Sun, annular if Sun > Moon
        eclipse_type = "total" if moon_rad >= sun_rad else "annular"
    else:
        eclipse_type = "partial"

    return SolarEclipseInfo(
        is_eclipse=(eclipse_type != "none"),
        eclipse_type=eclipse_type,
        sep_deg=float(sep_deg),
        obscuration=obscuration,
    )


def eclipse_factor_from_info(info: Optional[SolarEclipseInfo]) -> float:
    """Return dimming factor in [0,1]. 1=no eclipse (no dim)."""
    if not info or not info.is_eclipse or info.obscuration <= 1e-3:
        return 1.0
    obsc = max(0.0, min(1.0, info.obscuration))

    # Perceptual model: gets dark suddenly after ~90%
    f0 = 0.92
    is_total = info.eclipse_type == "total"
    k = 10.0 if is_total else 7.5

    s = 1.0 / (1.0 + math.exp(k * (obsc - f0)))  # 1->0
    # Totality is darker
    min_l = 0.02 if is_total else 0.15
    return min_l + (1.0 - min_l) * s


def calculate_planets(
    lat: float,
    lon: float,
    observer_height_m: float,
    astropy_time: astropy.time.Time,
    view_center: Tuple[float, float],
    planets: Any,
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> List[PlanetBody]:
    """Calculate all planetary bodies (Sun, Moon, planets)."""
    ts = skyfield.api.load.timescale()
    t = ts.from_astropy(astropy_time)
    observer = planets["earth"] + Topos(
        latitude_degrees=lat,
        longitude_degrees=lon,
        elevation_m=observer_height_m,
    )

    bodies: List[PlanetBody] = []
    for name, symbol in PLANET_SYMBOLS.items():
        planet_id = PLANET_IDS[name]
        planet = planets[planet_id]
        astrometric = observer.at(t).observe(planet).apparent()
        alt, az, _ = astrometric.altaz()
        is_visible = is_in_fov(alt.degrees, az.degrees, view_center, fov_deg=content_fov_deg)
        vmag = None
        if name not in ("sun", "moon"):
            try:
                value = float(planetary_magnitude(astrometric))
                if np.isfinite(value):
                    vmag = value
            except Exception:
                vmag = None

        pa = lei = None
        if name == "moon":
            pa = calculate_moon_phase_angle(observer, t, planets)
            lei = calculate_lunar_eclipse_data(t, observer, planets)
        sei = None
        if name == "sun":
            sei = calculate_solar_eclipse_data(t, observer, planets)

        bodies.append(
            PlanetBody(
                name=name,
                alt=alt.degrees,
                az=az.degrees,
                symbol=symbol,
                is_visible=is_visible,
                vmag=vmag,
                phase_angle=pa,
                lunar_eclipse_info=lei,
                solar_eclipse_info=sei,
            )
        )
    return bodies


def calculate_horizon_points() -> List[Tuple[float, float]]:
    """Generate altitude/azimuth samples along the horizon."""
    points: List[Tuple[float, float]] = []
    alt = 0.0
    for az in range(0, 360 + 5, 5):
        points.append((alt, float(az)))
    return points


def calculate_celestial_equator_points(location: EarthLocation, time: astropy.time.Time) -> List[Tuple[float, float]]:
    """Generate altitude/azimuth samples along the celestial equator."""
    points: List[Tuple[float, float]] = []
    for ra_deg in range(0, 360 + 5, 5):
        coord = SkyCoord(ra=ra_deg * u.deg, dec=0 * u.deg, frame="icrs")
        altaz = coord.transform_to(AltAz(obstime=time, location=location))
        points.append((float(altaz.alt.deg), float(altaz.az.deg)))
    return points


def calculate_ecliptic_points(location: EarthLocation, time: astropy.time.Time) -> List[Tuple[float, float]]:
    """Generate altitude/azimuth samples along the ecliptic."""
    points: List[Tuple[float, float]] = []
    for lon_deg in range(0, 360 + 5, 5):
        ecl = SkyCoord(lon=lon_deg * u.deg, lat=0 * u.deg, frame=GeocentricTrueEcliptic(obstime=time))
        icrs = ecl.transform_to("icrs")
        altaz = icrs.transform_to(AltAz(obstime=time, location=location))
        points.append((float(altaz.alt.deg), float(altaz.az.deg)))
    return points


def calculate_moon_render_data(
    sun_altaz: Optional[Tuple[float, float]],
    moon_altaz: Optional[Tuple[float, float]],
    view_center: Tuple[float, float],
    *,
    edge_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> Tuple[np.ndarray, float]:
    """
    Calculates all necessary data for rendering the moon.

    Returns a tuple of:
    - sun_dir_in_moon_frame: 3D vector of the sun's direction in the moon's local frame.
    - screen_rotation_deg: The angle to rotate the moon image on the screen.
    """

    # Helper to convert Alt/Az to a 3D vector in the observer's reference frame.
    # Y-up, X-East, Z-North (Right-handed system)
    def altaz_to_cartesian(alt: float, az: float) -> np.ndarray:
        alt_rad = math.radians(alt)
        az_rad = math.radians(az)  # Azimuth is measured from North (0) towards East (90)
        x = math.cos(alt_rad) * math.sin(az_rad)
        y = math.sin(alt_rad)
        z = math.cos(alt_rad) * math.cos(az_rad)
        return np.array([x, y, z])

    sun_dir_in_observer_frame = None
    if sun_altaz:
        sun_dir_in_observer_frame = altaz_to_cartesian(sun_altaz[0], sun_altaz[1])

    sun_dir_in_moon_frame = np.array([0.0, 0.0, 0.0])
    if sun_dir_in_observer_frame is not None and moon_altaz is not None:
        m_alt, m_az = moon_altaz

        # If the moon is very close to zenith or nadir (alt ~ +/-90 deg),
        # the calculation for the local coordinate frame can become unstable.
        # To avoid this singularity, we perturb the altitude slightly.
        if abs(m_alt) > 89.9999:
            m_alt = math.copysign(89.9999, m_alt)

        moon_vec = altaz_to_cartesian(m_alt, m_az)
        z_axis = -moon_vec / np.linalg.norm(moon_vec)
        zenith_vec = np.array([0.0, 1.0, 0.0])
        y_axis = zenith_vec - np.dot(zenith_vec, z_axis) * z_axis
        y_axis /= np.linalg.norm(y_axis)
        x_axis = np.cross(y_axis, z_axis)

        sun_dir_in_moon_frame = np.array(
            [
                np.dot(sun_dir_in_observer_frame, x_axis),
                np.dot(sun_dir_in_observer_frame, y_axis),
                np.dot(sun_dir_in_observer_frame, z_axis),
            ]
        )

    # Calculate the screen rotation for the moon image
    screen_rotation_deg = 0.0
    if moon_altaz:
        # Get screen coordinates for the moon's center and a point just "above" it
        m_alt, m_az = moon_altaz
        nx_center, ny_center = altaz_to_normalized_xy(
            m_alt,
            m_az,
            view_center,
            edge_fov_deg=edge_fov_deg,
        )
        nx_up, ny_up = altaz_to_normalized_xy(
            m_alt + 0.1,
            m_az,
            view_center,
            edge_fov_deg=edge_fov_deg,
        )

        # Calculate the angle of the "up" direction on the screen
        dx = nx_up - nx_center
        dy = ny_up - ny_center
        if abs(dx) > 1e-6 or abs(dy) > 1e-6:
            screen_up_angle_rad = math.atan2(dy, dx)
            # The image's "up" is vertical (points to -Y), so rotate to match the screen's "up"
            screen_rotation_deg = math.degrees(screen_up_angle_rad) - 90

    return sun_dir_in_moon_frame, screen_rotation_deg
