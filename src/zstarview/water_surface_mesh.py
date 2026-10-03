from __future__ import annotations

from pyproj import CRS, Transformer

def make_local_transformer(lat_deg: float, lon_deg: float) -> Transformer:
    local_crs = CRS.from_proj4(
        f"+proj=aeqd +lat_0={lat_deg} +lon_0={lon_deg} +datum=WGS84 +units=m +no_defs"
    )
    return Transformer.from_crs("EPSG:4326", local_crs, always_xy=True)


def project_ring_xy(
    ring_lonlat: tuple[tuple[float, float], ...],
    transformer: Transformer,
) -> list[tuple[float, float]]:
    lon = [point[0] for point in ring_lonlat]
    lat = [point[1] for point in ring_lonlat]
    x, y = transformer.transform(lon, lat)
    return [(float(px), float(py)) for px, py in zip(x, y)]
