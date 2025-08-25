import numpy as np
import math

EARTH_R_KM = 6371.0

def deg2rad(d): return d * math.pi / 180.0
def rad2deg(r): return r * 180.0 / math.pi

def enu_basis(lat_deg, lon_deg):
    lat, lon = deg2rad(lat_deg), deg2rad(lon_deg)
    cl, sl = math.cos(lat), math.sin(lat)
    co, so = math.cos(lon), math.sin(lon)
    u = np.array([cl*co, cl*so, sl])       # Up
    e = np.array([-so,   co,   0.0])       # East
    n = np.array([-sl*co,-sl*so,cl])       # North
    return e, n, u

def geodetic_to_ecef(lat_deg, lon_deg, R_km=EARTH_R_KM):
    lat, lon = deg2rad(lat_deg), deg2rad(lon_deg)
    cl, sl = math.cos(lat), math.sin(lat)
    co, so = math.cos(lon), math.sin(lon)
    return np.array([R_km*cl*co, R_km*cl*so, R_km*sl], dtype=float)

def azalt_to_dir_ecef(az_deg, alt_deg, lat0, lon0):
    e,n,u = enu_basis(lat0, lon0)
    az, alt = deg2rad(az_deg), deg2rad(alt_deg)
    ca, sa = math.cos(alt), math.sin(alt)
    saz, caz = math.sin(az), math.cos(az)
    d = ca*saz*e + ca*caz*n + sa*u
    return d / (np.linalg.norm(d) or 1.0)

def line_sphere_intersection(O, d, radius_km):
    # Solve |O + t d| = R
    b = 2.0*np.dot(O, d)
    c = np.dot(O,O) - radius_km*radius_km
    disc = b*b - 4.0*c
    if disc < 0: return None
    s = math.sqrt(disc)
    t1 = (-b - s)/2.0
    t2 = (-b + s)/2.0
    cands = [t for t in (t1,t2) if t>1e-6]
    return min(cands) if cands else None

def az_project_lonlat_grid(
        lat0, lon0, alt0, az0, radius_px, cloud_shell_km,
        alt_min_deg=-2.0, fov_deg: float = 90.0
    ):
    S = 2*radius_px
    xs = (np.arange(S)+0.5) - radius_px
    ys = radius_px - (np.arange(S)+0.5)
    X, Y = np.meshgrid(xs, ys)
    r = np.hypot(X, Y)

    # 画素半径→角距離（度）：Rpx が fov_deg に対応
    rho_deg = (90.0 * r / radius_px).astype(np.float32)

    # 視線ベクトル計算（従来通り）
    v0 = azalt_to_dir_ecef(az0, alt0, lat0, lon0)
    e_obs, n_obs, u_obs = enu_basis(lat0, lon0)
    u_obs = np.array(u_obs)
    O = geodetic_to_ecef(lat0, lon0)

    psi = np.arctan2(Y, X)
    def proj_tan(a):
        t = a - np.dot(a, v0)*v0
        n = np.linalg.norm(t) or 1.0
        return t/n
    tY = proj_tan(u_obs)
    tX = np.cross(v0, tY); tX /= (np.linalg.norm(tX) or 1.0)

    cd = np.cos(np.deg2rad(rho_deg)); sd = np.sin(np.deg2rad(rho_deg))
    b = (np.cos(psi)[...,None]*tX + np.sin(psi)[...,None]*tY)
    d = cd[...,None]*v0 + sd[...,None]*b
    d = d / (np.linalg.norm(d, axis=2, keepdims=True)+1e-12)

    # 地平線下マスク
    alt_rad = np.arcsin(np.dot(d, u_obs))
    visible = (np.degrees(alt_rad) >= alt_min_deg)
    mask_inside = (rho_deg <= fov_deg + 1e-6) & visible

    # 雲殻との交点（従来通り）
    O3 = O[None,None,:]
    H, W = S, S
    lon = np.full((H,W), np.nan, dtype=np.float32)
    lat = np.full((H,W), np.nan, dtype=np.float32)

    R = cloud_shell_km
    for y in range(H):
        dy = d[y]
        oy = O3[0,0]
        bq = 2.0*np.sum(oy*dy, axis=1)
        c = np.sum(oy*oy) - R*R
        disc = bq*bq - 4.0*c
        ok = disc >= 0.0
        if not np.any(ok):
            continue
        s = np.sqrt(np.maximum(disc, 0.0))
        t1 = (-bq - s)/2.0
        t2 = (-bq + s)/2.0
        t = np.where((t1>1e-6), t1, np.where((t2>1e-6), t2, np.nan))
        P = oy + dy * t[:,None]
        x, yv, z = P[:,0], P[:,1], P[:,2]
        lon_row = np.degrees(np.arctan2(yv, x))
        hyp = np.hypot(x, yv)
        lat_row = np.degrees(np.arctan2(z, hyp))
        lon[y,:] = lon_row.astype(np.float32)
        lat[y,:] = lat_row.astype(np.float32)

    return lon, lat, mask_inside
