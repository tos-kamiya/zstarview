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

def az_project_lonlat_grid(lat0, lon0, alt0, az0, radius_px, cloud_shell_km, alt_min_deg=-2.0):
    """AZ投影で画素中心に対応する地表(雲殻)上の lon/lat を返す(MVP)"""
    S = 2*radius_px
    # 画素中心座標
    xs = (np.arange(S)+0.5) - radius_px
    ys = radius_px - (np.arange(S)+0.5)
    X, Y = np.meshgrid(xs, ys)
    r = np.hypot(X, Y)
    # 円内マスク（境界規則 r <= R-0.5）
    mask_inside = r <= (radius_px - 0.5)

    # 視線ベクトル
    v0 = azalt_to_dir_ecef(az0, alt0, lat0, lon0)
    e_obs, n_obs, u_obs = enu_basis(lat0, lon0)
    u_obs = np.array(u_obs)
    O = geodetic_to_ecef(lat0, lon0)

    # 画素→方向
    # 90°に対応する画素半径 R なので、中心→縁で 0→90 deg
    delta_deg = 90.0 * (r / radius_px)
    psi = np.arctan2(Y, X)  # 0→+X

    # 直交基底（接平面）
    def proj_tan(a):
        t = a - np.dot(a, v0)*v0
        n = np.linalg.norm(t) or 1.0
        return t/n
    tY = proj_tan(u_obs)
    tX = np.cross(v0, tY); tX /= (np.linalg.norm(tX) or 1.0)

    cd = np.cos(np.deg2rad(delta_deg)); sd = np.sin(np.deg2rad(delta_deg))
    b = (np.cos(psi)[...,None]*tX + np.sin(psi)[...,None]*tY)  # (H,W,3)
    d = cd[...,None]*v0 + sd[...,None]*b
    d = d / (np.linalg.norm(d, axis=2, keepdims=True)+1e-12)

    # 地平線下を除外
    alt_rad = np.arcsin(np.dot(d, u_obs))
    mask_inside &= (np.degrees(alt_rad) >= alt_min_deg)

    # 雲殻との交点
    O3 = O[None,None,:]
    # 解をベクトル化で解くのは式が長いのでループ版に落とす（MVP）
    H, W = S, S
    lon = np.full((H,W), np.nan, dtype=np.float32)
    lat = np.full((H,W), np.nan, dtype=np.float32)

    R = cloud_shell_km
    for y in range(H):
        dy = d[y]
        oy = O3[0,0]
        # 行単位で計算
        b = 2.0*np.sum(oy*dy, axis=1)
        c = np.sum(oy*oy) - R*R
        disc = b*b - 4.0*c
        ok = disc >= 0.0
        if not np.any(ok): 
            continue
        s = np.sqrt(np.maximum(disc, 0.0))
        t1 = (-b - s)/2.0
        t2 = (-b + s)/2.0
        t = np.where((t1>1e-6), t1, np.where((t2>1e-6), t2, np.nan))
        P = oy + dy * t[:,None]
        # ECEF→緯度経度
        x, yv, z = P[:,0], P[:,1], P[:,2]
        lon_row = np.degrees(np.arctan2(yv, x))
        hyp = np.hypot(x, yv)
        lat_row = np.degrees(np.arctan2(z, hyp))
        lon[y,:] = lon_row.astype(np.float32)
        lat[y,:] = lat_row.astype(np.float32)

    return lon, lat, mask_inside
