import math
from typing import Tuple

# サブ衛星経度
SAT_LON = {"G16": -75.2, "G18": -137.0, "HIMAWARI": 140.7}

def central_angle_deg(lat_deg: float, lon_deg: float, sub_lon_deg: float) -> float:
    φ = math.radians(lat_deg)
    Δλ = math.radians((lon_deg - sub_lon_deg + 540) % 360 - 180)
    cosψ = max(-1.0, min(1.0, math.cos(φ) * math.cos(Δλ)))
    return math.degrees(math.acos(cosψ))

def pick_satellite(lat: float, lon: float, priority=("AUTO",)) -> str:
    if "AUTO" in priority:
        cand = []
        for sat, slon in SAT_LON.items():
            ψ = central_angle_deg(lat, lon, slon)
            if ψ <= 81.3:  # 可視域
                cand.append((ψ, sat))
        if not cand:
            # 最小角度を採用（視野外でも最も近い）
            cand = [(central_angle_deg(lat, lon, slon), s) for s, slon in SAT_LON.items()]
        cand.sort()
        return cand[0][1]
    # 明示優先順
    for s in priority:
        if s=="AUTO": continue
        return s
    return "HIMAWARI"
