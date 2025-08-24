import numpy as np
from PIL import Image

def _bt_to_weight(bt, bt_warm, bt_cold):
    # w = 0..1 （暖かい→0, 冷たい→1）
    with np.errstate(invalid="ignore"):
        w = (bt_warm - bt) / max(1e-6, (bt_warm - bt_cold))
    w = np.clip(w, 0.0, 1.0)
    w = np.nan_to_num(w, nan=0.0, posinf=1.0, neginf=0.0)
    return w

def bt_to_LA(bt: np.ndarray, mask_inside: np.ndarray,
             bt_warm: float, bt_cold: float, gamma: float,
             edge_antialias: bool=False):
    """BT[K]→グレイスケール(L)と円形アルファ(A)のLA画像"""
    w = _bt_to_weight(bt, bt_warm, bt_cold)
    if gamma is not None and gamma != 1.0:
        with np.errstate(invalid="ignore"):
            w = np.power(w, gamma)
        w = np.nan_to_num(w, nan=0.0)
    L = (w * 255.0).astype(np.uint8)

    # 本体は不変（不透明）
    A = np.zeros_like(L, dtype=np.uint8)
    if not edge_antialias:
        A[mask_inside] = 255
    else:
        m = mask_inside.astype(np.uint8)
        h, w_ = m.shape
        mp = np.pad(m, 1, mode="constant", constant_values=0)

        # 3x3 ガウシアン近似（整数カーネル）
        # [[1, 2, 1],
        #  [2, 4, 2],
        #  [1, 2, 1]] / 16
        neigh = (
            1*mp[0:h,   0:w_] + 2*mp[0:h,   1:w_+1] + 1*mp[0:h,   2:w_+2] +
            2*mp[1:h+1, 0:w_] + 4*mp[1:h+1, 1:w_+1] + 2*mp[1:h+1, 2:w_+2] +
            1*mp[2:h+2, 0:w_] + 2*mp[2:h+2, 1:w_+1] + 1*mp[2:h+2, 2:w_+2]
        ).astype(np.float32)  # 0..16

        gauss = neigh / 16.0  # 0..1 の連続値

        # 内側は255固定、外側の“接しているピクセル”だけフェザー
        A[mask_inside] = 255
        ring = (~mask_inside) & (gauss > 0)

        # 強調が強すぎるなら MAX_ALPHA を 220 前後に下げると自然
        MAX_ALPHA = 220
        A[ring] = np.minimum(gauss[ring] * MAX_ALPHA, 255.0).astype(np.uint8)

        la = np.dstack([L, A])
    return Image.fromarray(la, mode="LA")
