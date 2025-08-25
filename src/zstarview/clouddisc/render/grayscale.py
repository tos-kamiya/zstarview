import warnings
warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in log",
    category=RuntimeWarning,
    module="dask._task_spec",  # Satpyのdask実行パスで出る行をピンポイント抑制
)

import numpy as np
from PIL import Image

def _bt_to_weight(bt, bt_warm, bt_cold):
    # w = 0..1 （暖かい→0, 冷たい→1）
    with np.errstate(invalid="ignore"):
        w = (bt_warm - bt) / max(1e-6, (bt_warm - bt_cold))
    w = np.clip(w, 0.0, 1.0)
    w = np.nan_to_num(w, nan=0.0, posinf=1.0, neginf=0.0)
    return w

def bt_to_LA(
    bt: np.ndarray,
    mask_inside: np.ndarray,
    bt_warm: float, bt_cold: float,
    gamma: float,
    brightness_as_alpha: bool = False,
):
    """BT[K]→グレイスケール(L)と円形アルファ(A)のLA画像（フェードなし）"""
    # 輝度
    w = _bt_to_weight(bt, bt_warm, bt_cold)
    if gamma is not None and gamma != 1.0:
        with np.errstate(invalid="ignore"):
            w = np.power(w, gamma)
        w = np.nan_to_num(w, nan=0.0)
    L = (w * 255.0).astype(np.uint8)

    # アルファ（内側=255, 外側=0）
    A = np.zeros_like(L, dtype=np.uint8)
    A[mask_inside] = 255

    if brightness_as_alpha:
        la = np.dstack([A, L])
    else:
        la = np.dstack([L, A])
    return Image.fromarray(la, mode="LA")
