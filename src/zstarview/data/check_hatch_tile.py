import numpy as np
import math
from PIL import Image

def make_hatch_tile_pil(W: int, H: int, line_px: int, strength: int) -> Image.Image:
    """
    ハッチ用のタイル画像を生成し、PIL.Image として返す (グレイスケール).
    """
    norm = math.sqrt(W*W + H*H)
    P = W * H
    band_u = max(1, int(round(line_px * norm)))

    xs = np.arange(W, dtype=np.int32)[None, :]
    ys = np.arange(H, dtype=np.int32)[:, None]

    u = H * xs - W * ys
    u_mod = np.mod(u, P)
    dist = np.minimum(u_mod, P - u_mod)
    mask = dist <= (band_u / 2)

    arr = np.zeros((H, W), dtype=np.uint8)
    arr[mask] = np.uint8(np.clip(strength, 0, 255))

    return Image.fromarray(arr, mode="L")

# タイルを生成して保存
tile = make_hatch_tile_pil(32, 20, 5, 200)
tile.save("hatch_tile_debug.png")
