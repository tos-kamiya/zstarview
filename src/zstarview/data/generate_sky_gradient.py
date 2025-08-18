import math
import numpy as np
from PIL import Image

# --- ▼ パラメータ設定 ▼ ---
IMG_WIDTH = 300
IMG_HEIGHT = 181
SUN_ALTITUDE_DEG = 60.0   # 太陽の高度 (-5度で日没、90度で天頂)
EXPOSURE = 1.0          # 全体の明るさ調整
SATURATION = 1.2        # 色の鮮やかさ
GROUND_COLOR = (0.1, 0.1, 0.1) # 地面の色（地平線下の色）
OUTPUT_FILENAME = "sky_gradient_chart.png"
# --- ▲ パラメータ設定 ▲ ---


# --- ▼ ヘルパー関数 ▼ ---
def _deg2rad(d): return d * math.pi / 180.0
def _clamp01(x): return max(0.0, min(1.0, x))

def _angle_between(alt1_deg, az1_deg, alt2_deg, az2_deg) -> float:
    a1, z1 = _deg2rad(alt1_deg), _deg2rad(az1_deg)
    a2, z2 = _deg2rad(alt2_deg), _deg2rad(az2_deg)
    cos_g = (math.sin(a1)*math.sin(a2) + math.cos(a1)*math.cos(a2)*math.cos(z2 - z1))
    cos_g = max(-1.0, min(1.0, cos_g))
    return math.acos(cos_g)

# 色を混ぜる関数
def _lerp_color(c1, c2, t):
    t = _clamp01(t)
    return (
        c1[0] * (1-t) + c2[0] * t,
        c1[1] * (1-t) + c2[1] * t,
        c1[2] * (1-t) + c2[2] * t,
    )

# --- ▼ 現象模倣モデル ▼ ---

def get_sun_color(sun_alt_deg: float) -> tuple[float, float, float]:
    """
    ステップ1: 太陽の高度に基づいて太陽光の色を決める
    """
    # 色を定義
    zenith_color = (0.3, 0.5, 1.0)  # 天頂にあるときの色 (青)
    horizon_color = (1.0, 0.8, 0.4) # 地平線にあるときの色 (オレンジ)
    night_color = (0.01, 0.02, 0.05) # 夜の色 (濃紺)

    # 太陽高度 -5度(日没) から 90度(天頂) の範囲を 0-1 に正規化
    t = _clamp01((sun_alt_deg + 5.0) / 95.0)

    # 昼の色 (地平線〜天頂)
    day_color = _lerp_color(horizon_color, zenith_color, t)
    
    # 昼と夜の色を混ぜる (太陽が地平線近くで急激に暗くなるのを表現)
    fade = _clamp01((sun_alt_deg) / 10.0) # 0〜10度でフェード
    
    return _lerp_color(night_color, day_color, fade)


def get_sky_color(view_alt_deg: float, view_az_deg: float, sun_alt_deg: float, sun_az_deg: float):
    """
    ステップ2: 視線と太陽の角度に基づいて最終的な空の色を決める
    """
    # 1. 基本となる太陽光の色を取得
    sun_color = get_sun_color(sun_alt_deg)
    
    # 夜は単純に暗くする
    if sun_alt_deg < -5:
        return (0, 0, 0)

    # 2. 視線と太陽のなす角から、明るさ(不透明度)を決める
    gamma_rad = _angle_between(view_alt_deg, view_az_deg, sun_alt_deg, sun_az_deg)
    
    # 太陽方向を1, 反対側を0とする
    # cos(gamma) は -1〜1 なので、(1+cos)/2 で 0〜1 に変換
    brightness = (1.0 + math.cos(gamma_rad)) / 2.0
    
    # 太陽の周りをより明るく強調する (べき乗でカーブを調整)
    brightness = brightness ** 2.0

    # 3. 視線の高度による色の変化
    # 天頂は暗く(色が濃く)、地平線は明るく(白っぽく)なる効果
    horizon_fade = _clamp01(view_alt_deg / 90.0)
    zenith_darkness = 1.0 - (1.0 - horizon_fade) * 0.5 # 天頂を50%暗くする
    horizon_whiteness = (1.0 - horizon_fade) * 0.3 # 地平線を30%白っぽくする
    
    # 4. 全てを合成
    final_color = (
        sun_color[0] * brightness * zenith_darkness + horizon_whiteness,
        sun_color[1] * brightness * zenith_darkness + horizon_whiteness,
        sun_color[2] * brightness * zenith_darkness + horizon_whiteness,
    )
    
    return final_color

# --- ▼ メイン処理 ▼ ---
def main():
    print("空のグラデーション画像を生成します...")
    print(f"設定: 太陽高度={SUN_ALTITUDE_DEG}°")

    img = Image.new('RGB', (IMG_WIDTH, IMG_HEIGHT))
    sun_az = 0.0
    view_azimuths = {"sun": sun_az, "side": 90.0, "opposite": 180.0}

    for y in range(IMG_HEIGHT):
        view_alt = 90.0 - y

        for i, (name, view_az) in enumerate(view_azimuths.items()):
            
            # 現象模倣モデルで色を取得
            rgb = get_sky_color(view_alt, view_az, SUN_ALTITUDE_DEG, sun_az)
            
            # 彩度と露出を調整
            # (簡易的な彩度調整: グレースケールとの混合)
            gray = rgb[0] * 0.299 + rgb[1] * 0.587 + rgb[2] * 0.114
            r = _lerp_color((gray,gray,gray), rgb, SATURATION)[0]
            g = _lerp_color((gray,gray,gray), rgb, SATURATION)[1]
            b = _lerp_color((gray,gray,gray), rgb, SATURATION)[2]
            
            r = _clamp01(r * EXPOSURE)
            g = _clamp01(g * EXPOSURE)
            b = _clamp01(b * EXPOSURE)

            # 地平線より下は地面の色にフェード
            if view_alt < 0:
                t = _clamp01((view_alt + 5.0) / 5.0)
                r, g, b = _lerp_color(GROUND_COLOR, (r, g, b), t)

            color_int = (int(r * 255), int(g * 255), int(b * 255))
            
            strip_width = IMG_WIDTH // 3
            for x_offset in range(strip_width):
                x = i * strip_width + x_offset
                img.putpixel((x, y), color_int)

    img.save(OUTPUT_FILENAME)
    print(f"画像を '{OUTPUT_FILENAME}' として保存しました。")

if __name__ == '__main__':
    main()