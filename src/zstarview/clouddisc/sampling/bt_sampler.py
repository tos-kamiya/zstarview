import numpy as np
import xarray as xr
from pyproj import Transformer, CRS

def build_bt_sampler(da: xr.DataArray):
    """
    geostationary 固定格子の DataArray(BT[K]) に対して、
    (lon, lat) -> BT[K] を返すベクトル化サンプラを構築する。

    前提:
      - DataArray に satpy が付与する area 定義 (da.attrs["area"]) が存在すること
        (GOES: abi_l2_nc / Himawari: ahi_hsd, ahi_l2_nc)
      - da は 2D（y, x）であること
    方式:
      1) 入力 lon/lat を area の投影系へ変換 (pyproj)
      2) area_extent と shape から最近傍の配列インデックス (iy, ix) を算出
      3) 範囲外は NaN、範囲内は data[iy, ix] を返す
    """
    area = da.attrs.get("area", None)
    if area is None:
        # 互換フォールバック: 2D lon/lat がある場合のみ、簡易最近傍
        if "lon" in da.coords and "lat" in da.coords and da.lon.ndim == 2 and da.lat.ndim == 2:
            lon2 = da.lon.values
            lat2 = da.lat.values
            data = np.asarray(da.values, dtype=np.float32)
            flat = np.column_stack([lon2.ravel(), lat2.ravel()])
            data_flat = data.ravel()

            def sampler(lon, lat):
                xy = np.column_stack([lon.ravel(), lat.ravel()])
                out = np.full(lon.size, np.nan, dtype=np.float32)
                step = 20000
                for i0 in range(0, xy.shape[0], step):
                    blk = xy[i0:i0+step]
                    d2 = (flat[:, 0, None] - blk[:, 0])**2 + (flat[:, 1, None] - blk[:, 1])**2
                    idx = np.argmin(d2, axis=0)
                    out[i0:i0+step] = data_flat[idx]
                return out.reshape(lon.shape)
            return sampler

        raise RuntimeError(
            "BT sampler: 'area' attribute not present and no 2D lon/lat coords; "
            "ensure you open data with Satpy readers (abi_l2_nc / ahi_hsd / ahi_l2_nc)."
        )

    # --- area から投影と格子情報を取得 ---
    # extent の順序は (min_x, min_y, max_x, max_y)
    min_x, min_y, max_x, max_y = area.area_extent
    height, width = area.shape  # (rows, cols)
    dx = (max_x - min_x) / width
    dy = (max_y - min_y) / height

    # DataArray 実体
    data = np.asarray(da.compute().values, dtype=np.float32)

    # pyproj CRS: 新しめの pyresample は area.crs を持つ。無ければ proj_dict から作る。
    target_crs = getattr(area, "crs", None)
    if target_crs is None:
        # proj_dict は PROJ のパラメータ dict
        target_crs = CRS.from_dict(area.proj_dict)

    # lon/lat(deg) -> 投影 (x,y)
    transformer = Transformer.from_crs("EPSG:4326", target_crs, always_xy=True)

    def sampler(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
        # 1) 投影
        x, y = transformer.transform(lon, lat)

        # 2) まず有限判定（投影で NaN/inf になった点は即座に無効化）
        finite = np.isfinite(x) & np.isfinite(y)

        # 3) 浮動小数の画素位置を計算（ここでは丸めない）
        #    ※ まだ整数化しない（NaN を cast しないため）
        ix_f = (x - min_x) / dx
        iy_f = (max_y - y) / dy  # 行0が上なので y を反転

        # 4) 画像範囲内かどうか（浮動小数の段階で判定）
        in_bounds = (
            (ix_f >= 0.0) & (ix_f <= (width - 1)) &
            (iy_f >= 0.0) & (iy_f <= (height - 1))
        )

        # 5) 双一次補間（有効点のみ）
        out = np.full(lon.shape, np.nan, dtype=np.float32)
        valid = finite & in_bounds
        if np.any(valid):
            # 有効点の浮動小数座標
            ixv = ix_f[valid]
            iyv = iy_f[valid]

            ix0 = np.floor(ixv).astype(np.int32)
            iy0 = np.floor(iyv).astype(np.int32)
            ix1 = ix0 + 1
            iy1 = iy0 + 1

            wx = ixv - ix0  # 0..1
            wy = iyv - iy0

            # 範囲外（右端/下端で +1 がはみ出す）を除外
            ok = (
                (ix0 >= 0) & (iy0 >= 0) &
                (ix1 < width) & (iy1 < height)
            )

            # まず最近傍のフォールバックを用意（補間NG点向け）
            ix_nn = np.rint(ixv).astype(np.int32)
            iy_nn = np.rint(iyv).astype(np.int32)
            nn_ok = (
                (ix_nn >= 0) & (ix_nn < width) &
                (iy_nn >= 0) & (iy_nn < height)
            )
            tmp = np.full(ixv.shape, np.nan, dtype=np.float32)
            tmp[nn_ok] = data[iy_nn[nn_ok], ix_nn[nn_ok]]

            if np.any(ok):
                # 4隅の値を引く
                v00 = data[iy0[ok], ix0[ok]]
                v10 = data[iy0[ok], ix1[ok]]
                v01 = data[iy1[ok], ix0[ok]]
                v11 = data[iy1[ok], ix1[ok]]

                # 重み
                wx_ok = wx[ok]; wy_ok = wy[ok]
                w00 = (1 - wx_ok) * (1 - wy_ok)
                w10 = (     wx_ok) * (1 - wy_ok)
                w01 = (1 - wx_ok) * (     wy_ok)
                w11 = (     wx_ok) * (     wy_ok)

                # NaN をもつ画素は重み0にして再正規化
                m00 = np.isfinite(v00); m10 = np.isfinite(v10)
                m01 = np.isfinite(v01); m11 = np.isfinite(v11)
                w00 *= m00; w10 *= m10; w01 *= m01; w11 *= m11
                ws = w00 + w10 + w01 + w11

                # 補間値（ws==0 はそのまま NaN → 後で最近傍が入る）
                interp = (np.where(m00, v00*w00, 0.0) +
                        np.where(m10, v10*w10, 0.0) +
                        np.where(m01, v01*w01, 0.0) +
                        np.where(m11, v11*w11, 0.0))
                # 正規化
                good = ws > 0
                interp[good] /= ws[good]

                # ok のところは補間値で上書き
                tmp[ok] = interp

            out[valid] = tmp
        return out

    return sampler
