import datetime as dt
from typing import Tuple, List
import numpy as np
from PIL import Image

from .config import CloudDiscConfig
from .types import CloudMeta, VisibilityError, DataNotFoundError, DownloadError, RenderError
from .providers.select import pick_satellite
from .providers.goes import GoesProvider
from .providers.hima import HimaProvider
from .sampling.bt_sampler import build_bt_sampler
from .projectors.az import az_project_lonlat_grid
from .render.grayscale import bt_to_LA

class CloudDisc:
    def __init__(self, cfg: CloudDiscConfig):
        self.cfg = cfg
        self.goes = GoesProvider(cfg)
        self.hima  = HimaProvider(cfg)

    def _now_rounded(self) -> dt.datetime:
        t = dt.datetime.now(dt.timezone.utc).replace(second=0, microsecond=0)
        return t.replace(minute=(t.minute // 10) * 10)

    def render_now(self, lat: float, lon: float, alt: float, az: float, radius_px: int, brightness_as_alpha: bool = False):
        """現在UTC(10分丸め)の雲量LA画像を 2R×2R で返す"""
        when = self._now_rounded()

        EPS = 1e-3
        alt = max(-(90.0 - EPS), min(90.0 - EPS, alt))

        # 1) 衛星の自動選択 or 優先順
        sat = pick_satellite(lat, lon, priority=self.cfg.sat_priority)

        # 2) プロバイダでデータ取得（BT[K]の DataArray と投影情報を得る）
        sat_used = sat
        if sat in ("G16", "G18"):
            res, sat_used = self.goes.fetch_bt_c13_with_failover(sat=sat, when_utc=when)
            da, used_time, src_paths = res
            product = "CMIPF-C13"
            sub_lon = -75.2 if (sat_used == "G16") else -137.0
        elif sat == "HIMAWARI":
            da, used_time, src_paths, sub_lon = self.hima.fetch_bt_c13(when_utc=when)  # sub_lon ~ 140.7
            product = "HSD-B13" if len(src_paths) > 1 else "ISatSS-B13"
        else:
            raise VisibilityError("No suitable satellite")

        # 3) (lon,lat)->BT[K] のサンプラ作成
        sampler = build_bt_sampler(da, sub_lon)

        # 4) AZ投影で (2R×2R) の lon/lat グリッドを作る
        try:
            lon_grid, lat_grid, mask_inside = az_project_lonlat_grid(
                lat0=lat, lon0=lon, alt0=alt, az0=az,
                radius_px=radius_px + 1, cloud_shell_km=6371.0+5.0,  # Earth+5km
                alt_min_deg=self.cfg.alt_min_deg
            )
        except Exception as e:
            raise RenderError(f"Projection failed: {e}") from e

        # 5) サンプル→BT→LA
        #    サンプルは NaN を許容。mask_inside で円内のみ使う
        bt = sampler(lon_grid, lat_grid)  # np.float32 / NaN含む

        img = bt_to_LA(
            bt=bt, mask_inside=mask_inside,
            bt_warm=self.cfg.bt_warm_k, bt_cold=self.cfg.bt_cold_k,
            brightness_as_alpha=brightness_as_alpha,
            gamma=self.cfg.gamma, edge_antialias=self.cfg.edge_antialias
        )

        meta = CloudMeta(
            satellite=sat_used, product=product, time_utc=used_time, src_paths=src_paths
        )
        return img, meta
