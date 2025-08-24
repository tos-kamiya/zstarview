# clouddisc/providers/goes.py
import datetime as dt
from typing import Tuple, List
from pathlib import Path
import re

import boto3
from botocore import UNSIGNED
from botocore.config import Config
from satpy import Scene
import numpy as np
import xarray as xr

from ..config import CloudDiscConfig

_S3CFG = Config(signature_version=UNSIGNED, retries={"max_attempts": 3, "mode": "standard"})

_GOES_BUCKET = {"G16": "noaa-goes16", "G18": "noaa-goes18"}
_RE_NC = re.compile(r".*_C13_G1[68]_s\d{13}_e\d{13}_c\d{13}\.nc$")

def _doy(dt_utc: dt.datetime) -> int:
    return dt_utc.timetuple().tm_yday

class GoesProvider:
    def __init__(self, cfg: CloudDiscConfig):
        self.cfg = cfg
        self.root = (cfg.cache_root() / "goes_cmipf")
        self.root.mkdir(parents=True, exist_ok=True)
        self.s3 = boto3.client("s3", config=_S3CFG)

    def _list_hour(self, bucket: str, t: dt.datetime) -> list[str]:
        prefix = f"ABI-L2-CMIPF/{t.year:04d}/{_doy(t):03d}/{t.hour:02d}/"
        page = self.s3.get_paginator("list_objects_v2").paginate(Bucket=bucket, Prefix=prefix)
        keys: list[str] = []
        for p in page:
            for obj in p.get("Contents", []):
                k = obj["Key"]
                if _RE_NC.match(Path(k).name):
                    keys.append(k)
        return keys

    def _download(self, bucket: str, key: str) -> Path:
        dst = self.root / bucket / key
        if dst.exists():
            return dst
        dst.parent.mkdir(parents=True, exist_ok=True)
        tmp = dst.with_suffix(dst.suffix + ".tmp")
        with tmp.open("wb") as f:
            self.s3.download_fileobj(bucket, key, f)
        tmp.replace(dst)
        return dst

    def fetch_bt_c13(self, sat: str, when_utc: dt.datetime):
        """
        CMIPF/C13 を Satpy で開き、lon/lat 付き DataArray を返す。
        戻り値: (xarray.DataArray[np.float32 K], used_time(10min丸めUTC), [Path])
        """
        bucket = _GOES_BUCKET[sat]
        base = when_utc

        for mback in range(0, self.cfg.search_back_minutes + 1, 10):
            t = base - dt.timedelta(minutes=mback)
            keys = self._list_hour(bucket, t)
            if not keys:
                continue

            keys.sort()                # 最新 c... を末尾から
            key = keys[-1]
            path = self._download(bucket, key)

            # ★ Satpy で読み込み（abi_l2_nc）
            try:
                scn = Scene(reader="abi_l2_nc", filenames=[str(path)])
                # Dataset 名は "C13"（輝度温度）。Satpyが lon/lat 座標を付与してくれる。
                scn.load(["C13"])
                da = scn["C13"].astype(np.float32).compute()
                used_time = t.replace(minute=(t.minute//10)*10, second=0, microsecond=0, tzinfo=dt.timezone.utc)
                return da, used_time, [path]
            except Exception:
                # 万一 Satpy が読めなければ xarray フォールバック（この場合 lon/lat が無いので後段で失敗する）
                ds = xr.open_dataset(path, engine="netcdf4")
                for var in ("CMI_C13", "CMI"):
                    if var in ds:
                        da = ds[var].astype(np.float32)
                        ds.close()
                        used_time = t.replace(minute=(t.minute//10)*10, second=0, microsecond=0, tzinfo=dt.timezone.utc)
                        return da, used_time, [path]
                ds.close()
                continue

        raise RuntimeError("GOES CMIPF C13 not found in search window")
