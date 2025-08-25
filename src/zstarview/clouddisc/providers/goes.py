# goes.py
from botocore import UNSIGNED
from botocore.config import Config
import boto3, datetime as dt
from pathlib import Path
import xarray as xr
from satpy import Scene

_GOES_BUCKET = {"G16": "noaa-goes16", "G18": "noaa-goes18"}
_GOES_REGION = {"noaa-goes16": "us-east-1", "noaa-goes18": "us-west-2"}

def _doy(dt_utc: dt.datetime) -> int:
    return dt_utc.timetuple().tm_yday

class GoesProvider:
    def __init__(self, cfg):
        self.cfg = cfg
        self.root = (cfg.cache_root() / "goes_cmipf")
        self.root.mkdir(parents=True, exist_ok=True)
        self._list_cache = {}
        # TODO: self._list_cacheを定期的に（24hくらいで）クリーンアップする処理をつける

    def _s3(self, bucket: str):
        return boto3.client(
            "s3",
            region_name=_GOES_REGION[bucket],
            config=Config(signature_version=UNSIGNED, retries={"max_attempts": 3})
        )

    def _list_hour(self, bucket: str, t: dt.datetime) -> list[str]:
        key = (bucket, t.year, _doy(t), t.hour)
        if key in self._list_cache:
            return self._list_cache[key]

        prefix = f"ABI-L2-CMIPF/{t.year:04d}/{_doy(t):03d}/{t.hour:02d}/"
        s3 = self._s3(bucket)
        print(f"[DBG] client_region={s3.meta.region_name}")
        print(f"[DBG] listing s3://{bucket}/{prefix}")
        page = s3.get_paginator("list_objects_v2").paginate(Bucket=bucket, Prefix=prefix)
        keys: list[str] = []
        for p in page:
            for obj in p.get("Contents", []) or []:
                keys.append(obj["Key"])
        print(f"[DBG] objects={len(keys)}")
        self._list_cache[key] = keys
        return keys

    def _download(self, bucket: str, key: str) -> Path:
        dst = self.root / bucket / key
        if dst.exists():
            return dst
        dst.parent.mkdir(parents=True, exist_ok=True)
        tmp = dst.with_suffix(dst.suffix + ".tmp")
        s3 = self._s3(bucket)
        with tmp.open("wb") as f:
            s3.download_fileobj(bucket, key, f)
        tmp.replace(dst)
        return dst

    def _fetch_bt_c13_once(self, sat: str, when_utc: dt.datetime, search_back_minutes: int):
        """単一衛星・指定の検索窓でC13を探す（見つかれば (da, used_time, [path]) を返す）"""
        bucket = _GOES_BUCKET[sat]
        base = when_utc
        for mback in range(0, search_back_minutes + 1, 10):
            t = base - dt.timedelta(minutes=mback)
            keys = self._list_hour(bucket, t)
            if not keys:
                continue
            # C13 のみ抽出
            keys_c13 = [k for k in keys if ("-M6C13_" in k or "-C13_" in k)]
            if not keys_c13:
                continue
            keys_c13.sort()
            key = keys_c13[-1]
            print(f"[DBG] pick C13: {Path(key).name}")
            path = self._download(bucket, key)
            try:
                scn = Scene(reader="abi_l2_nc", filenames=[str(path)])
                scn.load(["C13"])
                da = scn["C13"].astype("float32").compute()
                used_time = t.replace(minute=(t.minute//10)*10, second=0, microsecond=0, tzinfo=dt.timezone.utc)
                return da, used_time, [path]
            except Exception as e:
                print(f"[WARN] satpy load failed for {Path(key).name}: {e}")
                # フォールバック（変数名は CMI のことが多い）
                with xr.open_dataset(path, engine="netcdf4") as ds:
                    if "CMI" in ds.variables:
                        da = ds["CMI"].astype("float32").compute()
                        used_time = t.replace(minute=(t.minute//10)*10, second=0, microsecond=0, tzinfo=dt.timezone.utc)
                        return da, used_time, [path]
                # 読めなければ次の候補へ
        return None  # 見つからず

    def fetch_bt_c13_with_failover(self, sat: str, when_utc: dt.datetime, extra_back_minutes: int = 30):
        """satを第一候補として試し、失敗時にもう片方へフェイルオーバー。
           それでもダメなら検索窓を extra_back_minutes だけ広げて再試行する。"""
        primary = sat
        secondary = "G18" if sat == "G16" else "G16"

        # 1st pass: 指定の検索窓
        res = self._fetch_bt_c13_once(primary, when_utc, self.cfg.search_back_minutes)
        if res:
            return res, primary
        print(f"[INFO] primary {primary} no C13; trying {secondary}")
        res = self._fetch_bt_c13_once(secondary, when_utc, self.cfg.search_back_minutes)
        if res:
            return res, secondary

        # 2nd pass: 検索窓を少し広げる
        widen = self.cfg.search_back_minutes + extra_back_minutes
        print(f"[INFO] widen search window to {widen} minutes")
        res = self._fetch_bt_c13_once(primary, when_utc, widen)
        if res:
            return res, primary
        res = self._fetch_bt_c13_once(secondary, when_utc, widen)
        if res:
            return res, secondary

        raise RuntimeError("GOES CMIPF C13 not found (after failover and widened window)")

