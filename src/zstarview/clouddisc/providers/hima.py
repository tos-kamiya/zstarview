from __future__ import annotations
import datetime as dt
from typing import Tuple, List
from pathlib import Path
import re

import boto3
from botocore import UNSIGNED
from botocore.config import Config
import numpy as np
import xarray as xr
from satpy import Scene

from ..config import CloudDiscConfig

_S3CFG = Config(signature_version=UNSIGNED, retries={"max_attempts": 3, "mode": "standard"})
_HIMA_BUCKETS = ["noaa-himawari9", "noaa-himawari8"]

# HSD B13 セグメント: AHI-L1b-FLDK/YYYY/MM/DD/HHMM/HS_H0[89]_YYYYMMDD_HHMM_B13_FLDK_R.._S####.DAT(.bz2)
_RE_HSD = re.compile(r"HS_H0[89]_\d{8}_\d{4}_B13_FLDK_R\d{2}_S\d{4}\.DAT(\.bz2)?$", re.I)

# ISatSS L2 B13: AHI-L2-FLDK-ISatSS/YYYY/MM/DD/HHMM/....nc（簡易）
_RE_IS = re.compile(r".*B13.*\.nc$", re.I)

def _ymdhm10(t: dt.datetime):
    t = t.astimezone(dt.timezone.utc).replace(second=0, microsecond=0)
    return t.year, t.month, t.day, f"{t.hour:02d}{(t.minute//10)*10:02d}"

class HimaProvider:
    def __init__(self, cfg: CloudDiscConfig):
        self.cfg = cfg
        self.root_hsd = cfg.cache_root() / "hima_hsd"
        self.root_is  = cfg.cache_root() / "hima_isatss"
        self.root_hsd.mkdir(parents=True, exist_ok=True)
        self.root_is.mkdir(parents=True, exist_ok=True)
        self.s3 = boto3.client("s3", config=_S3CFG)

    def _download(self, bucket: str, key: str, root: Path) -> Path:
        dst = root / bucket / key
        if dst.exists():
            return dst
        dst.parent.mkdir(parents=True, exist_ok=True)
        tmp = dst.with_suffix(dst.suffix + ".tmp")
        with tmp.open("wb") as f:
            self.s3.download_fileobj(bucket, key, f)
        tmp.replace(dst)
        return dst

    def _list(self, bucket: str, prefix: str) -> list[str]:
        page = self.s3.get_paginator("list_objects_v2").paginate(Bucket=bucket, Prefix=prefix)
        keys = []
        for p in page:
            for obj in p.get("Contents", []):
                keys.append(obj["Key"])
        return keys

    def _find_hsd(self, when_utc: dt.datetime):
        y,m,d,hm = _ymdhm10(when_utc)
        for slot in range(0, self.cfg.search_back_minutes+1, 10):
            t = when_utc - dt.timedelta(minutes=slot)
            y,m,d,hm = _ymdhm10(t)
            for b in _HIMA_BUCKETS:
                prefix = f"AHI-L1b-FLDK/{y:04d}/{m:02d}/{d:02d}/{hm}/"
                keys = [k for k in self._list(b, prefix) if _RE_HSD.search(Path(k).name)]
                if keys:
                    keys.sort()
                    return b, keys, t
        return None, None, None

    def _find_isatss(self, when_utc: dt.datetime):
        y,m,d,hm = _ymdhm10(when_utc)
        for slot in range(0, self.cfg.search_back_minutes+1, 10):
            t = when_utc - dt.timedelta(minutes=slot)
            y,m,d,hm = _ymdhm10(t)
            for b in _HIMA_BUCKETS:
                prefix = f"AHI-L2-FLDK-ISatSS/{y:04d}/{m:02d}/{d:02d}/{hm}/"
                keys = [k for k in self._list(b, prefix) if _RE_IS.search(Path(k).name)]
                if keys:
                    keys.sort()
                    return b, keys[-1], t
        return None, None, None

    def fetch_bt_c13(self, when_utc: dt.datetime):
        """HSD(B13)優先→無ければISatSS(B13)。戻り: (DataArray[K], used_time, [paths], sub_lon)"""
        # HSD (10枚程度)
        b, keys, t = self._find_hsd(when_utc)
        if b and keys:
            paths = [self._download(b, k, self.root_hsd) for k in keys]
            scn = Scene(reader="ahi_hsd", filenames=[str(p) for p in paths])
            scn.load(["B13"], calibration="brightness_temperature")
            da = scn["B13"].astype(np.float32).compute()
            # Himawari サブ衛星経度 ~ 140.7
            return da, t.replace(minute=(t.minute//10)*10, second=0, microsecond=0, tzinfo=dt.timezone.utc), paths, 140.7

        # ISatSS
        b2, key_nc, t2 = self._find_isatss(when_utc)
        if b2 and key_nc:
            path = self._download(b2, key_nc, self.root_is)
            scn = Scene(reader="ahi_l2_nc", filenames=[str(path)])
            scn.load(["B13"])
            da = scn["B13"].astype(np.float32).compute()
            return da, t2.replace(minute=(t2.minute//10)*10, second=0, microsecond=0, tzinfo=dt.timezone.utc), [path], 140.7

        raise RuntimeError("Himawari B13 not found in search window")

