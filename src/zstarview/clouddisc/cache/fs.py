from pathlib import Path
from typing import Optional
import os

def ensure_local_file(root: Path, bucket: str, key: str, fetch_func) -> Path:
    """
    簡易キャッシュ: 存在すれば返す。無ければ fetch_func で取得して atomic replace。
    fetch_func(fp) は fp(=一時ファイルパス) に書き込むコールバック。
    """
    dst = root / bucket / key
    if dst.exists():
        return dst
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    with tmp.open("wb") as f:
        fetch_func(f)
    tmp.replace(dst)
    return dst
