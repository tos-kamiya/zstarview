# -*- coding: utf-8 -*-
"""
Provides utility for cleaning up cache directories.
"""

from pathlib import Path
from datetime import datetime, timedelta, timezone


def cleanup_satellite_cache(root: Path, *, hours: int = 24, dry_run: bool = False) -> None:
    """
    Clean up satellite cache directories (goes_cmipf and hima_hsd).
    Deletes files older than `hours` while keeping the most recent file
    in each subdirectory and skipping files with an active .inprogress marker.
    """
    now = datetime.now(timezone.utc)
    ttl = timedelta(hours=hours)

    # 衛星データディレクトリを固定
    targets = ["goes_cmipf", "hima_hsd"]

    for kind in targets:
        base = root / kind
        if not base.is_dir():
            continue

        # ディレクトリごとにファイルを集めて新しい順にソート
        per_dir = {}
        for f in base.rglob("*"):
            if f.is_file():
                per_dir.setdefault(f.parent, []).append(f)

        for dir_path, files in per_dir.items():
            files.sort(key=lambda p: p.stat().st_mtime, reverse=True)

            for idx, file_path in enumerate(files):
                # 直近1個は残す
                if idx == 0:
                    continue

                # ダウンロード中は残す
                if file_path.with_suffix(file_path.suffix + ".inprogress").exists():
                    continue

                try:
                    mtime = datetime.fromtimestamp(file_path.stat().st_mtime, tz=timezone.utc)
                except Exception:
                    continue  # 取得失敗時は安全側で残す

                # TTLを超えていたら削除
                if (now - mtime) > ttl:
                    if dry_run:
                        print(f"[dry-run] delete {file_path}")
                    else:
                        try:
                            file_path.unlink()
                            print(f"deleted {file_path}")
                        except OSError as e:
                            print(f"error deleting {file_path}: {e}")

        # 空ディレクトリを削除
        for d in sorted(base.rglob("*"), key=lambda p: len(p.parts), reverse=True):
            if d.is_dir():
                try:
                    if not any(d.iterdir()):
                        if dry_run:
                            print(f"[dry-run] rmdir {d}")
                        else:
                            d.rmdir()
                            print(f"removed empty dir {d}")
                except OSError:
                    pass
