"""
Provides utility for cleaning up cache directories.
"""
import logging
from datetime import datetime, timedelta, timezone
from pathlib import Path

logger = logging.getLogger(__name__)


def cleanup_satellite_cache(root: Path, *, hours: int = 24, dry_run: bool = False) -> None:
    """
    Clean up satellite data cache directories.

    This function removes old files from specified cache subdirectories
    (goes_cmipf, hima_hsd) to free up space.

    The cleanup logic is as follows:
    - It targets files older than a specified number of hours (`ttl`).
    - In each directory, it always preserves the most recently modified file,
      regardless of its age.
    - It skips any file that has a corresponding `.inprogress` file, indicating
      it is part of an ongoing download.
    - After deleting files, it removes any empty subdirectories.

    Args:
        root (Path): The root directory of the cache.
        hours (int, optional): The time-to-live for cache files in hours.
                               Files older than this will be deleted. Defaults to 24.
        dry_run (bool, optional): If True, print the actions that would be taken
                                  without actually deleting anything. Defaults to False.
    """
    now = datetime.now(timezone.utc)
    ttl = timedelta(hours=hours)

    # Target specific satellite data directories
    targets = ["goes_cmipf", "hima_hsd"]

    for kind in targets:
        base = root / kind
        if not base.is_dir():
            continue

        # Group files by parent directory
        per_dir: dict[Path, list[Path]] = {}
        for f in base.rglob("*"):
            if f.is_file():
                per_dir.setdefault(f.parent, []).append(f)

        for dir_path, files in per_dir.items():
            # Sort files by modification time, newest first
            files.sort(key=lambda p: p.stat().st_mtime, reverse=True)

            for idx, file_path in enumerate(files):
                # Always keep the most recent file in each directory
                if idx == 0:
                    continue

                # Skip files currently being downloaded (marked with .inprogress)
                if file_path.with_suffix(file_path.suffix + ".inprogress").exists():
                    continue

                try:
                    mtime = datetime.fromtimestamp(file_path.stat().st_mtime, tz=timezone.utc)
                except Exception:
                    # Failsafe: if mtime cannot be read, keep the file
                    continue

                # If the file is older than the TTL, delete it
                if (now - mtime) > ttl:
                    if dry_run:
                        print(f"[dry-run] delete {file_path}")
                    else:
                        try:
                            file_path.unlink()
                            logger.info(f"deleted {file_path}")
                        except OSError as e:
                            logger.warning(f"error deleting {file_path}: {e}")

        # Clean up empty directories
        for d in sorted(base.rglob("*"), key=lambda p: len(p.parts), reverse=True):
            if d.is_dir():
                try:
                    if not any(d.iterdir()):
                        if dry_run:
                            print(f"[dry-run] rmdir {d}")
                        else:
                            d.rmdir()
                            logger.info(f"removed empty dir {d}")
                except OSError:
                    # Ignore errors when removing directories, as other processes
                    # might be accessing them.
                    pass
