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

    This function removes old files from the GOES and Himawari cache
    subdirectories (goes_cmipf, hima_isatss) to free up space.

    The cleanup logic is as follows:
    - It targets files older than a specified number of hours (`ttl`).
    - In each satellite/product cache tree, it preserves the most recently
      modified file, regardless of its age.
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
    targets = ["goes_cmipf", "hima_isatss"]

    for kind in targets:
        base = root / kind
        if not base.is_dir():
            continue

        # Cache paths are partitioned by bucket and product before the
        # date/time components. Group at that stable level so a one-file
        # timestamp directory does not exempt every old file from cleanup.
        per_stream: dict[tuple[str, ...], list[Path]] = {}
        for f in base.rglob("*"):
            if f.is_file():
                relative_parts = f.relative_to(base).parts
                stream_key = relative_parts[:2]
                per_stream.setdefault(stream_key, []).append(f)

        for files in per_stream.values():
            # Snapshot modification times once. A download may atomically rename
            # its temporary file while cleanup is running, so paths from rglob()
            # can disappear before stat().
            files_with_mtime: list[tuple[Path, float]] = []
            for file_path in files:
                try:
                    files_with_mtime.append((file_path, file_path.stat().st_mtime))
                except OSError:
                    # Keep cleanup resilient to files removed or made unreadable
                    # by concurrent download activity.
                    continue

            files_with_mtime.sort(key=lambda item: item[1], reverse=True)

            for idx, (file_path, mtime_seconds) in enumerate(files_with_mtime):
                # Always keep the most recent file in each satellite/product stream
                if idx == 0:
                    continue

                # Skip files currently being downloaded (marked with .inprogress)
                if file_path.with_suffix(file_path.suffix + ".inprogress").exists():
                    continue

                mtime = datetime.fromtimestamp(mtime_seconds, tz=timezone.utc)

                # If the file is older than the TTL, delete it
                if (now - mtime) > ttl:
                    if dry_run:
                        print(f"[dry-run] delete {file_path}")
                    else:
                        try:
                            file_path.unlink()
                            logger.debug("deleted %s", file_path)
                        except FileNotFoundError:
                            # Another download or cleanup may have removed it.
                            continue
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
                            logger.debug("removed empty dir %s", d)
                except OSError:
                    # Ignore errors when removing directories, as other processes
                    # might be accessing them.
                    pass
