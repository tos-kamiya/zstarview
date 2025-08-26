"""
Filesystem-based caching utilities.
"""

from pathlib import Path
from typing import Callable
import os


def ensure_local_file(root: Path, bucket: str, key: str, fetch_func: Callable[[object], None]) -> Path:
    """
    Ensures a file exists in the local cache, downloading it if necessary.

    This function provides a simple file-based cache. It checks if the desired file
    (specified by `bucket` and `key`) exists. If it does, the path is returned.
    If not, it calls `fetch_func` to retrieve the file.

    The download is performed atomically by first writing to a temporary file
    and then renaming it to the final destination.

    Args:
        root: The root directory of the cache.
        bucket: The "bucket" or subdirectory within the cache.
        key: The "key" or filename for the data.
        fetch_func: A callback function that takes a file-like object (opened in binary
                    write mode) and writes the content of the file to it.

    Returns:
        The path to the cached file.
    """
    dst = root / bucket / key
    if dst.exists():
        return dst

    # Ensure parent directory exists
    dst.parent.mkdir(parents=True, exist_ok=True)

    # Download to a temporary file first for atomicity
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    try:
        with tmp.open("wb") as f:
            fetch_func(f)
        # Atomically move the file to the final destination
        tmp.replace(dst)
    finally:
        # Clean up the temporary file if it still exists (e.g., on error)
        if tmp.exists():
            os.remove(tmp)

    return dst
