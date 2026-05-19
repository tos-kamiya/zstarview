# -*- coding: utf-8 -*-
"""
Provides a simple, filesystem-based caching utility.

This module contains a function to ensure a file exists in a local cache,
fetching it via a callback if it is not present. It handles atomic writes
to prevent corrupted files from being created in case of an interruption.
"""

import os
from pathlib import Path
from typing import Callable


def ensure_local_file(root: Path, bucket: str, key: str, fetch_func: Callable[[object], None]) -> Path:
    """
    Ensures a file exists in the local cache, downloading it if necessary.

    This function provides a simple file-based cache. It checks if the desired file
    (specified by `bucket` and `key`) exists. If it does, the path is returned.
    If not, it calls `fetch_func` to retrieve the file.

    The download is performed atomically by first writing to a temporary file
    and then renaming it to the final destination. This prevents race conditions
    and ensures that the final file is never in a partially downloaded state.

    Args:
        root: The root directory of the cache.
        bucket: The subdirectory within the cache (e.g., a data source name).
        key: The filename for the data.
        fetch_func: A callback function that takes a file-like object (opened in
                    binary write mode) and writes the content of the file to it.

    Returns:
        The path to the cached file.
    """
    dst = root / bucket / key
    if dst.exists():
        return dst

    # Ensure the parent directory for the file exists.
    dst.parent.mkdir(parents=True, exist_ok=True)

    # Download to a temporary file first for atomicity.
    # This prevents other processes from reading a partially written file.
    tmp_path = dst.with_suffix(dst.suffix + ".tmp")
    try:
        with tmp_path.open("wb") as f:
            fetch_func(f)
        # Atomically move the completed file to its final destination.
        tmp_path.replace(dst)
    finally:
        # Clean up the temporary file if it still exists (e.g., on an error).
        if tmp_path.exists():
            try:
                os.remove(tmp_path)
            except OSError:
                pass  # Ignore errors if another process already cleaned it up.

    return dst
