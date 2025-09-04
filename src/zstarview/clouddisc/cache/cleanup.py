# -*- coding: utf-8 -*-
"""
Provides a generic, rule-based utility for cleaning up cache directories.

This module defines a flexible cache cleaning mechanism. You can define
`CleanupSpec` objects to specify which files to keep or delete based on
glob patterns, regular expressions, or custom predicate functions. It is used
to prevent the satellite data cache from growing indefinitely.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Callable, Optional, Pattern


@dataclass
class CleanupSpec:
    """
    Defines a set of rules for cleaning up a specific subdirectory in the cache.

    Attributes:
        kind: The name of the subdirectory in the cache root to which this spec applies.
        keep_globs: A list of glob patterns. Files matching any of these will be kept.
        keep_regexes: A list of compiled regex patterns to match against the full file path.
        keep_pred: An optional function that takes a Path and returns True to keep the file.
        skip_if_inprogress: If True, skips cleanup for a file if a corresponding
                            `.inprogress` file exists, indicating a download is active.
    """

    kind: str
    keep_globs: List[str]
    keep_regexes: List[Pattern] = field(default_factory=list)
    keep_pred: Optional[Callable[[Path], bool]] = None
    skip_if_inprogress: bool = True


@dataclass
class CleanupReport:
    """
    Summarizes the results of a cache cleanup operation.

    Attributes:
        scanned: The total number of files scanned.
        kept: The number of files that were kept.
        deleted: The number of files that were deleted.
        emptied_dirs: The number of empty directories that were removed.
        errors: The number of errors encountered during the process.
        dry_run: True if no actual file deletions were performed.
    """

    scanned: int
    kept: int
    deleted: int
    emptied_dirs: int
    errors: int
    dry_run: bool


def cleanup_cache(root: Path, specs: List[CleanupSpec], *, dry_run=False) -> CleanupReport:
    """
    Cleans up the cache directory based on a list of cleanup specifications.

    Args:
        root: The root directory of the cache.
        specs: A list of CleanupSpec objects defining the cleanup rules.
        dry_run: If True, simulates the cleanup without actually deleting files or directories.

    Returns:
        A CleanupReport object summarizing the results of the operation.
    """
    scanned = kept = deleted = emptied = errors = 0
    kinds = {s.kind: s for s in specs}

    for kind, spec in kinds.items():
        base = root / kind
        if not base.is_dir():
            continue

        # --- Phase 1: File Deletion ---
        # Iterate through all files in the subdirectory.
        all_files = [p for p in base.rglob("*") if p.is_file()]
        for file_path in all_files:
            scanned += 1
            should_keep = False

            # Apply rules in order to decide whether to keep the file.
            if spec.skip_if_inprogress and file_path.with_suffix(file_path.suffix + ".inprogress").exists():
                should_keep = True

            if not should_keep and any(file_path.match(p) for p in spec.keep_globs):
                should_keep = True

            if not should_keep and spec.keep_regexes and any(r.search(str(file_path)) for r in spec.keep_regexes):
                should_keep = True

            if not should_keep and spec.keep_pred:
                try:
                    if spec.keep_pred(file_path):
                        should_keep = True
                except Exception:
                    errors += 1

            # Perform the final action.
            if should_keep:
                kept += 1
            else:
                try:
                    if not dry_run:
                        file_path.unlink(missing_ok=True)
                    deleted += 1
                except OSError:
                    errors += 1

        # --- Phase 2: Empty Directory Cleanup ---
        # Iterate from the deepest directories upwards to remove empty ones.
        all_dirs = sorted([d for d in base.rglob("*") if d.is_dir()], key=lambda p: len(p.parts), reverse=True)
        for dir_path in all_dirs:
            try:
                if not any(dir_path.iterdir()):
                    if not dry_run:
                        dir_path.rmdir()
                    emptied += 1
            except OSError:
                errors += 1

    return CleanupReport(scanned, kept, deleted, emptied, errors, dry_run)
