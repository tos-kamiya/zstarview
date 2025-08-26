"""
Utilities for cleaning up the cache directory.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Callable, Optional, Pattern


@dataclass
class CleanupSpec:
    """
    Specifies the rules for cleaning up a certain kind of data in the cache.

    Attributes:
        kind: The type of data, corresponding to a subdirectory in the cache root.
        keep_globs: A list of glob patterns. Files matching any of these will be kept.
        keep_regexes: A list of compiled regular expressions. Files with paths matching
                      any of these will be kept.
        keep_pred: An optional function that takes a Path object and returns True if
                   the file should be kept.
        skip_if_inprogress: If True, skips cleanup for a file if a corresponding
                            '.inprogress' file exists.
    """

    kind: str
    keep_globs: List[str]
    keep_regexes: List[Pattern] = field(default_factory=list)
    keep_pred: Optional[Callable[[Path], bool]] = None
    skip_if_inprogress: bool = True


@dataclass
class CleanupReport:
    """
    A report summarizing the results of a cache cleanup operation.

    Attributes:
        scanned: The total number of files scanned.
        kept: The number of files that were kept.
        deleted: The number of files that were deleted.
        emptied_dirs: The number of empty directories that were removed.
        errors: The number of errors encountered during cleanup.
        dry_run: True if the cleanup was a dry run (no actual deletions).
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
        dry_run: If True, simulates the cleanup without actually deleting anything.

    Returns:
        A CleanupReport object summarizing the operation.
    """
    scanned = kept = deleted = emptied = errors = 0
    kinds = {s.kind: s for s in specs}
    for kind, spec in kinds.items():
        base = root / kind
        if not base.exists():
            continue

        # --- File Deletion Phase ---
        all_files = [p for p in base.rglob("*") if p.is_file()]
        for file_path in all_files:
            scanned += 1
            should_keep = False

            # Rule 1: Skip if an '.inprogress' file exists
            if spec.skip_if_inprogress and file_path.with_suffix(file_path.suffix + ".inprogress").exists():
                should_keep = True

            # Rule 2: Check against glob patterns
            if not should_keep:
                for glob_pattern in spec.keep_globs:
                    if file_path.match(glob_pattern):
                        should_keep = True
                        break

            # Rule 3: Check against regex patterns
            if not should_keep and spec.keep_regexes:
                for regex_pattern in spec.keep_regexes:
                    if regex_pattern.search(str(file_path)):
                        should_keep = True
                        break

            # Rule 4: Check against predicate function
            if not should_keep and spec.keep_pred:
                try:
                    if spec.keep_pred(file_path):
                        should_keep = True
                except Exception:
                    errors += 1  # Error in predicate function

            # Perform action based on the decision
            if should_keep:
                kept += 1
            else:
                # Delete the file
                try:
                    if not dry_run:
                        file_path.unlink(missing_ok=True)
                    deleted += 1
                except Exception:
                    errors += 1

        # --- Empty Directory Cleanup Phase ---
        # Iterate from deepest to shallowest to remove empty directories
        all_dirs = sorted([d for d in base.rglob("*") if d.is_dir()], key=lambda p: len(p.parts), reverse=True)

        for dir_path in all_dirs:
            if dir_path.is_dir():
                try:
                    # Check if the directory is empty
                    if not any(dir_path.iterdir()):
                        if not dry_run:
                            dir_path.rmdir()
                        emptied += 1
                except Exception:
                    # This can happen if the directory is deleted in another process, etc.
                    errors += 1

    return CleanupReport(scanned, kept, deleted, emptied, errors, dry_run)
