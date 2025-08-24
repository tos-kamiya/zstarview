from dataclasses import dataclass
from pathlib import Path
from typing import List, Callable, Optional, Pattern
import re
import os

@dataclass
class CleanupSpec:
    kind: str
    keep_globs: List[str]
    keep_regexes: List[Pattern] = None
    keep_pred: Optional[Callable[[Path], bool]] = None
    skip_if_inprogress: bool = True

@dataclass
class CleanupReport:
    scanned: int
    kept: int
    deleted: int
    emptied_dirs: int
    errors: int
    dry_run: bool

def cleanup_cache(root: Path, specs: List[CleanupSpec], *, dry_run=False) -> CleanupReport:
    scanned=kept=deleted=emptied=errors=0
    kinds = {s.kind: s for s in specs}
    for kind, spec in kinds.items():
        base = root / kind
        if not base.exists(): continue
        # 収集
        for p in base.rglob("*"):
            if not p.is_file(): continue
            scanned += 1
            name = str(p)
            if spec.skip_if_inprogress and (p.with_suffix(p.suffix + ".inprogress")).exists():
                kept += 1; continue
            keep = False
            for g in spec.keep_globs:
                if p.match(g):
                    keep = True; break
            if not keep and spec.keep_regexes:
                for rg in spec.keep_regexes:
                    if rg.search(name):
                        keep = True; break
            if not keep and spec.keep_pred:
                try:
                    keep = bool(spec.keep_pred(p))
                except Exception:
                    pass
            if keep:
                kept += 1; continue
            # 削除
            try:
                if not dry_run:
                    p.unlink(missing_ok=True)
                deleted += 1
            except Exception:
                errors += 1
        # 空ディレクトリ掃除
        for d in sorted(base.rglob("*"), key=lambda x: len(str(x)), reverse=True):
            if d.is_dir():
                try:
                    if not any(d.iterdir()):
                        if not dry_run:
                            d.rmdir()
                        emptied += 1
                except Exception:
                    pass
    return CleanupReport(scanned, kept, deleted, emptied, errors, dry_run)
