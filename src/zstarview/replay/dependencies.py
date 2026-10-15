"""Conservative read-only checks of local rendering dependencies."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from ..paths import EPHEMERIS_FILENAME, STARS_CSV_FILE


def check_dependencies(
    root: Path, args: Any, issues: list[str]
) -> tuple[list[dict[str, str]], list[tuple[datetime, datetime]] | None]:
    checks = [
        {
            "name": "observing-height",
            "status": "unverified" if issues else "metadata-available",
            "reason": ", ".join(issues) or "Local observing height resolved",
        }
    ]
    checks.append(
        {
            "name": "star-catalog",
            "status": "unverified" if Path(STARS_CSV_FILE).is_file() else "missing",
            "reason": "Bundled star catalog; file existence only",
        }
    )
    path = root / EPHEMERIS_FILENAME
    checks.append(
        {
            "name": "ephemeris",
            "status": "unverified" if path.is_file() else "missing",
            "reason": "File existence only; coverage and contents are not inspected",
        }
    )
    return checks, None
