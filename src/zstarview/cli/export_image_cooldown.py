"""Best-effort startup cooldown for the one-shot image exporter."""

from __future__ import annotations

import json
import logging
import math
import os
import tempfile
import time
import uuid
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

logger = logging.getLogger(__name__)
EXPORT_COOLDOWN_SECONDS = 60.0
EXPORT_COOLDOWN_FILE = "export-image-start.json"


def _read_record(path: Path) -> tuple[float, str] | None:
    try:
        contents = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return None
    except UnicodeDecodeError:
        logger.warning("Invalid export cooldown record; replacing it.")
        return None
    try:
        data = json.loads(contents)
        timestamp = data["started_at"]
        token = data["token"]
        if isinstance(timestamp, bool) or not isinstance(timestamp, (int, float)):
            raise ValueError
        if not math.isfinite(timestamp) or not isinstance(token, str) or not token:
            raise ValueError
        return float(timestamp), token
    except (ValueError, TypeError, KeyError):
        logger.warning("Invalid export cooldown record; replacing it.")
        return None


def _write_record(path: Path, timestamp: float, token: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent, delete=False
        ) as temp:
            temp_path = Path(temp.name)
            json.dump({"started_at": timestamp, "token": token}, temp)
        os.replace(temp_path, path)
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)


@contextmanager
def export_cooldown(cache_path: Path, *, bypass: bool = False) -> Iterator[None]:
    """Throttle starts without claiming to provide mutual exclusion."""
    if bypass:
        yield
        return

    path = cache_path / EXPORT_COOLDOWN_FILE
    try:
        record = _read_record(path)
        now = time.time()
        if record is not None:
            age = now - record[0]
            if 0 <= age < EXPORT_COOLDOWN_SECONDS:
                remaining = math.ceil(EXPORT_COOLDOWN_SECONDS - age)
                raise SystemExit(
                    f"Image export started recently; retry in {remaining} seconds. "
                    "For debugging, use --debug-ignore-export-cooldown."
                )
            if age < 0:
                logger.warning(
                    "Export cooldown timestamp is in the future; replacing it."
                )
        token = uuid.uuid4().hex
        _write_record(path, now, token)
    except OSError as exc:
        raise SystemExit(
            "Cannot access the export cooldown record. "
            "For debugging, use --debug-ignore-export-cooldown."
        ) from exc

    try:
        yield
    finally:
        try:
            current = _read_record(path)
            if current is not None and current[1] != token:
                logger.warning(
                    "The export startup record changed; another image export "
                    "may have overlapped this run."
                )
        except OSError:
            logger.warning("Could not check the export startup record for overlap.")
