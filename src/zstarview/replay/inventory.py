"""Read-only metadata inventory; layer payloads are never opened."""

from __future__ import annotations

import hashlib
import json
import math
import os
from collections import Counter
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

from ..aircraft.opensky import build_observer_bbox
from ..satellite_constants import SATELLITE_GROUP_VALIDITY_SECONDS
from ..satellites.cache import (
    _cache_has_all_expected_targets,
    satellite_cache_scope_key,
)
from ..satellites.cache import (
    _cached_set_from_payload as decode_satellites,
)
from .models import UTC, HistoryRecord, parse_utc

LAYERS = (
    "cloud-native",
    "cloud-geo",
    "aircraft",
    "satellite-iss",
    "satellite-horizons",
)


class InvalidRecord(ValueError):
    pass


class IncompatibleRecord(ValueError):
    pass


class ScopeMismatch(ValueError):
    pass


def json_paths(root: Path) -> list[Path]:
    """Missing roots are empty; permission and traversal errors are fatal."""
    try:
        root.stat()
        if not root.is_dir():
            raise NotADirectoryError(str(root))
    except FileNotFoundError:
        return []

    def fail(error: OSError) -> None:
        raise error

    result = []
    for directory, _, files in os.walk(root, onerror=fail, followlinks=False):
        for name in sorted(files):
            if name.endswith(".json") and not name.startswith("."):
                result.append(Path(directory) / name)
    return sorted(result)


def contained(root: Path, path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise InvalidRecord("Path escapes cache root")
    return resolved


def check_scope(
    lat: float, lon: float, expected_lat: float, expected_lon: float
) -> None:
    if not all(math.isfinite(v) for v in (lat, lon)):
        raise InvalidRecord("Non-finite observer coordinates")
    if abs(lat - expected_lat) > 0.00005 or abs(lon - expected_lon) > 0.00005:
        raise ScopeMismatch("Observer coordinates differ")


def check_asset_reference(root: Path, relative: str, digest: object) -> None:
    """Check reference metadata and file existence without opening the asset."""
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(c not in "0123456789abcdefABCDEF" for c in digest)
    ):
        raise InvalidRecord("Missing or invalid SHA-256 metadata")
    path = contained(root, root / relative)
    if not path.is_file():
        raise InvalidRecord("Referenced asset missing")


class Inventory:
    def __init__(
        self, root: Path, lat: float, lon: float, height: float | None
    ) -> None:
        self.root, self.lat, self.lon, self.height = root, lat, lon, height
        self.records: dict[str, list[HistoryRecord]] = {layer: [] for layer in LAYERS}
        self.diagnostics: dict[str, Counter[str]] = {
            layer: Counter() for layer in LAYERS
        }
        self._seen: set[tuple[str, str]] = set()

    def scan(self) -> None:
        for path in json_paths(self.root / "source_history"):
            self._read("cloud-native", path, self._native)
        for path in json_paths(self.root / "geosatellite" / "raw"):
            if path.name in {"infrared.json", "latest.json"}:
                self._read("cloud-geo", path, self._geo)
        for path in json_paths(self.root / "aircraft" / "opensky"):
            if "rate_limit" not in path.name and "rate-limit" not in path.name:
                self._read("aircraft", path, self._aircraft)
        satellite_root = self.root / "satellites" / "celestrak"
        for path in json_paths(satellite_root):
            for group in ("iss", "horizons"):
                if path.name == f"{group}.json" or path.name.startswith(f"{group}-"):
                    self._read(f"satellite-{group}", path, self._satellite)

    def _read(self, layer: str, path: Path, decoder: Any) -> None:
        try:
            raw = contained(self.root, path).read_bytes()
            payload = json.loads(raw)
            if not isinstance(payload, dict):
                raise InvalidRecord("Expected JSON object")
            record = decoder(path, payload)
            if record.layer != layer:
                raise InvalidRecord("Record group differs from filename")
            parts = path.relative_to(self.root).parts
            marker = "source_history" if "source_history" in parts else "history"
            if marker in parts:
                directory_time = (
                    record.source_time if layer == "cloud-geo" else record.fetched_at
                )
                if list(parts[parts.index(marker) + 1 : -1]) != directory_time.strftime(
                    "%Y/%m/%d/%H/%M"
                ).split("/"):
                    raise InvalidRecord("History directory time differs from record")
            identity = (
                f"{record.source_time.isoformat()}:{record.fetched_at.isoformat()}:{record.digest}"
                if record.digest
                else json.dumps(payload, sort_keys=True, separators=(",", ":"))
            )
            key = (layer, hashlib.sha256(identity.encode()).hexdigest())
            if key in self._seen:
                return
            self._seen.add(key)
            self.records[layer].append(record)
        except PermissionError:
            raise
        except ScopeMismatch:
            self.diagnostics[layer]["scope-mismatch"] += 1
        except IncompatibleRecord:
            self.diagnostics[layer]["incompatible"] += 1
        except (
            ValueError,
            TypeError,
            KeyError,
            AttributeError,
            EOFError,
            OverflowError,
            OSError,
        ):
            self.diagnostics[layer]["corrupt-or-missing"] += 1

    def _native(self, path: Path, payload: dict[str, Any]) -> HistoryRecord:
        if payload.get("schema_version") != 1:
            raise IncompatibleRecord("Unsupported cloud schema")
        check_scope(
            float(payload["observer_latitude_deg"]),
            float(payload["observer_longitude_deg"]),
            self.lat,
            self.lon,
        )
        source_time = parse_utc(payload["source_time_utc"])
        fetched = parse_utc(payload["fetched_at_utc"])
        if source_time > fetched:
            raise InvalidRecord("Cloud observation is after acquisition")
        check_asset_reference(
            self.root, payload["artifact_path"], payload["artifact_sha256"]
        )
        return HistoryRecord(
            "cloud-native",
            path,
            source_time,
            fetched,
            fetched + timedelta(minutes=30),
            digest=payload["artifact_sha256"],
        )

    def _geo(self, path: Path, payload: dict[str, Any]) -> HistoryRecord:
        if payload.get("schema_version") != 1:
            raise IncompatibleRecord("Unsupported Geo-satellite schema")
        if not (32 <= self.lat <= 73 and -15 <= self.lon <= 35):
            raise ScopeMismatch("Outside Europe coverage")
        if (payload.get("kind"), payload.get("area"), payload.get("size")) != (
            "infrared",
            "europe",
            "normal",
        ):
            raise ScopeMismatch("Different image kind, area or size")
        geo_root = self.root / "geosatellite"
        image_path = payload.get(
            "image_path", str(path.with_suffix(".png").relative_to(geo_root))
        )
        resolved = contained(geo_root, geo_root / image_path)
        if resolved != path.with_suffix(".png").resolve():
            raise InvalidRecord("Image reference differs from sidecar")
        check_asset_reference(geo_root, image_path, payload["sha256"])
        source = parse_utc(payload["captured_at_utc"])
        fetched = parse_utc(payload["fetched_at_utc"])
        if source > fetched:
            raise InvalidRecord("Cloud observation is after acquisition")
        return HistoryRecord(
            "cloud-geo",
            path,
            source,
            fetched,
            fetched + timedelta(minutes=60),
            digest=payload["sha256"],
        )

    def _aircraft(self, path: Path, payload: dict[str, Any]) -> HistoryRecord:
        if any(not isinstance(item, dict) for item in payload["snapshots"]):
            raise InvalidRecord("Invalid aircraft snapshot")
        from ..aircraft.opensky import AircraftBoundingBox

        expected = build_observer_bbox(observer_lat=self.lat, observer_lon=self.lon)
        if AircraftBoundingBox(**payload["bbox"]) != expected:
            raise ScopeMismatch("Different aircraft bbox")
        fetched = parse_utc(payload["fetched_at_utc"])
        contacts = tuple(
            datetime.fromtimestamp(int(item["last_contact_unix"]), UTC)
            for item in payload["snapshots"]
            if item.get("last_contact_unix") is not None
        )
        source = max(contacts, default=fetched)
        return HistoryRecord(
            "aircraft",
            path,
            source,
            fetched,
            source + timedelta(minutes=10),
            contacts,
            not payload["snapshots"],
        )

    def _satellite(self, path: Path, payload: dict[str, Any]) -> HistoryRecord:
        group = str(payload["group_key"])
        if group not in {"iss", "horizons"}:
            raise IncompatibleRecord("Unknown satellite group")
        expected_scope = satellite_cache_scope_key(
            observer_lat=self.lat,
            observer_lon=self.lon,
            observer_height_m=self.height or 0,
        )
        relative = path.relative_to(self.root / "satellites" / "celestrak")
        if group == "horizons" and (
            self.height is None or relative.parts[0] != expected_scope
        ):
            raise ScopeMismatch("Different or unknown satellite observer scope")
        if group == "iss" and relative.parts[0].startswith("earth_"):
            raise ScopeMismatch("Unexpected ISS scope")
        if not payload.get("fetched_at_utc"):
            raise InvalidRecord("Missing acquisition time")
        cached = decode_satellites(payload, group_key=group)
        if not cached.records or not _cache_has_all_expected_targets(group, cached):
            raise IncompatibleRecord("Incomplete satellite group")
        epochs = [_element_epoch(r["EPOCH"]) for r in cached.records if r.get("EPOCH")]
        start = max(epochs, default=cached.element_epoch_utc)
        end = min(epochs, default=cached.element_epoch_utc) + timedelta(
            seconds=SATELLITE_GROUP_VALIDITY_SECONDS[group]
        )
        return HistoryRecord(
            f"satellite-{group}", path, start, cached.fetched_at_utc, end
        )


def _element_epoch(value: object) -> datetime:
    """OMM epochs without a suffix are UTC by the satellite cache convention."""
    parsed = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    return (
        parsed.replace(tzinfo=UTC) if parsed.tzinfo is None else parsed.astimezone(UTC)
    )
