# -*- coding: utf-8 -*-
"""One-shot subprocess worker for cloud source fetches."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import os
import pickle
import subprocess
import sys
import time
import traceback
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from .. import CloudDisc, CloudDiscConfig
from ..types import (
    CloudDiscError,
    CloudSourceData,
    DataNotFoundError,
    DownloadCancelledError,
    DownloadError,
    RenderError,
    TimeoutError,
    VisibilityError,
)
from .cloud_source import (
    CloudSourceFetchRequest,
    build_cloud_source_fetch_request,
    fetch_cloud_source,
)
from .constants import DEFAULT_CLOUD_SHELLS_KM

logger = logging.getLogger(__name__)

WORKER_RESULT_FILENAME = "result.json"
WORKER_ARTIFACT_FILENAME = "source.pkl"
WORKER_LOG_FILENAME = "worker.log"
DEFAULT_WORKER_TIMEOUT_S = 900.0
DEFAULT_POLL_INTERVAL_S = 0.25


@dataclass(frozen=True, slots=True)
class CloudSourceWorkerManifest:
    """JSON result written by the worker subprocess."""

    request_id: int
    status: str
    started_at_utc: str
    finished_at_utc: str
    work_dir: str
    artifact_path: str | None = None
    satellite: str | None = None
    product: str | None = None
    time_utc: str | None = None
    source_key: dict[str, object] | None = None
    src_paths: list[str] | None = None
    source_expected_count: int | None = None
    source_available_count: int | None = None
    source_completeness_ratio: float | None = None
    error_type: str | None = None
    error_message: str | None = None
    error_traceback: str | None = None


def _isoformat_utc(value: dt.datetime) -> str:
    value = value.astimezone(dt.timezone.utc)
    return value.isoformat().replace("+00:00", "Z")


def _parse_utc_datetime(value: str | None) -> dt.datetime | None:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    parsed = dt.datetime.fromisoformat(text)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed.astimezone(dt.timezone.utc)


def _write_text_atomic(path: Path, payload: str) -> None:
    tmp_path = path.with_name(path.name + ".tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tmp_path.open("w", encoding="utf-8") as fp:
        fp.write(payload)
    tmp_path.replace(path)


def _write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    _write_text_atomic(
        path, json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True)
    )


def _write_pickle_atomic(path: Path, payload: object) -> None:
    tmp_path = path.with_name(path.name + ".tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tmp_path.open("wb") as fp:
        pickle.dump(payload, fp, protocol=pickle.HIGHEST_PROTOCOL)
    tmp_path.replace(path)


def _encode_source_key(source_key) -> dict[str, object]:
    return {
        "satellite": str(source_key.satellite),
        "provider": None if source_key.provider is None else str(source_key.provider),
        "timeslot_utc": _isoformat_utc(source_key.timeslot_utc),
        "sat_priority": [str(v) for v in source_key.sat_priority],
    }


def _decode_error_type(error_type: str | None) -> type[CloudDiscError]:
    if error_type == "DownloadCancelledError":
        return DownloadCancelledError
    if error_type == "DataNotFoundError":
        return DataNotFoundError
    if error_type == "DownloadError":
        return DownloadError
    if error_type == "TimeoutError":
        return TimeoutError
    if error_type == "RenderError":
        return RenderError
    if error_type == "VisibilityError":
        return VisibilityError
    return CloudDiscError


def _build_manifest_payload(
    *,
    request_id: int,
    status: str,
    started_at_utc: dt.datetime,
    finished_at_utc: dt.datetime,
    work_dir: Path,
    artifact_path: Path | None = None,
    source: CloudSourceData | None = None,
    error: BaseException | None = None,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "request_id": int(request_id),
        "status": str(status),
        "started_at_utc": _isoformat_utc(started_at_utc),
        "finished_at_utc": _isoformat_utc(finished_at_utc),
        "work_dir": str(work_dir),
        "artifact_path": None if artifact_path is None else str(artifact_path),
    }
    if source is not None:
        payload.update(
            {
                "satellite": str(getattr(source, "satellite", "")),
                "product": str(getattr(source, "product", "")),
                "time_utc": _isoformat_utc(getattr(source, "time_utc")),
                "source_key": _encode_source_key(source.source_key),
                "src_paths": [
                    str(Path(path)) for path in getattr(source, "src_paths", [])
                ],
                "source_expected_count": getattr(source, "source_expected_count", None),
                "source_available_count": getattr(
                    source, "source_available_count", None
                ),
                "source_completeness_ratio": getattr(
                    source, "source_completeness_ratio", None
                ),
            }
        )
    if error is not None:
        payload.update(
            {
                "error_type": type(error).__name__,
                "error_message": str(error),
                "error_traceback": traceback.format_exc(),
            }
        )
    return payload


def _build_cloud_disc_from_args(args: argparse.Namespace) -> CloudDisc:
    return CloudDisc(
        CloudDiscConfig(
            cache_dir=args.cache_dir,
            sat_priority=tuple(args.sat_priority),
            bt_warm_k=float(args.bt_warm_k),
            bt_cold_k=float(args.bt_cold_k),
            alt_min_deg=float(args.alt_min_deg),
            search_back_minutes=int(args.search_back_minutes),
            connect_timeout=float(args.connect_timeout),
            read_timeout=float(args.read_timeout),
        )
    )


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="zstarview.clouddisc.workers.cloud_source_worker",
        description="Run one cloud-source fetch request in a subprocess.",
    )
    parser.add_argument(
        "--work-dir", type=Path, required=True, help="Worker output directory."
    )
    parser.add_argument(
        "--request-id", type=int, required=True, help="Parent request identifier."
    )
    parser.add_argument(
        "--lat", type=float, required=True, help="Observer latitude in degrees."
    )
    parser.add_argument(
        "--lon", type=float, required=True, help="Observer longitude in degrees."
    )
    parser.add_argument(
        "--when-utc",
        type=str,
        default=None,
        help="Optional UTC timestamp for the source request (ISO-8601).",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        required=False,
        default=None,
        help="Cloud cache directory.",
    )
    parser.add_argument(
        "--sat-priority",
        nargs="+",
        default=("AUTO",),
        help="Satellite priority list used for selection.",
    )
    parser.add_argument(
        "--bt-warm-k", type=float, default=310.0, help="Warm BT threshold."
    )
    parser.add_argument(
        "--bt-cold-k", type=float, default=190.0, help="Cold BT threshold."
    )
    parser.add_argument(
        "--alt-min-deg", type=float, default=0.0, help="Minimum satellite altitude."
    )
    parser.add_argument(
        "--search-back-minutes",
        type=int,
        default=120,
        help="How many minutes to search backwards for source data.",
    )
    parser.add_argument(
        "--connect-timeout", type=float, default=5.0, help="Network connect timeout."
    )
    parser.add_argument(
        "--read-timeout", type=float, default=30.0, help="Network read timeout."
    )
    parser.add_argument(
        "--cloud-shells-km",
        nargs="+",
        type=float,
        default=DEFAULT_CLOUD_SHELLS_KM,
        help="Cloud shell radii used for Himawari selection.",
    )
    return parser


def _make_request_from_args(args: argparse.Namespace) -> CloudSourceFetchRequest:
    return build_cloud_source_fetch_request(
        lat=float(args.lat),
        lon=float(args.lon),
        when_utc=_parse_utc_datetime(args.when_utc),
        cloud_shells_km=tuple(float(v) for v in args.cloud_shells_km),
    )


def _build_worker_command(
    *,
    clouddisc: CloudDisc,
    request_id: int,
    lat: float,
    lon: float,
    work_dir: Path,
    when_utc: dt.datetime | None = None,
    cloud_shells_km: Sequence[float] = DEFAULT_CLOUD_SHELLS_KM,
) -> list[str]:
    cfg = clouddisc.cfg
    cmd = [
        sys.executable,
        "-m",
        "zstarview.clouddisc.workers.cloud_source_worker",
        "--work-dir",
        str(work_dir),
        "--request-id",
        str(int(request_id)),
        "--lat",
        str(float(lat)),
        "--lon",
        str(float(lon)),
        "--cache-dir",
        str(cfg.cache_root()),
        "--sat-priority",
        *[str(v) for v in cfg.sat_priority],
        "--bt-warm-k",
        str(float(cfg.bt_warm_k)),
        "--bt-cold-k",
        str(float(cfg.bt_cold_k)),
        "--alt-min-deg",
        str(float(cfg.alt_min_deg)),
        "--search-back-minutes",
        str(int(cfg.search_back_minutes)),
        "--connect-timeout",
        str(float(cfg.connect_timeout)),
        "--read-timeout",
        str(float(cfg.read_timeout)),
        "--cloud-shells-km",
        *[str(float(v)) for v in cloud_shells_km],
    ]
    if when_utc is not None:
        cmd.extend(["--when-utc", _isoformat_utc(when_utc)])
    return cmd


def _make_default_worker_work_dir(
    *,
    cache_root: Path,
    request_id: int,
    parent_pid: int,
    parent_started_at_utc: dt.datetime,
) -> Path:
    session_stamp = (
        _isoformat_utc(parent_started_at_utc).replace(":", "").replace("-", "")
    )
    root = Path(cache_root) / "cloud-worker" / f"parent-{parent_pid}-{session_stamp}"
    request_dir = root / f"request-{int(request_id):08d}"
    request_dir.mkdir(parents=True, exist_ok=True)
    return request_dir


def _load_manifest(path: Path) -> CloudSourceWorkerManifest:
    with path.open("r", encoding="utf-8") as fp:
        payload = json.load(fp)
    return CloudSourceWorkerManifest(
        request_id=int(payload["request_id"]),
        status=str(payload["status"]),
        started_at_utc=str(payload["started_at_utc"]),
        finished_at_utc=str(payload["finished_at_utc"]),
        work_dir=str(payload["work_dir"]),
        artifact_path=payload.get("artifact_path"),
        satellite=payload.get("satellite"),
        product=payload.get("product"),
        time_utc=payload.get("time_utc"),
        source_key=payload.get("source_key"),
        src_paths=payload.get("src_paths"),
        source_expected_count=payload.get("source_expected_count"),
        source_available_count=payload.get("source_available_count"),
        source_completeness_ratio=payload.get("source_completeness_ratio"),
        error_type=payload.get("error_type"),
        error_message=payload.get("error_message"),
        error_traceback=payload.get("error_traceback"),
    )


def _load_artifact(path: Path) -> CloudSourceData:
    with path.open("rb") as fp:
        source = pickle.load(fp)
    if not isinstance(source, CloudSourceData):
        raise TypeError(f"Unexpected worker artifact type: {type(source)!r}")
    return source


def _cleanup_worker_work_dir(work_dir: Path) -> None:
    try:
        for child in work_dir.iterdir():
            if child.is_file() or child.is_symlink():
                child.unlink(missing_ok=True)
            elif child.is_dir():
                import shutil

                shutil.rmtree(child, ignore_errors=True)
        work_dir.rmdir()
    except OSError:
        logger.debug("Failed to clean worker work dir: %s", work_dir, exc_info=True)


def load_cloud_source_worker_result(result_path: Path) -> CloudSourceData:
    """Load and validate a worker result manifest + artifact pair."""
    manifest = _load_manifest(result_path)
    if manifest.status != "succeeded":
        error_cls = _decode_error_type(manifest.error_type)
        message = manifest.error_message or "cloud source worker failed"
        raise error_cls(message)

    artifact_path_value = manifest.artifact_path
    if not artifact_path_value:
        raise DownloadError(
            "cloud source worker result did not include an artifact path"
        )
    artifact_path = Path(artifact_path_value)
    source = _load_artifact(artifact_path)
    source.sampler = None
    return source


def run_cloud_source_worker_process(
    clouddisc: CloudDisc,
    request: CloudSourceFetchRequest,
    *,
    request_id: int,
    work_dir: Path | None = None,
    abort_event: threading.Event | None = None,
    timeout_s: float = DEFAULT_WORKER_TIMEOUT_S,
    poll_interval_s: float = DEFAULT_POLL_INTERVAL_S,
) -> CloudSourceData:
    """Run the cloud source worker as a one-shot subprocess."""
    if abort_event is not None and abort_event.is_set():
        raise DownloadCancelledError("cloud source worker cancelled")

    parent_pid = os.getpid()
    parent_started_at_utc = dt.datetime.now(dt.timezone.utc)
    if work_dir is None:
        work_dir = _make_default_worker_work_dir(
            cache_root=clouddisc.cfg.cache_root(),
            request_id=request_id,
            parent_pid=parent_pid,
            parent_started_at_utc=parent_started_at_utc,
        )
    else:
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

    cmd = _build_worker_command(
        clouddisc=clouddisc,
        request_id=request_id,
        lat=request.lat,
        lon=request.lon,
        work_dir=work_dir,
        when_utc=request.when_utc,
        cloud_shells_km=request.cloud_shells_km,
    )
    log_path = work_dir / WORKER_LOG_FILENAME
    result_path = work_dir / WORKER_RESULT_FILENAME
    logger.info("Launching cloud source worker: %s", " ".join(cmd))

    with log_path.open("ab") as log_fp:
        process = subprocess.Popen(
            cmd,
            cwd=str(work_dir),
            stdout=log_fp,
            stderr=log_fp,
        )
        deadline = time.monotonic() + max(0.0, float(timeout_s))
        terminated_by_abort = False
        try:
            while True:
                if abort_event is not None and abort_event.is_set():
                    terminated_by_abort = True
                    process.terminate()
                    try:
                        process.wait(timeout=5.0)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait(timeout=5.0)
                    raise DownloadCancelledError("cloud source worker cancelled")
                remaining = deadline - time.monotonic()
                if remaining <= 0.0:
                    process.terminate()
                    try:
                        process.wait(timeout=5.0)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait(timeout=5.0)
                    raise TimeoutError("cloud source worker timed out")
                try:
                    returncode = process.wait(timeout=min(poll_interval_s, remaining))
                except subprocess.TimeoutExpired:
                    continue
                break
        finally:
            if process.poll() is None and not terminated_by_abort:
                process.kill()

    if not result_path.exists():
        raise DownloadError(
            f"cloud source worker exited with code {returncode} without writing a result"
        )

    manifest = _load_manifest(result_path)
    if manifest.status != "succeeded":
        error_cls = _decode_error_type(manifest.error_type)
        message = manifest.error_message or "cloud source worker failed"
        raise error_cls(message)

    source = _load_artifact(Path(manifest.artifact_path or ""))
    source.sampler = None

    if returncode != 0:
        raise DownloadError(f"cloud source worker exited with code {returncode}")

    _cleanup_worker_work_dir(work_dir)
    return source


def _run_one_shot_worker(
    *,
    args: argparse.Namespace,
) -> int:
    started_at_utc = dt.datetime.now(dt.timezone.utc)
    work_dir = Path(args.work_dir)
    result_path = work_dir / WORKER_RESULT_FILENAME
    artifact_path = work_dir / WORKER_ARTIFACT_FILENAME
    cloud_disc = _build_cloud_disc_from_args(args)
    request = _make_request_from_args(args)

    try:
        logger.info(
            "Worker started (request_id=%s, lat=%.4f, lon=%.4f, work_dir=%s)",
            args.request_id,
            request.lat,
            request.lon,
            work_dir,
        )
        source = fetch_cloud_source(cloud_disc, request)
        source.sampler = None
        logger.info("Building alt/az grid in worker...")
        source.altaz_grid = cloud_disc.build_altaz_grid_from_source(
            source=source,
            lat=request.lat,
            lon=request.lon,
            cloud_shells_km=request.cloud_shells_km,
        )
        logger.info("Alt/az grid built.")
        _write_pickle_atomic(artifact_path, source)
        finished_at_utc = dt.datetime.now(dt.timezone.utc)
        _write_json_atomic(
            result_path,
            _build_manifest_payload(
                request_id=int(args.request_id),
                status="succeeded",
                started_at_utc=started_at_utc,
                finished_at_utc=finished_at_utc,
                work_dir=work_dir,
                artifact_path=artifact_path,
                source=source,
            ),
        )
        logger.info(
            "Worker finished successfully (request_id=%s, satellite=%s, product=%s)",
            args.request_id,
            getattr(source, "satellite", "?"),
            getattr(source, "product", "?"),
        )
        return 0
    except Exception as exc:
        finished_at_utc = dt.datetime.now(dt.timezone.utc)
        try:
            _write_json_atomic(
                result_path,
                _build_manifest_payload(
                    request_id=int(args.request_id),
                    status="failed",
                    started_at_utc=started_at_utc,
                    finished_at_utc=finished_at_utc,
                    work_dir=work_dir,
                    artifact_path=None,
                    error=exc,
                ),
            )
        except Exception:
            logger.error("Failed to write worker failure manifest", exc_info=True)
        logger.error("Cloud source worker failed: %s", exc, exc_info=True)
        return 1


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)
    from ...logging_utils import setup_root_logger

    setup_root_logger()
    args.when_utc = args.when_utc or None
    args.work_dir = Path(args.work_dir)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    return _run_one_shot_worker(args=args)


if __name__ == "__main__":
    raise SystemExit(main())
