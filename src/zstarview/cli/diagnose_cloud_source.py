"""Diagnostic CLI for cloud source acquisition failures."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from ..clouddisc.diagnostics import FileDiagnosticSink, emit_diagnostic
from ..clouddisc.providers._goes_abi import load_cmi_with_area
from ..clouddisc.workers.cloud_source_worker import (
    WORKER_LOG_FILENAME,
    WORKER_RESULT_FILENAME,
)

DIAGNOSTIC_JSONL_FILENAME = "diagnostics.jsonl"
DIAGNOSTIC_LOG_FILENAME = "diagnose.log"
DIAGNOSTIC_SUMMARY_FILENAME = "diagnose.json"


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="zstarview-diagnose-cloud-source",
        description="Diagnose satellite cloud source acquisition failures.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where diagnostic logs, cache, and results are written.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow an existing non-empty output directory to be reused.",
    )
    parser.add_argument(
        "--source-file",
        type=Path,
        default=None,
        help="Diagnose an existing GOES source file without network access.",
    )
    parser.add_argument(
        "--satellite",
        default=None,
        help="Satellite name for --source-file diagnostics, such as G19 or G18.",
    )
    parser.add_argument(
        "--no-grid",
        action="store_true",
        help="Fetch and load source data without building the alt/az grid.",
    )
    _add_worker_compatible_arguments(parser)
    return parser


def _add_worker_compatible_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--request-id", type=int, default=0, help="Request identifier.")
    parser.add_argument("--lat", type=float, default=None, help="Observer latitude in degrees.")
    parser.add_argument("--lon", type=float, default=None, help="Observer longitude in degrees.")
    parser.add_argument(
        "--when-utc",
        type=str,
        default=None,
        help="Optional UTC timestamp for the source request (ISO-8601).",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Cloud cache directory. Defaults to OUTPUT_DIR/cache.",
    )
    parser.add_argument(
        "--sat-priority",
        nargs="+",
        default=("AUTO",),
        help="Satellite priority list used for selection.",
    )
    parser.add_argument("--bt-warm-k", type=float, default=310.0, help="Warm BT threshold.")
    parser.add_argument("--bt-cold-k", type=float, default=190.0, help="Cold BT threshold.")
    parser.add_argument("--alt-min-deg", type=float, default=1.0, help="Minimum satellite altitude.")
    parser.add_argument(
        "--search-back-minutes",
        type=int,
        default=120,
        help="How many minutes to search backwards for source data.",
    )
    parser.add_argument("--connect-timeout", type=float, default=5.0, help="Network connect timeout.")
    parser.add_argument("--read-timeout", type=float, default=30.0, help="Network read timeout.")
    parser.add_argument(
        "--cloud-shells-km",
        nargs="+",
        type=float,
        default=(6374.0, 6376.0, 6378.0),
        help="Cloud shell radii used for Himawari selection.",
    )


def _prepare_output_dir(path: Path, *, force: bool) -> Path:
    output_dir = path.expanduser().resolve()
    if output_dir.exists():
        if not output_dir.is_dir():
            raise SystemExit(f"Output path exists but is not a directory: {output_dir}")
        if any(output_dir.iterdir()) and not force:
            raise SystemExit(
                f"Output directory is not empty: {output_dir}\n"
                "Use --force to reuse it, or choose a new --output-dir."
            )
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def _build_worker_command(args: argparse.Namespace, *, output_dir: Path) -> list[str]:
    if args.lat is None or args.lon is None:
        raise SystemExit("--lat and --lon are required unless --source-file is used")
    cache_dir = (args.cache_dir or output_dir / "cache").expanduser().resolve()
    work_dir = output_dir / "worker"
    work_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        "-m",
        "zstarview.clouddisc.workers.cloud_source_worker",
        "--work-dir",
        str(work_dir),
        "--request-id",
        str(int(args.request_id)),
        "--lat",
        str(float(args.lat)),
        "--lon",
        str(float(args.lon)),
        "--cache-dir",
        str(cache_dir),
        "--sat-priority",
        *[str(v) for v in args.sat_priority],
        "--bt-warm-k",
        str(float(args.bt_warm_k)),
        "--bt-cold-k",
        str(float(args.bt_cold_k)),
        "--alt-min-deg",
        str(float(args.alt_min_deg)),
        "--search-back-minutes",
        str(int(args.search_back_minutes)),
        "--connect-timeout",
        str(float(args.connect_timeout)),
        "--read-timeout",
        str(float(args.read_timeout)),
        "--cloud-shells-km",
        *[str(float(v)) for v in args.cloud_shells_km],
        "--diagnostic-jsonl",
        str(output_dir / DIAGNOSTIC_JSONL_FILENAME),
        "--diagnostic-log",
        str(output_dir / DIAGNOSTIC_LOG_FILENAME),
    ]
    if args.when_utc:
        cmd.extend(["--when-utc", str(args.when_utc)])
    if args.no_grid:
        cmd.append("--skip-altaz-grid")
    return cmd


def _run_worker_diagnostic(args: argparse.Namespace, *, output_dir: Path) -> int:
    cmd = _build_worker_command(args, output_dir=output_dir)
    print("Running cloud source worker diagnostic:")
    print(" ".join(cmd))
    env = dict(os.environ)
    env["ZSTARVIEW_DISABLE_FILE_LOGGING"] = "1"
    process = subprocess.run(cmd, cwd=str(output_dir), check=False, env=env)
    result_path = output_dir / "worker" / WORKER_RESULT_FILENAME
    worker_log_path = output_dir / "worker" / WORKER_LOG_FILENAME
    summary = _build_summary(
        mode="worker",
        output_dir=output_dir,
        returncode=process.returncode,
        result_path=result_path,
        worker_log_path=worker_log_path,
    )
    _write_summary(output_dir, summary)
    _print_summary(summary)
    return _exit_code_for_summary(summary)


def _run_source_file_diagnostic(args: argparse.Namespace, *, output_dir: Path) -> int:
    source_file = Path(args.source_file).expanduser()
    jsonl_path = output_dir / DIAGNOSTIC_JSONL_FILENAME
    log_path = output_dir / DIAGNOSTIC_LOG_FILENAME
    sink = FileDiagnosticSink(jsonl_path, log_path)
    satellite = _infer_satellite(args.satellite, source_file)
    emit_diagnostic(
        sink,
        "resolve_source",
        "ok",
        "Using existing source file",
        satellite=satellite,
        path=source_file,
    )
    status = "succeeded"
    error_message = None
    try:
        if satellite not in {"G19", "G18"}:
            raise RuntimeError(
                "--source-file diagnostics currently support GOES G19/G18 files"
            )
        da = load_cmi_with_area(source_file, diagnostic_sink=sink)
        emit_diagnostic(
            sink,
            "build_altaz_grid",
            "info",
            "Alt/az grid skipped for source-file diagnostic",
            shape=tuple(int(v) for v in da.shape),
            dtype=str(da.dtype),
        )
    except Exception as exc:
        status = "failed"
        error_message = str(exc)
        emit_diagnostic(
            sink,
            "open_source_file",
            "failed",
            "Source file diagnostic failed",
            error_type=type(exc).__name__,
            error_message=str(exc),
        )
    summary = {
        "mode": "source-file",
        "status": status,
        "output_dir": str(output_dir),
        "source_file": str(source_file),
        "satellite": satellite,
        "diagnostic_jsonl": str(jsonl_path),
        "diagnostic_log": str(log_path),
        "error_message": error_message,
        "last_event": _load_last_event(jsonl_path),
    }
    _write_summary(output_dir, summary)
    _print_summary(summary)
    return _exit_code_for_summary(summary)


def _infer_satellite(satellite: str | None, source_file: Path) -> str | None:
    if satellite:
        return satellite.upper()
    name = source_file.name.upper()
    for candidate in ("G19", "G18", "HIMAWARI"):
        if candidate in name:
            return candidate
    return None


def _build_summary(
    *,
    mode: str,
    output_dir: Path,
    returncode: int,
    result_path: Path,
    worker_log_path: Path,
) -> dict[str, Any]:
    manifest = _load_json(result_path)
    jsonl_path = output_dir / DIAGNOSTIC_JSONL_FILENAME
    status = str(manifest.get("status", "failed")) if manifest else "failed"
    if returncode != 0 and status == "succeeded":
        status = "failed"
    return {
        "mode": mode,
        "status": status,
        "returncode": returncode,
        "output_dir": str(output_dir),
        "result_path": str(result_path),
        "worker_log": str(worker_log_path),
        "diagnostic_jsonl": str(jsonl_path),
        "diagnostic_log": str(output_dir / DIAGNOSTIC_LOG_FILENAME),
        "manifest": manifest,
        "last_event": _load_last_event(jsonl_path),
    }


def _load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as fp:
        payload = json.load(fp)
    return payload if isinstance(payload, dict) else {}


def _load_last_event(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    last = ""
    with path.open("r", encoding="utf-8") as fp:
        for line in fp:
            if line.strip():
                last = line
    if not last:
        return {}
    payload = json.loads(last)
    return payload if isinstance(payload, dict) else {}


def _write_summary(output_dir: Path, summary: dict[str, Any]) -> None:
    payload = dict(summary)
    payload["finished_at_utc"] = dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")
    path = output_dir / DIAGNOSTIC_SUMMARY_FILENAME
    with path.open("w", encoding="utf-8") as fp:
        json.dump(payload, fp, ensure_ascii=True, indent=2, sort_keys=True)
        fp.write("\n")


def _print_summary(summary: dict[str, Any]) -> None:
    status = str(summary.get("status", "failed"))
    print(f"Diagnostic status: {status}")
    print(f"Output directory: {summary.get('output_dir')}")
    last_event = summary.get("last_event")
    if isinstance(last_event, dict) and last_event:
        print(
            "Last diagnostic event: "
            f"{last_event.get('status')} {last_event.get('step')} - {last_event.get('message')}"
        )
    manifest = summary.get("manifest")
    if isinstance(manifest, dict) and manifest.get("error_message"):
        print(f"Worker error: {manifest.get('error_type')}: {manifest.get('error_message')}")
    elif summary.get("error_message"):
        print(f"Error: {summary.get('error_message')}")
    print(f"Diagnostic log: {summary.get('diagnostic_log')}")
    if summary.get("worker_log"):
        print(f"Worker log: {summary.get('worker_log')}")


def _exit_code_for_summary(summary: dict[str, Any]) -> int:
    if summary.get("status") == "succeeded":
        return 0
    last_event = summary.get("last_event")
    if not isinstance(last_event, dict):
        return 10
    step = str(last_event.get("step", ""))
    return {
        "list_s3_prefix": 2,
        "select_product": 3,
        "download_source": 4,
        "open_source_file": 5,
        "load_brightness_temperature": 5,
        "build_projection_area": 6,
        "build_altaz_grid": 7,
    }.get(step, 10)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)
    output_dir = _prepare_output_dir(args.output_dir, force=bool(args.force))
    args.cache_dir = (args.cache_dir.expanduser().resolve() if args.cache_dir is not None else None)
    if args.source_file is not None:
        args.source_file = Path(args.source_file).expanduser().resolve()
    if args.source_file is not None:
        return _run_source_file_diagnostic(args, output_dir=output_dir)
    return _run_worker_diagnostic(args, output_dir=output_dir)


if __name__ == "__main__":
    raise SystemExit(main())
