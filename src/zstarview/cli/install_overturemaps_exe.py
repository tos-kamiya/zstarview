"""Stage a downloaded OvertureMaps Windows executable into the app cache."""

from __future__ import annotations

import argparse
import shutil
import sys
from collections.abc import Sequence
from pathlib import Path

from ..data.import_overture_buildings import staged_overturemaps_executable_path


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="zstarview-install-overturemaps-exe-cli",
        description="Copy a downloaded OvertureMaps Windows executable into the zstarview cache.",
    )
    parser.add_argument(
        "source_exe",
        type=Path,
        help="Path to the downloaded OvertureMaps Windows executable.",
    )
    return parser


def stage_overturemaps_executable(source_exe: Path) -> Path:
    source_path = source_exe.expanduser()
    if not source_path.exists():
        raise FileNotFoundError(f"Source executable does not exist: {source_path}")

    destination_path = staged_overturemaps_executable_path()
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_path, destination_path)
    return destination_path


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)
    try:
        destination = stage_overturemaps_executable(args.source_exe)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except PermissionError as exc:
        print(f"Failed to copy OvertureMaps executable: {exc}", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"Failed to copy OvertureMaps executable: {exc}", file=sys.stderr)
        return 1

    print(f"Copied OvertureMaps executable to {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
