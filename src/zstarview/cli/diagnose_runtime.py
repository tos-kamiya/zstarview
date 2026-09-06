"""Print runtime information needed for crash reports."""

from __future__ import annotations

import argparse
import sys

from ..runtime_diagnostics import collect_runtime_diagnostics, format_runtime_diagnostics


def _print_ascii(value: object) -> None:
    print(str(value).encode("ascii", "backslashreplace").decode("ascii"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="zstarview-diagnose-runtime",
        description="Print zstarview runtime and dependency diagnostics.",
    )
    parser.add_argument("--json", action="store_true", help="Print formatted JSON.")
    args = parser.parse_args(argv)
    diagnostics = collect_runtime_diagnostics()
    if args.json:
        print(format_runtime_diagnostics(diagnostics))
    else:
        for key, value in diagnostics.items():
            if isinstance(value, dict):
                for child_key, child_value in value.items():
                    _print_ascii(f"{key}.{child_key}={child_value}")
            else:
                _print_ascii(f"{key}={value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
