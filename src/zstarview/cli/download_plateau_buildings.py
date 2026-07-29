"""Console entry point for preparing PLATEAU building cache data."""

from __future__ import annotations

from collections.abc import Sequence

from ..data.plateau_buildings import main as plateau_main


def main(argv: Sequence[str] | None = None) -> int:
    return plateau_main(argv)
