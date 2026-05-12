# SPDX-FileCopyrightText: 2025-present Toshihiro Kamiya <kamiya@mbj.nifty.com>
#
# SPDX-License-Identifier: MIT
"""Current GUI application entrypoint.

This module intentionally keeps a thin wrapper around the legacy viewer entrypoint
so the command-line surface can move to ``zstarview-gui`` without deleting the
older implementation file yet.
"""

from __future__ import annotations


def main() -> None:
    from .viewer import main as _main

    _main()


if __name__ == "__main__":
    main()
