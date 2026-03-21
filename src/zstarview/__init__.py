# SPDX-FileCopyrightText: 2025-present Toshihiro Kamiya <kamiya@mbj.nifty.com>
#
# SPDX-License-Identifier: MIT


def main() -> None:
    from .zstarview import main as _main

    _main()


def make_desktop_file_main() -> None:
    from .make_desktop_file import main as _main

    _main()


def export_image_main() -> None:
    from .export_image import main as _main

    _main()
