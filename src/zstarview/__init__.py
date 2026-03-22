# SPDX-FileCopyrightText: 2025-present Toshihiro Kamiya <kamiya@mbj.nifty.com>
#
# SPDX-License-Identifier: MIT


def main() -> None:
    from .gui.viewer import main as _main

    _main()


def debug_main() -> None:
    from .gui.viewer import main as _main

    _main()


def make_desktop_file_main() -> None:
    from .cli.make_desktop_file import main as _main

    _main()


def export_image_main() -> None:
    from .cli.export_image import main as _main

    _main()
