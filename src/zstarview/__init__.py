# SPDX-FileCopyrightText: 2025-present Toshihiro Kamiya <kamiya@mbj.nifty.com>
#
# SPDX-License-Identifier: MIT

from .offline import install_from_env as _install_offline_guard_from_env

# Worker subprocesses (``python -m zstarview...``) inherit offline mode
# through the environment; apply the guards as soon as the package loads.
_install_offline_guard_from_env()


def main() -> None:
    from .gui.viewer import main as _main

    _main()


def make_desktop_file_main() -> None:
    from .cli.make_desktop_file import main as _main

    _main()


def export_image_main() -> None:
    from .cli.export_image import main as _main

    _main()
