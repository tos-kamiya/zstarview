import os
import argparse
import os.path

_dir = os.path.dirname(os.path.abspath(__file__))

from .paths import APP_ID

APP_ICON_FILE = os.path.join(_dir, "data", "icon-256.png")
APP_COMMAND = "zstarview"
DESKTOP_FILE = f"{APP_ID}.desktop"

DESKTOP_TEMPLATE = """[Desktop Entry]
Type=Application
Version=1.0
Name=Zenith Star View
Comment=Interactive night sky viewer
Exec={exec_cmd}
Icon={icon_path}
Terminal=false
Categories=Education;Science;Astronomy;
StartupWMClass={wmclass}
"""


def main():
    parser = argparse.ArgumentParser(description=".desktop file generator for Zenith Star View")
    parser.add_argument(
        "--write",
        action="store_true",
        help=f"Write the desktop file to ~/.local/share/applications/{APP_ID} (default: output the file in the current directory)",
    )
    args = parser.parse_args()

    exec_cmd = APP_COMMAND

    out_dir = os.path.expanduser("~/.local/share/applications") if args.write else os.getcwd()
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, DESKTOP_FILE)

    content = DESKTOP_TEMPLATE.format(
        exec_cmd=exec_cmd,
        icon_path=APP_ICON_FILE,
        wmclass=APP_ID,
    )

    print(f"{out_path=}")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(content)
    os.chmod(out_path, 0o644)

    print(f"Wrote: {out_path}")

    if args.write:
        os.system(f'update-desktop-database "{out_dir}" >/dev/null 2>&1')


if __name__ == "__main__":
    main()
