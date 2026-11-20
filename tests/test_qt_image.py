from __future__ import annotations

import numpy as np
from PySide6.QtGui import QImage

from zstarview.render.qt_image import (
    np_rgba_to_qimage,
    qimage_to_np_rgba,
    save_qimage_png,
)


def test_np_rgba_to_qimage_round_trip_for_contiguous_array() -> None:
    arr = np.array(
        [
            [[1, 2, 3, 4], [5, 6, 7, 8]],
            [[9, 10, 11, 12], [13, 14, 15, 16]],
        ],
        dtype=np.uint8,
    )

    qimg = np_rgba_to_qimage(arr)

    assert qimage_to_np_rgba(qimg).tolist() == arr.tolist()


def test_np_rgba_to_qimage_handles_non_contiguous_view() -> None:
    base = np.arange(4 * 4 * 4, dtype=np.uint8).reshape(4, 4, 4)
    view = base[::2, :, :]

    qimg = np_rgba_to_qimage(view)
    round_trip = qimage_to_np_rgba(qimg)

    assert view.flags.c_contiguous is False
    assert round_trip.tolist() == np.ascontiguousarray(view).tolist()


def test_save_qimage_png_writes_png_regardless_of_extension(tmp_path) -> None:
    # replay/sequence.py writes to ".<name>.png.tmp" before renaming.
    arr = np.zeros((3, 5, 4), dtype=np.uint8)
    arr[..., 0] = 200
    arr[..., 3] = 255
    path = tmp_path / ".frame.png.tmp"

    assert save_qimage_png(np_rgba_to_qimage(arr), path)

    assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    loaded = QImage(str(path), "PNG")
    assert not loaded.isNull()
    assert qimage_to_np_rgba(loaded).tolist() == arr.tolist()


def test_save_qimage_png_reports_failure(tmp_path) -> None:
    arr = np.zeros((2, 2, 4), dtype=np.uint8)

    assert not save_qimage_png(np_rgba_to_qimage(arr), tmp_path / "missing" / "a.png")
