from __future__ import annotations

import numpy as np

from zstarview.render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba


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
