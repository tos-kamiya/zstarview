import numpy as np

from zstarview.clouddisc.render.grayscale import convert_bt_to_rgba_image


def test_convert_bt_to_rgba_image_encodes_cloudiness_in_alpha() -> None:
    bt = np.array(
        [
            [310.0, 250.0, 190.0],
            [310.0, 220.0, 190.0],
        ],
        dtype=np.float32,
    )
    mask_inside = np.array(
        [
            [True, True, True],
            [False, True, True],
        ],
        dtype=bool,
    )
    rgba = convert_bt_to_rgba_image(bt, mask_inside, bt_warm=310.0, bt_cold=190.0)
    rgb = rgba[..., :3]
    a = rgba[..., 3]

    # RGB is white only inside visible region.
    assert tuple(int(v) for v in rgb[0, 0]) == (255, 255, 255)
    assert tuple(int(v) for v in rgb[1, 0]) == (0, 0, 0)

    # Alpha channel carries cloud amount.
    assert int(a[0, 0]) == 0
    assert int(a[0, 2]) == 255
    assert 0 < int(a[0, 1]) < 255
    assert int(a[1, 0]) == 0


def test_convert_bt_to_rgba_image_suppresses_very_low_cloud_amount() -> None:
    bt = np.array([[307.6, 304.0, 300.4, 250.0]], dtype=np.float32)
    mask_inside = np.array([[True, True, True, True]], dtype=bool)

    a = convert_bt_to_rgba_image(bt, mask_inside, bt_warm=310.0, bt_cold=190.0)[..., 3]

    # weight ~= 0.02 should be cut to zero.
    assert int(a[0, 0]) == 0
    # Around knee region, alpha should ramp gently and stay below dense cloud.
    assert 0 < int(a[0, 1]) <= int(a[0, 2])
    assert int(a[0, 2]) < int(a[0, 3])
