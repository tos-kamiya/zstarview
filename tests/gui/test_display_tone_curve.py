from __future__ import annotations

import numpy as np
import pytest
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QImage
from PySide6.QtTest import QTest

from zstarview.cli.args import parse_args
from zstarview.gui.display_tone_curve import (
    CALIBRATION_BLACK_VALUES,
    CALIBRATION_WHITE_VALUES,
    DisplayToneCalibrationDialog,
    DisplayToneCalibrationPattern,
    apply_display_tone_curve,
    build_display_tone_curve_lut,
)
from zstarview.gui.window_render_cache import SkyWindowRenderCacheMixin


def test_cli_display_tone_curve_default_and_pair() -> None:
    assert parse_args([]).display_tone_curve is None
    assert parse_args(["--display-tone-curve", "12,247"]).display_tone_curve == (
        12,
        247,
    )
    assert parse_args(["--display-tone-curve", "off"]).display_tone_curve is None


@pytest.mark.parametrize("value", ["12", "x,247", "-1,247", "12,12", "12,256"])
def test_cli_display_tone_curve_rejects_invalid_values(value: str) -> None:
    with pytest.raises(SystemExit):
        parse_args(["--display-tone-curve", value])


def test_display_tone_curve_lut_is_monotone_and_preserves_endpoints() -> None:
    lut = build_display_tone_curve_lut(12, 247)
    assert lut.dtype == np.uint8
    assert lut.shape == (256,)
    assert int(lut[0]) == 0
    assert int(lut[255]) == 255
    assert np.all(np.diff(lut.astype(np.int16)) >= 0)
    assert int(lut[1]) >= 12
    assert int(lut[254]) <= 247
    assert int(lut[64]) == 64
    assert int(lut[191]) == 191


@pytest.mark.parametrize("curve", [(0, 1), (64, 191), (200, 220), (254, 255)])
def test_display_tone_curve_lut_is_monotone_for_valid_extremes(
    curve: tuple[int, int],
) -> None:
    lut = build_display_tone_curve_lut(*curve)
    assert int(lut[0]) == 0
    assert int(lut[255]) == 255
    assert np.all(np.diff(lut.astype(np.int16)) >= 0)


def test_apply_display_tone_curve_preserves_alpha_and_endpoints() -> None:
    image = QImage(3, 1, QImage.Format.Format_RGBA8888)
    image.setPixelColor(0, 0, QColor(0, 0, 0, 17))
    image.setPixelColor(1, 0, QColor(1, 254, 100, 91))
    image.setPixelColor(2, 0, QColor(255, 255, 255, 203))

    corrected = apply_display_tone_curve(image, (12, 247))

    assert image.pixelColor(1, 0).getRgb() == (1, 254, 100, 91)
    assert corrected.pixelColor(0, 0).getRgb() == (0, 0, 0, 17)
    middle = corrected.pixelColor(1, 0)
    assert middle.red() >= 12
    assert middle.green() <= 247
    assert middle.alpha() == 91
    assert corrected.pixelColor(2, 0).getRgb() == (255, 255, 255, 203)


def test_display_frame_cache_reuses_conversion_and_off_uses_source() -> None:
    owner = SkyWindowRenderCacheMixin()
    owner.display_tone_curve = (12, 247)
    owner._display_frame_cache_key = None
    owner._display_frame_cache_image = None
    source = QImage(2, 2, QImage.Format.Format_ARGB32_Premultiplied)
    source.fill(QColor(1, 1, 1))

    first = owner._display_frame_image(source)
    second = owner._display_frame_image(source)

    assert first is second
    assert first is not source
    owner.display_tone_curve = None
    assert owner._display_frame_image(source) is source
    assert owner._display_frame_cache_key is None
    assert owner._display_frame_cache_image is None


def test_calibration_dialog_defaults_to_compact_window() -> None:
    dialog = DisplayToneCalibrationDialog(show_copy=True)

    assert dialog.size().width() == 800
    assert dialog.size().height() == 600
    assert not bool(dialog.windowState() & Qt.WindowState.WindowFullScreen)


def test_calibration_pattern_click_selects_numbered_patches() -> None:
    pattern = DisplayToneCalibrationPattern((12, 247))
    pattern.resize(800, 500)
    pattern.show()

    black_rects = pattern._patch_rects(
        pattern.rect().adjusted(0, 0, 0, -pattern.height() // 2),
        CALIBRATION_BLACK_VALUES,
    )
    assert black_rects[0].top() == 90
    assert black_rects[1].left() - black_rects[0].right() - 1 == 8
    QTest.mouseClick(pattern, Qt.MouseButton.LeftButton, pos=black_rects[5].center())
    assert pattern.black == 10
    assert pattern._last_selected == "black"

    white_section = pattern.rect().adjusted(0, pattern.height() // 2, 0, 0)
    white_rects = pattern._patch_rects(white_section, CALIBRATION_WHITE_VALUES)
    QTest.mouseClick(pattern, Qt.MouseButton.LeftButton, pos=white_rects[-3].center())
    assert pattern.white == 251
    assert pattern._last_selected == "white"
