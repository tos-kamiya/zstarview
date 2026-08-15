from __future__ import annotations

import json
from datetime import datetime, timezone

from PySide6.QtCore import QBuffer, QIODevice
from PySide6.QtGui import QColor, QImage

from zstarview.solar_hover import (
    apply_solar_black_to_alpha,
    closest_image_url,
    fetch_solar_hover_image,
    normalize_solar_hover_time,
    screenshot_url,
)


class _Response:
    def __init__(self, payload: bytes) -> None:
        self._payload = payload

    def __enter__(self) -> _Response:
        return self

    def __exit__(self, *args: object) -> None:
        return None

    def read(self) -> bytes:
        return self._payload


def _png_bytes() -> bytes:
    image = QImage(1024, 1024, QImage.Format.Format_RGB32)
    image.fill(0)
    buffer = QBuffer()
    buffer.open(QIODevice.OpenModeFlag.WriteOnly)
    assert image.save(buffer, "PNG")
    return bytes(buffer.data())


def test_solar_hover_time_uses_ten_minute_bucket() -> None:
    target = datetime(2026, 8, 15, 6, 29, 59, tzinfo=timezone.utc)
    assert normalize_solar_hover_time(target) == datetime(
        2026, 8, 15, 6, 20, tzinfo=timezone.utc
    )
    assert "sourceId=18" in closest_image_url(target)


def test_screenshot_url_is_centered_hmi_continuum() -> None:
    target = datetime(2026, 8, 15, 6, 30, tzinfo=timezone.utc)
    url = screenshot_url(target)
    assert "%5BSDO%2CHMI%2CHMI%2Ccontinuum%2C1%2C100%5D" in url
    assert "x0=0" in url
    assert "y0=0" in url
    assert "watermark=false" in url


def test_black_to_alpha_preserves_disc_and_reveals_outer_emission() -> None:
    image = QImage(21, 21, QImage.Format.Format_RGB32)
    image.fill(0)
    emission = QColor(128, 12, 0)
    image.setPixelColor(10, 10, emission)
    image.setPixelColor(17, 10, emission)
    image.setPixelColor(20, 10, emission)

    converted = apply_solar_black_to_alpha(
        image,
        5.0,
        feather_width_px=4.0,
    )

    disc = converted.pixelColor(10, 10)
    assert (disc.red(), disc.green(), disc.blue(), disc.alpha()) == (128, 12, 0, 255)

    transition = converted.pixelColor(17, 10)
    assert 128 < transition.alpha() < 255

    outside = converted.pixelColor(20, 10)
    assert outside.alpha() == 128
    assert outside.red() == 255
    assert 23 <= outside.green() <= 24
    assert outside.blue() == 0
    assert converted.pixelColor(0, 0).alpha() == 0


def test_fetch_solar_hover_image_uses_metadata_and_cache(tmp_path) -> None:
    response = json.dumps(
        {
            "id": "154849463",
            "date": "2024-03-20 11:59:53",
            "scale": 0.5976089064757952,
            "rsun": 1605.8075,
        }
    ).encode()
    calls: list[str] = []
    image_bytes = _png_bytes()

    def opener(request, timeout):
        del timeout
        url = request.full_url
        calls.append(url)
        return _Response(response if "getClosestImage" in url else image_bytes)

    target = datetime(2024, 3, 20, 12, tzinfo=timezone.utc)
    loaded = fetch_solar_hover_image(target, cache_root=tmp_path, opener=opener)
    assert loaded.image.size().width() == 1024
    assert loaded.time_utc == datetime(2024, 3, 20, 11, 59, 53, tzinfo=timezone.utc)
    assert loaded.image_id == 154849463
    assert 399.0 < loaded.source_radius_px < 401.0
    assert loaded.stale is False

    cached = fetch_solar_hover_image(target, cache_root=tmp_path, opener=opener)
    assert cached.image.size().height() == 1024
    assert calls == [closest_image_url(target), screenshot_url(loaded.time_utc)]


def test_fetch_solar_hover_image_falls_back_to_stale_cache(tmp_path) -> None:
    metadata = json.dumps(
        {
            "id": "154849463",
            "date": "2024-03-20 11:59:53",
            "scale": 0.5976089064757952,
            "rsun": 1605.8075,
        }
    ).encode()
    image_bytes = _png_bytes()

    def first_opener(request, timeout):
        del timeout
        return _Response(metadata if "getClosestImage" in request.full_url else image_bytes)

    fetch_solar_hover_image(
        datetime(2024, 3, 20, 12, tzinfo=timezone.utc),
        cache_root=tmp_path,
        opener=first_opener,
    )

    def failed_opener(request, timeout):
        del request, timeout
        raise OSError("offline")

    stale = fetch_solar_hover_image(
        datetime(2024, 3, 20, 13, tzinfo=timezone.utc),
        cache_root=tmp_path,
        opener=failed_opener,
    )
    assert stale.image_id == 154849463
    assert stale.stale is True
