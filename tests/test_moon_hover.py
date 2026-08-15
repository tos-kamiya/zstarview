from __future__ import annotations

import json
from datetime import datetime, timezone

from PySide6.QtCore import QBuffer, QIODevice
from PySide6.QtGui import QImage

from zstarview.moon_hover import (
    dialamoon_url,
    fetch_moon_hover_image,
    normalize_dialamoon_time,
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


def _jpeg_bytes() -> bytes:
    image = QImage(730, 730, QImage.Format.Format_RGB32)
    image.fill(0)
    buffer = QBuffer()
    buffer.open(QIODevice.OpenModeFlag.WriteOnly)
    assert image.save(buffer, "JPG")
    return bytes(buffer.data())


def test_dialamoon_time_rounds_to_nearest_hour() -> None:
    target = datetime(2026, 8, 15, 6, 29, 59, tzinfo=timezone.utc)
    assert normalize_dialamoon_time(target) == datetime(
        2026, 8, 15, 6, tzinfo=timezone.utc
    )
    assert dialamoon_url(target).endswith("/2026-08-15T06:00")


def test_fetch_moon_hover_image_uses_api_metadata_and_cache(tmp_path) -> None:
    image_url = "https://example.test/moon.5431.jpg"
    response = json.dumps(
        {
            "image": {"url": image_url},
            "time": "2026-08-15T06:00",
            "phase": 8.07,
            "age": 2.516,
            "diameter": 1894.2,
            "posangle": 21.956,
        }
    ).encode()
    calls: list[str] = []
    image_bytes = _jpeg_bytes()

    def opener(request, timeout):
        del timeout
        url = request.full_url
        calls.append(url)
        return _Response(response if url.endswith("T06:00") else image_bytes)

    target = datetime(2026, 8, 15, 6, 29, tzinfo=timezone.utc)
    loaded = fetch_moon_hover_image(target, cache_root=tmp_path, opener=opener)
    assert loaded.image.size().width() == 730
    assert loaded.time_utc.hour == 6
    assert loaded.diameter_arcsec == 1894.2
    assert loaded.posangle_deg == 21.956

    cached = fetch_moon_hover_image(target, cache_root=tmp_path, opener=opener)
    assert cached.image.size().height() == 730
    assert calls == [dialamoon_url(target), image_url]
